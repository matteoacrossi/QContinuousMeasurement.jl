using BlockDiagonalMatrices
using SparseArrays
using LinearAlgebra
using TimerOutputs

struct CollectiveDephasingModel <: Model
    params::ModelParameters
    Jx::SparseMatrixCSC
    Jy::SparseMatrixCSC
    Jz::SparseMatrixCSC
    Jx2::SparseMatrixCSC
    Jy2::SparseMatrixCSC
    Jz2::SparseMatrixCSC
    second_terms::Array{SparseMatrixCSC}
    M0::SparseMatrixCSC
    dM::SparseMatrixCSC
    measurement::Eigen


    function CollectiveDephasingModel(modelparams::ModelParameters, liouvillianfile::Union{String, Nothing}=nothing)
        Nj = modelparams.Nj
        dt = modelparams.dt
        kcoll = modelparams.kcoll
        Gamma = modelparams.Gamma
        ω = modelparams.omega
        η = modelparams.eta

        firstblock = block_sizes(Nj)[1]

        # Spin operators (only the first block)
        (Jx, Jy, Jz) = map(x -> sparse(x[1:firstblock, 1:firstblock]), jspin(Nj))

        # TODO: Find better name
        second_terms = [sqrt((1 - η) * dt * Gamma) * Jy, sqrt(dt * kcoll) * Jz]

        Jx2 = Jx^2
        Jy2 = Jy^2
        Jz2 = Jz^2

        H = ω * Jz
        dH = Jz

        # Kraus-like operator, trajectory-independent part
        M0 = I - 1im * H * dt - (Gamma/2) * Jy2 * dt - (kcoll/2) * Jz2 * dt

        # Derivative of the Kraus-like operator wrt to ω
        dM = -1im * dH * dt

        measurement = eigen(Matrix(Jy))

        new(modelparams, Jx, Jy, Jz, Jx2, Jy2, Jz2, second_terms, M0, dM, measurement)
    end
end

function measure_current(state::FixedjState, model::CollectiveDephasingModel)
    @inline dW() = sqrt(model.params.dt) * randn() # Define the Wiener increment
    # Homodyne current (Eq. 35)
    # dy = 2 sqrt(Gamma * eta) * tr(ρ * Jy) * dt + dW
    mul!(state._tmp1, model.Jy, state.ρ)
    return 2 * sqrt(model.params.Gamma * model.params.eta) * real(tr(state._tmp1)) * model.params.dt + dW()
end

function get_kraus(model::CollectiveDephasingModel, dy::Real)
    p = model.params
    kraus = similar(model.M0)

    # Kraus operator Eq. (36)
    kraus = (model.M0 + sqrt(p.eta * p.Gamma) * model.Jy * dy +
            p.eta * (p.Gamma / 2) * model.Jy2 * (dy^2 - p.dt))

    return kraus
end

function updatestate!(state::FixedjState, model::CollectiveDephasingModel, dy::Real)
    # Get Kraus operator for measured current dy
    @timeit_debug "get_kraus" M = get_kraus(model, dy)

    # Non-allocating code for
    # new_ρ = Mpre * Mpost * ρ + second_term * ρ

    @timeit_debug "rho" begin
        mul!(state._tmp1, state.ρ, M')
        mul!(state._new_ρ, M, state._tmp1)

        # Apply the second term of Eq. (34)
        for second_term in model.second_terms
            mul!(state._tmp1, state.ρ, second_term')
            mul!(state._new_ρ, second_term, state._tmp1, 1., 1.)
        end

        zchop!(state._new_ρ) # Round off elements smaller than 1e-14
        tr_ρ = tr(state._new_ρ)
        # Evolve the unnormalized derivative wrt ω
    end

    # Non-allocating code for
    # τ = (M * (τ * M' + ρ * dM') + dM * ρ * M' +
    #     ∑ second_term[i] * τ * second_term[i]') / tr_ρ;
    @timeit_debug "tau" begin

        fill!(state._tmp2, 0.0)
        # ∑ second_term[i] * τ * second_term[i]')
        for second_term in model.second_terms
            mul!(state._tmp1, second_term, state.τ)
            mul!(state._tmp2, state._tmp1, second_term', 1., 1.)
        end

        # (τ * M' + ρ * dM')
        mul!(state._tmp1, state.ρ, model.dM')
        mul!(state._tmp1, state.τ, M', 1., 1.)

        copy!(state.τ, state._tmp2)

        #  (M * (τ * M' + ρ * dM')
        mul!(state.τ, M, state._tmp1, 1., 1.)

        #  ... + dM * ρ * M'
        mul!(state._tmp1, state.ρ, M')
        mul!(state.τ, model.dM, state._tmp1, 1., 1.)

        state.τ ./=  tr_ρ
    end

    zchop!(state.τ) # Round off elements smaller than 1e-14

    tr_τ = tr(state.τ)
    # Now we can renormalize ρ and its derivative wrt ω
    # TODO: Use broadcasting when it is implemented
    @timeit_debug "norms" begin
        state.ρ  .= state._new_ρ ./ tr_ρ
        state.dρ .= state.τ .- tr_τ .* state.ρ
    end
    return tr_ρ, tr_τ
end

coherentspinstate(model::CollectiveDephasingModel) = fixedj_css(model.params.Nj)