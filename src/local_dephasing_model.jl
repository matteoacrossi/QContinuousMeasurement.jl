using BlockDiagonalMatrices
using SparseArrays
using LinearAlgebra
using TimerOutputs

struct LocalDephasingModel <: Model
    params::ModelParameters
    Jx::BlockDiagonal
    Jy::BlockDiagonal
    Jz::BlockDiagonal
    Jx2::BlockDiagonal
    Jy2::BlockDiagonal
    Jz2::BlockDiagonal
    inefficient_measurement::BlockDiagonal
    second_term::SuperOperator
    M0::BlockDiagonal
    dM::BlockDiagonal
    measurement::Array{Eigen}

    function LocalDephasingModel(modelparams::ModelParameters, liouvillianfile::Union{String,Nothing}=nothing)
        Nj = modelparams.Nj
        dt = modelparams.dt
        Gamma = modelparams.Gamma
        kind = modelparams.kind
        ω = modelparams.omega
        η = modelparams.eta

        # Spin operators
        (Jx, Jy, Jz) = map(blockdiagonal, jspin(Nj))

        # TODO: Find better name
        inefficient_measurement = sqrt((1 - η) * dt * Gamma) * Jy

        second_term = dt * (kind / 2) * (isnothing(liouvillianfile) ?
                                            initliouvillian(Nj) :
                                            initliouvillian(Nj, liouvillianfile))

        dropzeros!(second_term)
        second_term = SuperOperator(second_term)

        Jx2 = Jx^2
        Jy2 = Jy^2
        Jz2 = Jz^2
        H = ω * Jz
        dH = Jz

        # Kraus-like operator, trajectory-independent part
        M0 = (I - 1im * H * dt -
              0.25 * dt * kind * Nj * I - # The Id comes from the squares of sigmaz_j
              (Gamma / 2) * Jy2 * dt)

        # Derivative of the Kraus-like operator wrt to ω
        dM = -1im * dH * dt

        measurement = map(x -> eigen(Matrix(x)), Jy.blocks)

        new(modelparams, Jx, Jy, Jz, Jx2, Jy2, Jz2, inefficient_measurement, second_term, M0, dM, measurement)
    end
end

    get_time(m::Model) = get_time(m.params)

function initliouvillian(Nj::Integer)
sys = piqs.Dicke(Nj)
    sys.dephasing = 4.

    liouvillian = tosparse(sys.liouvillian())
    return liouvillian + Nj * I
end

function initliouvillian(Nj::Integer, filename::String)
    liouvillian = sparse_fromfile(filename)
    return liouvillian + Nj * I
end

function measure_current(state::BlockDiagonalState, model::LocalDephasingModel)
    @inline dW() = sqrt(model.params.dt) * randn() # Define the Wiener increment
    # Homodyne current (Eq. 35)
    # dy = 2 sqrt(Gamma * eta) * tr(ρ * Jy) * dt + dW
    mul!(state._tmp1, model.Jy, state.ρ)
    return 2 * sqrt(model.params.Gamma * model.params.eta) * real(tr(state._tmp1)) * model.params.dt + dW()
end

function get_kraus(model::LocalDephasingModel, dy::Real)
    p = model.params
    kraus = similar(model.M0)
    # Kraus operator Eq. (36)
    @inbounds for i in eachindex(kraus.blocks)
        kraus.blocks[i] = model.M0.blocks[i] + sqrt(p.eta * p.Gamma) * model.Jy.blocks[i] * dy +
        p.eta * (p.Gamma / 2) * model.Jy2.blocks[i] * (dy^2 - p.dt)
    end
    return kraus
end

function updatestate!(state::BlockDiagonalState, model::LocalDephasingModel, dy::Real)
    # Get Kraus operator for measured current dy
    @timeit_debug "get_kraus" M = get_kraus(model, dy)

    # Non-allocating code for
    # new_ρ = Mpre * Mpost * ρ + second_term * ρ
        
    @timeit_debug "rho" begin
        mul!(state._tmp1, state.ρ, M')
        mul!(state._new_ρ, M, state._tmp1)

            @timeit_debug "superop" apply_superop!(state._tmp1, model.second_term, state.ρ)

        if model.params.eta < 1.0
            mul!(state._tmp2, model.inefficient_measurement, state.ρ)
            mul!(state._tmp1, state._tmp2, model.inefficient_measurement', 1.0, 1.0)
        end

        # TODO: Replace with broadcasting once implemented
        @inbounds for i in eachindex(state._new_ρ.blocks)
            state._new_ρ.blocks[i] .+= state._tmp1.blocks[i]
        end

        zchop!(state._new_ρ) # Round off elements smaller than 1e-14
        tr_ρ = tr(state._new_ρ)
        # Evolve the unnormalized derivative wrt ω
    end
    # Non-allocating code for
    # τ = (Mpre * (Mpost * τ  +  dMpost * ρ) + dMpre * Mpost * ρ +
    #      second_term * τ )/ tr_ρ;
    @timeit_debug "tau" begin
        mul!(state._tmp1, state.ρ, model.dM')
        mul!(state._tmp1, state.τ, M', 1., 1.)

        @timeit_debug "superop" apply_superop!(state._tmp2, model.second_term, state.τ)
        mul!(state._tmp2, M, state._tmp1, 1., 1.)

        if model.params.eta < 1.0
            mul!(state._tmp1, model.inefficient_measurement, state.τ)
            mul!(state._tmp2, state._tmp1, model.inefficient_measurement', 1.0, 1.0)
        end

        mul!(state._tmp1, state.ρ, M')
        mul!(state._tmp2, model.dM, state._tmp1, 1., 1.)

        # TODO: Use broadcasting when it is implemented
        @inbounds for i in eachindex(state.τ.blocks)
            state.τ.blocks[i] .= state._tmp2.blocks[i] ./ tr_ρ
        end
    end
    zchop!(state.τ) # Round off elements smaller than 1e-14

    tr_τ = tr(state.τ)
    # Now we can renormalize ρ and its derivative wrt ω
    # TODO: Use broadcasting when it is implemented
    @timeit_debug "norms" begin
        @inbounds for i = eachindex(state.ρ.blocks)
            state.ρ.blocks[i] .= state._new_ρ.blocks[i] ./ tr_ρ
        end
        @inbounds for i = eachindex(state.dρ.blocks)
            state.dρ.blocks[i] .= state.τ.blocks[i] .- tr_τ .* state.ρ.blocks[i]
        end
    end
    return tr_ρ, tr_τ
end

expectation_value!(state::State, op::AbstractArray) = real(tr(mul!(state._tmp1, op, density_matrix(state))))

coherentspinstate(model::LocalDephasingModel) = blockdiag_css(model.params.Nj)
