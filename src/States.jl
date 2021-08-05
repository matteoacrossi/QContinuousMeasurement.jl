#= 
Definition of a state =#
using LinearAlgebra

abstract type State end

struct BlockDiagonalState <: State
    ρ::BlockDiagonal
    dρ::BlockDiagonal
    τ::BlockDiagonal
    # Internal fields
    _tmp1::BlockDiagonal
    _tmp2::BlockDiagonal
    _new_ρ::BlockDiagonal

    BlockDiagonalState(ρ::BlockDiagonal) = new(ρ, zero(ρ), zero(ρ), similar(ρ), similar(ρ), similar(ρ))
end

density_matrix(state::BlockDiagonalState) = state.ρ

nspins(state::BlockDiagonalState) = nspins(size(density_matrix(state), 1))

# Copy constructor
BlockDiagonalState(state::BlockDiagonalState) = BlockDiagonalState(copy(density_matrix(state)))

function blockdiag_css(Nj, dense=true)
    ρ0 = blockdiagonal(css(Nj), dense=dense)
    BlockDiagonalState(ρ0)
end

function blockdiag_sss(Nj, mu::Real, dense=true)
    ρ0 = blockdiagonal(spin_squeezed_state(Nj, mu), dense=dense)
    BlockDiagonalState(ρ0)
end

struct FixedjState <: State
    ρ::Matrix
    dρ::Matrix
    τ::Matrix
    # Internal fields
    _tmp1::Matrix
    _tmp2::Matrix
    _new_ρ::Matrix

    FixedjState(ρ::Matrix) = new(ρ, zero(ρ), zero(ρ), similar(ρ), similar(ρ), similar(ρ))
end

density_matrix(state::FixedjState) = state.ρ
FixedjState(state::FixedjState) = FixedjState(copy(density_matrix(state)))

function fixedj_css(Nj::Int64)
    state = css(Nj)
    firstblock = block_sizes(Nj)[1]
    return FixedjState(Matrix(state[1:firstblock, 1:firstblock]))
end

function spin_squeezed_state(Nj::Int64, mu::Real, nu::Real)
    css_state = css(Nj)
    (Jx, Jy, Jz) = jspin(Nj)

    # We only need the first block, with size Nj + 1
    css_block = css_state[1:Nj + 1, 1:Nj + 1]

    Jx = Jx[1:Nj + 1, 1:Nj + 1]
    Jz = Jz[1:Nj + 1, 1:Nj + 1]

    U = exp(-1im * mu / 2 * Array(Jz^2))
    R = exp(-1im * nu * Array(Jx))

    css_block = R * U * css_block * U' * R'
    css_state[1:Nj + 1, 1:Nj + 1] = css_block
    
    return css_state
end

function spin_squeezed_state(Nj::Int64, mu::Real)
    # Applies the optimal rotation (Kitagawa, Ueda, above Eq. (4))
    function delta(J, μ)
        A = 1 - cos(μ)^(2 * J - 2)
        B = 4 * sin(μ / 2) * cos(μ / 2)^(2 * J - 2)
        if A ≈ 0
            return pi / 4
        else
            return 0.5 * atan(B / A)
        end

    end
    δ = delta(Nj / 2, mu)
    return spin_squeezed_state(Nj, mu, pi / 2 - δ)
end

# Older functions

"""
    coherent_state(n)

Return the state ``\\ket{+}^{\\otimes n}`` in the computational basis
"""
function coherent_state(n::Int)
    @assert n > 0 "n must be a positive integer"
    spinup = Vector{Complex{Float64}}([1., 0.])
    spindown = Vector{Complex{Float64}}([0., 1.])

    if n == 1
        return (spinup + spindown) / sqrt(2.)
    end
    return kron([(spinup + spindown) / sqrt(2.) for i in 1:n]...)
end

