#=
Definition of a state
=#

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

# Copy constructor
BlockDiagonalState(state::BlockDiagonalState) = BlockDiagonalState(copy(density_matrix(state)))

function blockdiag_css(Nj, dense=true)
    ρ0 = blockdiagonal(css(Nj), dense=dense)
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

# Older functions

"""
    coherent_state(n)

Return the state ``\\ket{+}^{\\otimes n}`` in the computational basis
"""
function coherent_state(n::Int)
    @assert n > 0 "n must be a positive integer"
    spinup = Vector{Complex{Float64}}([1., 0.])
    spindown = Vector{Complex{Float64}}([0., 1.])

    if n==1
        return (spinup + spindown) / sqrt(2.)
    end
    return kron([(spinup + spindown) / sqrt(2.) for i in 1:n]...)
end