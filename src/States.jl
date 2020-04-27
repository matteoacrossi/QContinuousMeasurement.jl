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