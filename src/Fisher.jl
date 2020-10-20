using ZChop

"""
    QFI(ρ, dρ [, abstol])

Numerically evaluate the quantum Fisher information for the matrix ρ given its derivative dρ wrt the parameter

This function is the implementation of Eq. (13) in Paris, Int. J. Quantum Inform. 7, 125 (2009).

# Arguments
    * `ρ`:  Density matrix
    * `dhro`: Derivative wrt the parameter to be estimated
    * `abstol = 1e-5`: tolerance in the denominator of the formula
"""
function QFI(ρ, dρ; abstol = 1e-5)

    # Get the eigenvalues and eigenvectors of the density matrix
    # We enforce its Hermiticity so that the algorithm is more efficient and returns real values

    eigval, eigvec = eigen(Hermitian(Matrix(zchop(ρ, 1e-10))))

    dim = length(eigval)
    tmpvec = Array{eltype(ρ)}(undef, size(ρ, 1))
    res = 0.
    tmp = 0.
    for m = 1:dim
        mul!(tmpvec, dρ, view(eigvec, :, m))
        for n = 1:dim
            tmp = eigval[n] + eigval[m]
            if tmp > abstol
                res += 2 * (1. / tmp) * abs2(dot(view(eigvec, :, n), tmpvec))
            end
        end
    end
    return res
end

"""
    QFI!(ρ, dρ [, abstol])

Numerically evaluate the quantum Fisher information for the matrix ρ given its derivative dρ wrt the parameter.
Can overwrite ρ.

This function is the implementation of Eq. (13) in Paris, Int. J. Quantum Inform. 7, 125 (2009).

# Arguments
    * `ρ`:  Density matrix (gets overwritten)
    * `dhro`: Derivative wrt the parameter to be estimated
    * `abstol = 1e-5`: tolerance in the denominator of the formula
"""
function QFI!(ρ, dρ; abstol = 1e-5)
    # Get the eigenvalues and eigenvectors of the density matrix
    # We enforce its Hermiticity so that the algorithm is more efficient and returns real values
    eigval, eigvec = eigen!(Hermitian(ρ))
    dim = length(eigval)
    tmpvec = Array{eltype(ρ)}(undef, size(ρ, 1))
    res = 0.
    tmp = 0.
    for m = 1:dim
        mul!(tmpvec, dρ, view(eigvec, :, m))
        for n = 1:dim
            tmp = eigval[n] + eigval[m]
            if tmp > abstol
                res += 2 * (1. / tmp) * abs2(dot(view(eigvec, :, n), tmpvec))
            end
        end
    end
    return res
end

function QFI(ρ::BlockDiagonal, dρ::BlockDiagonal; abstol = 1e-5)
    qfi = 0.
    for i in 1:nblocks(ρ)
        qfi += QFI(ρ.blocks[i], dρ.blocks[i])
    end
    return qfi
end

function QFI!(ρ::BlockDiagonal, dρ::BlockDiagonal; abstol = 1e-5)
    qfi = 0.
    for i in 1:nblocks(ρ)
        qfi += QFI!(ρ.blocks[i], dρ.blocks[i])
    end
    return qfi
end

"""
    QFI(state::State, [, abstol])

Numerically evaluate the quantum Fisher information for the State state.

This function is the implementation of Eq. (13) in Paris, Int. J. Quantum Inform. 7, 125 (2009).

# Arguments
    * `state`: A State object
    * `abstol = 1e-5`: tolerance in the denominator of the formula
"""
function QFI(state::State, abstol = 1e-5)
    QFI(state.ρ, state.dρ)
end

"""
    QFI!(state::State, [, abstol])

Numerically evaluate the quantum Fisher information for the State state,
overwriting one of its temporary variables.

This function is the implementation of Eq. (13) in Paris, Int. J. Quantum Inform. 7, 125 (2009).

# Arguments
    * `state`: A State object
    * `abstol = 1e-5`: tolerance in the denominator of the formula
"""
function QFI!(state::State, abstol = 1e-5)
    copy!(state._tmp1, state.ρ)
    QFI!(state._tmp1, state.dρ)
end

"""
    FI(ρ, dρ, op, [, tol])

Evaluates the Fisher information for the measurement operator op, for density matrix ρ, with
derivative dρ with respect to the parameter.

tol gives the tolerance to treat the denominator as 0 (and avoid NaNs). Default is `tol=1e-10`.
"""
function FI(ρ::BlockDiagonal, dρ::BlockDiagonal, op::BlockDiagonal)
    ed = eigen(op)
    return FI(ρ, dρ, ed)
end

function FI(ρ::Matrix, dρ::Matrix, measurement::Eigen; tol=1e-10)
    # We use Eq. (4) of Paris, Int. J. Quantum Inf. 7, 125 (2009)
    # We need the projectors onto the subspaces of the operator op
    fi = 0.0
    for i in 1:size(measurement.vectors, 1)
        evec = measurement.vectors[:, i]
        denominator = real(evec' * ρ * evec)

        # Empty blocks would cause a division by 0
        # We use a tolerance
        if denominator > tol
            fi += real(evec' * dρ * evec)^2 / denominator
        end
    end
    return fi
end

function FI(ρ::BlockDiagonal, dρ::BlockDiagonal, measurement::Array{Eigen}; tol=1e-10)
    # We use Eq. (4) of Paris, Int. J. Quantum Inf. 7, 125 (2009)
    # We need the projectors onto the subspaces of the operator op
    fi = 0.0
    for i in 1:nblocks(ρ)
        fi += FI(ρ.blocks[i], dρ.blocks[i], measurement[i], tol=tol)
    end
    return fi
end

FI(state::State, measurement::Eigen) = FI(state.ρ, state.dρ, measurement)
FI(state::State, measurement::Array{Eigen}) = FI(state.ρ, state.dρ, measurement)