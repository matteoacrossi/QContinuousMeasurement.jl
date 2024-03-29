using ZChop
function squeezing_param(N, ΔJ1, J2m, J3m)
    """
        ξ2 = squeezing_param(N, ΔJ1, J2m, J3m)

    Returns the squeezing parameter defined, e.g., as
    the inverse of Eq. (1) in Phys. Rev. A 65, 061801 (2002).

    The state is squeezed if greater than one.
    """
    return (J2m.^2 + J3m.^2) ./ (N * ΔJ1)
end


function squeezing_param(rho::QContinuousMeasurement.FixedjState)
    """
        ξ2 = squeezing_param(rho)

    Returns the squeezing parameter defined, e.g., as
    the inverse of Eq. (1) in Phys. Rev. A 65, 061801 (2002).

    The state is squeezed if greater than one.
    """
    Nj = nspins(rho)
    firstblock = block_sizes(Nj)[1]
    (Jx, Jy, Jz) =  map(x -> sparse(x[1:firstblock, 1:firstblock]), jspin(Nj))

    Jy2 = Jy^2
    Jzm = expectation_value!(rho, Jz)
    Jxm =  expectation_value!(rho, Jx)
    Jym =  expectation_value!(rho, Jy)
    Jy2m = expectation_value!(rho, Jy2)
    return (Jzm.^2 + Jxm.^2) ./ (Nj * (Jy2m - Jym.^2))
end

function squeezing_param(rho::QContinuousMeasurement.State)
    """
        ξ2 = squeezing_param(rho)

    Returns the squeezing parameter defined, e.g., as
    the inverse of Eq. (1) in Phys. Rev. A 65, 061801 (2002).

    The state is squeezed if greater than one.
    """
    Nj = nspins(rho)
    (Jx, Jy, Jz) = jspin(Nj)
    Jy2 = Jy^2
    Jzm = expectation_value!(rho, Jz)
    Jxm =  expectation_value!(rho, Jx)
    Jym =  expectation_value!(rho, Jy)
    Jy2m = expectation_value!(rho, Jy2)
    return (Jzm.^2 + Jxm.^2) ./ (Nj * (Jy2m - Jym.^2))
end

"""
    density(M)

Returns the density of a sparse matrix
"""
function density(s::SparseMatrixCSC)
    return length(s.nzval) / (s.n * s.m)
end

function Unconditional_QFI_Dicke(Nj::Int64, Tfinal::Real, dt::Real;
    κ::Real=1.,                    # Independent noise strength
    κcoll::Real=1.,                # Collective noise strength
    ω::Real=0.0                   # Frequency of the Hamiltonian
    )
    return Eff_QFI_HD_Dicke(Nj, 1, Tfinal, dt; κ=κ, κcoll=κcoll, ω=ω, η=0.0)
end

# Specialize zchop! to BlockDiagonal
function ZChop.zchop!(A::BlockDiagonal)
    for b in blocks(A)
        zchop!(b)
    end
end