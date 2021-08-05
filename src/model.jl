using BlockDiagonalMatrices
using SparseArrays
using LinearAlgebra
using TimerOutputs

abstract type Model end

struct ModelParameters
    Nj::Integer
    kind::Real
    kcoll::Real
    Gamma::Real
    omega::Real
    eta::Real
    dt::Real
    Tfinal::Real
    Ntime::Integer
    outpoints::Integer
    mu::Real
    _outsteps::Integer

    function ModelParameters(; Nj::Integer=1,
                               kind::Real=1.0,
                               kcoll::Real=1.0,
                               Gamma::Real=1.0,
                               omega::Real=1.0,
                               eta::Real=1.0,
                               dt::Real=0.0001,
                               Tfinal::Real=1.0,
                               outpoints::Integer=0
                               mu::Real=0.0)

        Ntime = Int(floor(Tfinal / dt)) # Number of timesteps

        outsteps = 1
        if outpoints > 0

            try
                outsteps = Int(round(Tfinal / dt / outpoints, digits=3))
            catch InexactError
                @warn "The requested $outpoints output points does not divide
                the total time steps. Using the full time output."
            end
        end
        # outpoints = Ntime

        new(Nj, kind, kcoll, Gamma, omega, eta, dt, Tfinal, Ntime, outpoints, outsteps, mu)
    end
end

function get_time(mp::ModelParameters)
    t = (1:mp.Ntime) * mp.dt
    t[mp._outsteps:mp._outsteps:end]
end

struct KrausOperator
    M0::BlockDiagonal
    M::BlockDiagonal
end
