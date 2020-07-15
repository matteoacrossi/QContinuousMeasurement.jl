using QContinuousMeasurement
using HDF5
using DelimitedFiles
using TimerOutputs
using LinearAlgebra

function random_rholike(N)
    blocksizes = QContinuousMeasurement.block_sizes(N)
    blocks = map(x -> x + x', map(x-> rand(ComplexF64, x, x), blocksizes))
    r = BlockDiagonal(blocks)
    return r / tr(r)
end

outsteps = 1
Tfinal = 1.0
dt = 0.01
outpoints = 10
Ntraj = 3

if outpoints > 0
    try
        outsteps = Int(round(Tfinal / dt / outpoints, digits=3))
    catch InexactError
        @warn "The requested $outpoints output points does not divide
        the total time steps. Using the full time output."
    end
end

outsteps = Int(round(Tfinal / dt / outpoints, digits=3))

modelparams = ModelParameters(Nj=40, kind=1.0, kcoll=1.0, omega=0.01, Tfinal=1.0, dt=dt)

Ntime = Int(floor(Tfinal/dt)) # Number of timesteps
t = (1 : Ntime) * dt
t = t[outsteps:outsteps:end]

# Precompile
model_p = InitializeModel(modelparams)
initial_state_p = blockdiag_css(modelparams.Nj)
dy = measure_current(initial_state_p, model_p)
#updatekraus!(model_p, dy)
tr_ρ, tr_τ = updatestate!(initial_state_p, model_p, dy)
QFI(initial_state_p)
expectation_value!(initial_state_p, model_p.Jx)

TimerOutputs.enable_debug_timings(QContinuousMeasurement)
reset_timer!()
@timeit "Model creation" model = InitializeModel(modelparams)
@timeit "State creation" initial_state = blockdiag_css(modelparams.Nj)

@info "Starting trajectories simulation..."
for ktraj = 1 : Ntraj
    @timeit "State copy" state = State(copy(initial_state.ρ))

    # Output variables
    jx = similar(t)
    jy = similar(t)
    jz = similar(t)

    jx2 = similar(t)
    jy2 = similar(t)
    jz2 = similar(t)

    FisherT = similar(t)
    QFisherT = similar(t)

    jto = 1 # Counter for the output
    for jt = 1 : Ntime
        @timeit "Current" dy = measure_current(state, model)
        #@timeit to "Kraus" updatekraus!(model, dy)
        @timeit "Dynamics" tr_ρ, tr_τ = updatestate!(state, model, dy)

        if jt % outsteps == 0
        @timeit "Output" begin
            jx[jto] = expectation_value!(state, model.Jx)
            jy[jto] = expectation_value!(state, model.Jy)
            jz[jto] = expectation_value!(state, model.Jz)

            jx2[jto] = expectation_value!(state, model.Jx2)
            jy2[jto] = expectation_value!(state, model.Jy2)
            jz2[jto] = expectation_value!(state, model.Jz2)

            # We evaluate the classical FI for the continuous measurement
            FisherT[jto] = real(tr_τ^2)

            # We evaluate the QFI for a final strong measurement done at time t
            @timeit "QFI" QFisherT[jto] = QFI(state)

            jto += 1
            end
        end
    end
end
print_timer()