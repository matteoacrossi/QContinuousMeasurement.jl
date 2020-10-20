using ZChop # For chopping small imaginary parts in ρ
using Distributed
using TimerOutputs
using ProgressMeter

function simulate_trajectory(model::Tm,
                             initial_state::Ts,
                             file_channel::Union{Channel,RemoteChannel,Nothing}=nothing,
                             progress_channel::Union{Channel,RemoteChannel,Nothing}=nothing;
                             callbacks=[]) where {Tm <: Model, Ts <:State}

    state = Ts(initial_state) # Copy the initial state

    jx = Array{Float64}(undef, length(get_time(model)))

    # Output variables
    jx = similar(jx)
    jy = similar(jx)
    jz = similar(jx)

    jx2 = similar(jx)
    jy2 = similar(jx)
    jz2 = similar(jx)

    FisherT = similar(jx)
    QFisherT = similar(jx)

    FisherStrong = similar(jx)


    jto = 1 # Counter for the output

    for jt = 1 : model.params.Ntime
        dy = measure_current(state, model)
        tr_ρ, tr_τ = updatestate!(state, model, dy)

        if jt % model.params._outsteps == 0

            jx[jto] = expectation_value!(state, model.Jx)
            jy[jto] = expectation_value!(state, model.Jy)
            jz[jto] = expectation_value!(state, model.Jz)

            jx2[jto] = expectation_value!(state, model.Jx2)
            jy2[jto] = expectation_value!(state, model.Jy2)
            jz2[jto] = expectation_value!(state, model.Jz2)

            # We evaluate the classical FI for the continuous measurement
            FisherT[jto] = real(tr_τ^2)

            # We evaluate the QFI for a final strong measurement done at time t
            QFisherT[jto] = QFI!(state)

            # Fisher of the strong final measurement
            FisherStrong[jto] = FI(state, model.measurement)

            # apply any callbacks
            isnothing(callbacks) || for cb in callbacks
                cb(state, model)
            end

            jto += 1
            isnothing(progress_channel) || put!(progress_channel, true)
        end
    end

    Δjx2 = jx2 - jx.^2
    Δjy2 = jy2 - jy.^2
    Δjz2 = jz2 - jz.^2

    xi2x = squeezing_param(model.params.Nj, Δjx2, jy, jz)
    xi2y = squeezing_param(model.params.Nj, Δjy2, jx, jz)
    xi2z = squeezing_param(model.params.Nj, Δjz2, jx, jy)

    result = (FI=FisherT, QFI=QFisherT,
                    Jx=jx, Jy=jy, Jz=jz,
                    Δjx2=Δjx2, Δjy2=Δjy2, Δjz2=Δjz2,
                    xi2x=xi2x, xi2y=xi2y, xi2z=xi2z, FIstrong=FisherStrong)
    isnothing(file_channel) || put!(file_channel, result)

    GC.gc()
    return result
end


function Eff_QFI_HD_Dicke(Nj::Int64, # Number of spins
    Ntraj::Int64,                    # Number of trajectories
    Tfinal::Real,                    # Final time
    dt::Real;                        # Time step
    κ::Real = 1.,                    # Independent noise strength
    κcoll::Real = 1.,                # Collective noise strength
    ω::Real = 0.0,                   # Frequency of the Hamiltonian
    η::Real = 1.,                    # Measurement efficiency
    outpoints = 0,                   # Number of output points
    to = TimerOutput(), file_channel=nothing)

    @info "Eff_QFI_HD_Dicke starting"
    @info "Parameters" Nj Ntraj Tfinal dt κ κcoll ω η outpoints

    modelparams = ModelParameters(Nj=Nj, kind=κ, Gamma=κcoll, kcoll=0., omega=ω, eta=η, dt=dt, Tfinal=Tfinal, outpoints=outpoints)

    model = LocalDephasingModel(modelparams)
    initial_state = blockdiag_css(Nj)

    # Run evolution for each trajectory, and build up the average
    # for FI and final strong measurement QFI
    @timeit_debug to "trajectories" begin
    #result = @showprogress 1 "Computing..." @distributed (+) for ktraj = 1 : Ntraj
    result = @distributed (+) for ktraj = 1 : Ntraj
        result = simulate_trajectory(model, initial_state)
        hcat(result.FI, result.QFI,
             result.Jx, result.Jy, result.Jz,
             result.Δjx2, result.Δjy2, result.Δjz2,
             result.xi2x, result.xi2y, result.xi2z, result.FIstrong)
    end
    end

    # Evaluate averages
    jx = result[:, 3] / Ntraj
    jy = result[:, 4] / Ntraj
    jz = result[:, 5] / Ntraj

    Δjx2 = result[:, 6] / Ntraj
    Δjy2 = result[:, 7] / Ntraj
    Δjz2 = result[:, 8] / Ntraj

    xi2x = result[:, 9] / Ntraj
    xi2y = result[:, 10] / Ntraj
    xi2z = result[:, 11] / Ntraj

    @info "Time details \n$to"
    return (t=get_time(model),
            FI=result[:,1] / Ntraj,
            QFI=result[:,2] / Ntraj,
            jx=jx, jy=jy, jz=jz,
            Δjx=Δjx2, Δjy=Δjy2, Δjz=Δjz2,
            xi2x=xi2x, xi2y=xi2y, xi2z=xi2z, FIstrong=result[:, 12] / Ntraj)
end