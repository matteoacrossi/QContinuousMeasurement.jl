"""
Fisher information for magnetometry with continuously monitored spin systems, with independent Markovian noise acting on each spin
"""
module QContinuousMeasurement
    using PyCall
    const qutip = PyNULL()
    const piqs = PyNULL()
    const sp = PyNULL()
    const pystuff = PyNULL()

    function __init__()
        # The commands below import the modules, and make sure that they Arguments
        # installed using Conda.jl
        copy!(qutip, pyimport_conda("matplotlib", "matplotlib"))
        copy!(qutip, pyimport_conda("qutip", "qutip", "conda-forge"))
        copy!(piqs, pyimport_conda("qutip.piqs", "qutip", "conda-forge"))
        copy!(sp, pyimport_conda("scipy.sparse", "scipy"))

        pushfirst!(PyVector(pyimport("sys")."path"), @__DIR__)

        copy!(pystuff, pyimport("pystuff"))
    end

    include("piqs.jl")

    using SparseArrays
    using LinearAlgebra

    export Molmer_QFI_GHZ, Molmer_qfi_transverse, uncond_QFI_transverse
    export QFI, QFI!, squeezing_param
    export Unconditional_QFI, Unconditional_QFI_Dicke
    export Eff_QFI_HD
    export Eff_QFI_HD_Dicke, Eff_QFI_HD_Dicke_0
    export simulate_trajectory
    export liouvillian, jspin, css
    export Model, ModelParameters, State
    export CollectiveDephasingModel, CollectiveDephasingModelParameters, FixedjState
    export CollectiveLocalDephasingModel, CollectiveLocalDephasingModelParameters, BlockDiagonalState
    export updatestate!, measure_current, expectation_value!
    export get_time
    export InitializeModel, blockdiag_css
    export FileWriter

    include("NoiseOperators.jl")
    include("States.jl")
    include("model.jl")
    include("Fisher.jl")
    include("Molmer_QFI.jl")
    include("Uncond_qfi_transverse.jl")

    include("filewriter.jl")
    include("Eff_QFI_HD_Dicke.jl")
    include("Eff_QFI_HD.jl")

    include("utils.jl")
    include("show.jl") # Pretty formatting

end
