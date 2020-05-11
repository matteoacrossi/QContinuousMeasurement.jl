using QContinuousMeasurement
using Test

@testset "Pretty printing" begin
    @testset "CollectiveLocalDephasingModel" begin
        Nj = 2
        modelparams = CollectiveLocalDephasingModelParameters(Nj=Nj)

        model = CollectiveLocalDephasingModel(modelparams)
        initial_state = blockdiag_css(Nj)

        println(model)
        println(initial_state)
    end
end