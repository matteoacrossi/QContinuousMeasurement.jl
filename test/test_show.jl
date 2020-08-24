using QContinuousMeasurement
using Test

@testset "Pretty printing" begin
    @testset "LocalDephasingModel" begin
        Nj = 2
        modelparams = ModelParameters(Nj=Nj)

        model = LocalDephasingModel(modelparams)
        initial_state = blockdiag_css(Nj)

        println(model)
        println(initial_state)
    end
end