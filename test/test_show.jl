using QContinuousMeasurement
using Test

@testset "Pretty printing" begin
Nj = 2
modelparams = ModelParameters(Nj=Nj)

model = InitializeModel(modelparams)
initial_state = blockdiag_css(Nj)

println(model)
println(initial_state)

end