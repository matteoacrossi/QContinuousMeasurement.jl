using QContinuousMeasurement
using Test
using Random

@testset "kind=0" begin
    params = ModelParameters(Nj=2, kind=0, kcoll=0)
    model = LocalDephasingModel(params)
    state = blockdiag_css(params.Nj)

    model2 = CollectiveDephasingModel(params)
    state2 = fixedj_css(params.Nj)

    seed = 1

    Random.seed!(seed)
    res = simulate_trajectory(model, state)

    Random.seed!(seed)
    res2 = simulate_trajectory(model2, state2)

    for k in keys(res)
        @test getfield(res, k) ≈ getfield(res2, k)
    end
end