using QContinuousMeasurement
using Test
using Random

@testset "kind=0" begin
    params = QContinuousMeasurement.LocalDephasingModelParameters(Nj = 2, kind = 0)
    model = QContinuousMeasurement.LocalDephasingModel(params)
    state = QContinuousMeasurement.blockdiag_css(params.Nj)

    params2 = QContinuousMeasurement.CollectiveDephasingModelParameters(Nj = 2)
    model2 = QContinuousMeasurement.CollectiveDephasingModel(params2)
    state2 = QContinuousMeasurement.fixedj_css(params2.Nj)

    seed = 1

    Random.seed!(seed)
    res = simulate_trajectory(model, state)

    Random.seed!(seed)
    res2 = simulate_trajectory(model2, state2)

    for k in keys(res)
        @test getfield(res, k) ≈ getfield(res2, k)
    end
end