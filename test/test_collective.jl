using QContinuousMeasurement
using Test
using Random
using QContinuousMeasurement
using Test
using Random

@testset "collective" begin
    @testset "Nj = 1" begin
        # We check that Local and Collective dephasing are the same
        # for a single spin.

        params = ModelParameters(Nj=1, kind=1, kcoll=0, Gamma=1)
        model = LocalDephasingModel(params)
        state = blockdiag_css(params.Nj)

        params_coll = ModelParameters(Nj=1, kind=0, kcoll=1, Gamma=1)
        model_coll = CollectiveDephasingModel(params_coll)
        state_coll = fixedj_css(params_coll.Nj)

        seed = 1

        Random.seed!(seed)
        res = simulate_trajectory(model, state)

        Random.seed!(seed)
        res_coll = simulate_trajectory(model_coll, state_coll)

        for k in keys(res)
            @testset "$k" begin
                @test getfield(res, k) ≈ getfield(res_coll, k)
            end
        end
    end

    @testset "Nj = 5" begin
        # We check that <Jx> is the same for local and collective
        # when there is no monitoring

        params = ModelParameters(Nj=5, kind=1, kcoll=0, Gamma=0)
        model = LocalDephasingModel(params)
        state = blockdiag_css(params.Nj)

        params_coll = ModelParameters(Nj=5, kind=0, kcoll=1, Gamma=0)
        model_coll = CollectiveDephasingModel(params_coll)
        state_coll = fixedj_css(params_coll.Nj)

        seed = 1

        Random.seed!(seed)
        res = simulate_trajectory(model, state)

        Random.seed!(seed)
        res_coll = simulate_trajectory(model_coll, state_coll)

        @testset "Jx" begin
            @test getfield(res, :Jx) ≈ getfield(res_coll, :Jx) rtol=10*params.dt
        end
    end
end