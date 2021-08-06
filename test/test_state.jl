using QContinuousMeasurement
using Test

@testset "css" begin

    @testset "LocalDephasing" begin
        params = ModelParameters(Nj=2, kind=0, kcoll=0)
        model = LocalDephasingModel(params)
        state = blockdiag_css(params.Nj)

        @test state isa BlockDiagonalState
        @test coherentspinstate(model).ρ == state.ρ
    end

    @testset "CollectiveDephasing" begin
        params = ModelParameters(Nj=2, kind=0, kcoll=0)
        model = CollectiveDephasingModel(params)
        state = fixedj_css(params.Nj)

        @test state isa FixedjState
        @test coherentspinstate(model).ρ == state.ρ
    end
end

@testset "SpinSqueezed" begin
    params = ModelParameters(Nj=2, kind=0, kcoll=0)
    model = LocalDephasingModel(params)
    state = blockdiag_sss(params.Nj, 0.0)
    css_state = blockdiag_css(params.Nj)

    @test state.ρ ≈ css_state.ρ
    @test squeezing_param(state) ≈ 1.0
    @test state isa BlockDiagonalState

    sq_state = blockdiag_sss(params.Nj, 0.1)
    @test squeezing_param(sq_state) > 1.0
end