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