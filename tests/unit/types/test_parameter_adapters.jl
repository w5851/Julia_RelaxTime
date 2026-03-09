using Test

if !isdefined(Main, :ParameterTypes)
    Base.include(Main, joinpath(@__DIR__, "..", "..", "..", "src", "types", "ParameterTypes.jl"))
end
if !isdefined(Main, :ParameterAdapters)
    Base.include(Main, joinpath(@__DIR__, "..", "..", "..", "src", "utils", "ParameterAdapters.jl"))
end

using Main.ParameterTypes: QuarkParams, ThermoParams, as_namedtuple
using Main.ParameterAdapters: normalize_quark_input, normalize_thermo_input

@testset "ParameterAdapters" begin
    q_nt = (m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.1))
    t_nt = (T=0.15, Φ=0.5, Φbar=0.45, ξ=0.0)

    q_struct = QuarkParams(q_nt)
    t_struct = ThermoParams(t_nt)

    @test normalize_quark_input(q_nt) === q_nt
    @test normalize_thermo_input(t_nt) === t_nt

    @test normalize_quark_input(q_struct) == as_namedtuple(q_struct)
    @test normalize_thermo_input(t_struct) == as_namedtuple(t_struct)

    @test normalize_quark_input(q_struct) == q_nt
    @test normalize_thermo_input(t_struct) == t_nt
end