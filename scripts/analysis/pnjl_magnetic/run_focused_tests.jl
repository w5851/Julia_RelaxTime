using Test
const ROOT=normpath(joinpath(@__DIR__,"..","..",".."))
include(joinpath(ROOT,"src","models","Models.jl"))
Models.pnjl_module()
files = isempty(ARGS) ? [
    "tests/unit/pnjl/test_magnetic_integrals.jl",
    "tests/unit/pnjl/test_magnetic_thermodynamics.jl",
    "tests/unit/pnjl/test_pnjl_magnetic_model.jl",
    "tests/unit/models/test_magnetic_scan.jl",
    "tests/regression/pnjl/test_magnetic_fixedpoint_regression.jl",
    "tests/validation/pnjl/reference/test_magnetic_pnjl_mag_fixedpoint_kernel.jl",
    "tests/validation/pnjl/reference/test_magnetic_pnjl_mag_equilibrium_anchors.jl",
] : ARGS
@testset "magnetic correctness focused tests" begin
    for file in files
        println("TEST ",file);flush(stdout)
        m=Module(gensym(:MagneticTest))
        Core.eval(m,:(const Models=Main.Models))
        Core.eval(m,:(const Constants_PNJL=Main.Constants_PNJL))
        Core.eval(m,:(include(path)=Base.include(@__MODULE__,path)))
        Base.include(m,joinpath(ROOT,file))
    end
end
