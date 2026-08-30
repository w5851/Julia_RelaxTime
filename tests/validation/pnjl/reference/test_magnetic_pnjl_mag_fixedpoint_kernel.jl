using Test
using StaticArrays

if !isdefined(Main, :Models)
    include(joinpath(@__DIR__, "..", "..", "..", "..", "src", "models", "Models.jl"))
end
if !isdefined(Main, :validation_targets_path)
    include(joinpath(@__DIR__, "..", "..", "common", "data_paths.jl"))
end

const MAGNETIC_KERNEL_TARGET_PATH = validation_targets_path(
    "pnjl", "reference", "pnjl_magnetic_external_kernel_targets_v1.csv",
)

function _load_magnetic_kernel_target(path::String)
    isfile(path) || error("magnetic kernel target CSV not found: $path")
    lines = readlines(path)
    length(lines) == 2 || error("expected one magnetic kernel target row, got $(length(lines) - 1)")
    header = split(lines[1], ',')
    values = split(lines[2], ',')
    length(header) == length(values) || error("invalid magnetic kernel target row")
    return Dict(header[i] => values[i] for i in eachindex(header))
end

@inline _magnetic_kernel_f(row, key) = parse(Float64, row[key])
@inline _magnetic_kernel_i(row, key) = parse(Int, row[key])

@testset "pnjl_mag external MFIR fixed-point kernel" begin
    row = _load_magnetic_kernel_target(MAGNETIC_KERNEL_TARGET_PATH)
    @test row["source_commit"] == "e1fc81d3c3c9d220c49972e54307b66a604cb9db"
    @test row["acceptance_status"] == "accepted_kernel_only"
    ext_hc = _magnetic_kernel_f(row, "external_hc_MeV_fm")
    current_hc = Main.Constants_PNJL.ħc_MeV_fm

    # Match the external source convention without changing the production
    # default: external hc=197.33 is used for input conversion, while the
    # IMC Lambda_QCD value is mapped into Julia's current internal fm units.
    Λ_MeV = 602.3
    Λ_inv_fm = Λ_MeV / ext_hc
    params = Models.PNJLCore.PNJLParams(
        profile="pnjl_mag_external_kernel",
        physics_profile="pnjl_mag_external_kernel",
        hbarc_MeV_fm=ext_hc,
        alpha_em=1 / 137.035999084,
        N_color=3,
        N_flavor=3,
        rho0_fm3=0.16,
        Λ_inv_fm=Λ_inv_fm,
        m_ud0_inv_fm=5.5 / ext_hc,
        m_s0_inv_fm=140.7 / ext_hc,
        G_fm2=1.835 / Λ_inv_fm^2,
        K_fm5=12.36 / Λ_inv_fm^5,
        thermal_p_max_inv_fm=40.0,
        T0_inv_fm=210.0 / ext_hc,
        polyakov_scheme=:log,
        a0=3.51,
        a1=-2.47,
        a2=15.2,
        a3=0.0,
        b3=-1.75,
        b4=7.555,
    )
    imc = Models.MagneticThermodynamics.MagneticIMCParams(
        a=0.0108805,
        b=-1.0133e-4,
        c=0.02228,
        d=1.84558e-4,
        Λ_QCD_MeV=300.0 * current_hc / ext_hc,
    )
    eB_fm2 = _magnetic_kernel_f(row, "eB_GeV2") * (1000.0 / ext_hc)^2
    conf = Models.MagneticThermodynamics.MagneticConfig(
        eB_fm2=eB_fm2,
        # External n_num is a count including n=0; Julia n_max is inclusive.
        n_max=_magnetic_kernel_i(row, "landau_levels") - 1,
        p_num=_magnetic_kernel_i(row, "p_num"),
        pz_max=_magnetic_kernel_f(row, "pz_max_fm_inv"),
        imc=imc,
        route=:mfir,
        zeta_num=_magnetic_kernel_i(row, "zeta_num"),
        params=params,
    )
    x = SVector{5, Float64}(
        _magnetic_kernel_f(row, "phi_u"),
        _magnetic_kernel_f(row, "phi_d"),
        _magnetic_kernel_f(row, "phi_s"),
        _magnetic_kernel_f(row, "Phi"),
        _magnetic_kernel_f(row, "PhiBar"),
    )
    mu = SVector{3, Float64}(fill(_magnetic_kernel_f(row, "muB_MeV") / ext_hc / 3.0, 3))
    components = Models.MagneticThermodynamics.calculate_magnetic_omega_components(
        x, mu, _magnetic_kernel_f(row, "T_MeV") / ext_hc, conf,
    )
    @test isfinite(components.omega)
    @test isapprox(
        components.omega,
        _magnetic_kernel_f(row, "expected_omega_fm4");
        rtol=_magnetic_kernel_f(row, "rtol"),
        atol=_magnetic_kernel_f(row, "atol"),
    )
end
