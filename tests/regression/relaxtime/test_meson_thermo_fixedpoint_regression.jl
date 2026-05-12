using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const BASELINE_PATH = joinpath(PROJECT_ROOT, "tests", "baselines", "relaxtime", "baseline_meson_thermo_fixedpoints_v1.csv")

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

function _load_rows(path::String)
    isfile(path) || error("baseline CSV not found: $path")
    rows = NamedTuple[] 
    lines = readlines(path)
    isempty(lines) && error("baseline CSV is empty: $path")

    header = [strip(col) for col in split(first(lines), ',')]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue
        cols = split(s, ','; keepempty=true)
        length(cols) == length(header) || error("invalid baseline row: $line")

        raw = Dict{String,String}(header[i] => strip(cols[i]) for i in eachindex(header))
        push!(rows, (
            label=raw["label"],
            T_MeV=parse(Float64, raw["T_MeV"]),
            muB_MeV=parse(Float64, raw["muB_MeV"]),
            xi=parse(Float64, raw["xi"]),
            workflow=raw["workflow"],
            channel_set=raw["channel_set"],
            primary_channel=raw["primary_channel"],
            secondary_channel=raw["secondary_channel"],
            phase_shift_variant=raw["phase_shift_variant"],
            thermo_derivation_mode=raw["thermo_derivation_mode"],
            qmax=parse(Float64, raw["qmax"]),
            q_nodes=parse(Int, raw["q_nodes"]),
            omega_min=parse(Float64, raw["omega_min"]),
            omega_max=parse(Float64, raw["omega_max"]),
            omega_nodes=parse(Int, raw["omega_nodes"]),
            eta=parse(Float64, raw["eta"]),
            ld_cutoff=parse(Float64, raw["ld_cutoff"]),
            ld_cutoff_mode=raw["ld_cutoff_mode"],
            ld_threshold_mode=raw["ld_threshold_mode"],
            equilibrium_converged=lowercase(raw["equilibrium_converged"]) == "true",
            P_meson=parse(Float64, raw["P_meson"]),
            P_meson_qp=parse(Float64, raw["P_meson_qp"]),
            P_meson_ld=parse(Float64, raw["P_meson_ld"]),
            P_quark_meanfield=parse(Float64, raw["P_quark_meanfield"]),
            P_total=parse(Float64, raw["P_total"]),
            entropy=parse(Float64, raw["entropy"]),
            epsilon=parse(Float64, raw["epsilon"]),
            trace_anomaly=parse(Float64, raw["trace_anomaly"]),
            P_primary=parse(Float64, raw["P_primary"]),
            P_secondary=parse(Float64, raw["P_secondary"]),
            P_primary_qp=parse(Float64, raw["P_primary_qp"]),
            P_primary_ld=parse(Float64, raw["P_primary_ld"]),
            P_secondary_qp=parse(Float64, raw["P_secondary_qp"]),
            P_secondary_ld=parse(Float64, raw["P_secondary_ld"]),
        ))
    end
    return rows
end

function _run_row(row)
    T_fm = row.T_MeV / Main.Constants_PNJL.ħc_MeV_fm
    mu_fm = row.muB_MeV / (3.0 * Main.Constants_PNJL.ħc_MeV_fm)
    point = Models.solve_gap_and_phase_shift_meson_thermo_point(
        T_fm,
        mu_fm;
        xi=row.xi,
        mesons=(Symbol(row.primary_channel), Symbol(row.secondary_channel)),
        mixed_branch_align=:strict_sign_binding,
        p_num=8,
        t_num=4,
        solver_kwargs=(iterations=20,),
        mass_kwargs=(iterations=20,),
        thermo_kwargs=(;
            pi_channel=Symbol(row.primary_channel),
            k_channel=Symbol(row.secondary_channel),
            scheme=:current,
            qmax=row.qmax,
            q_nodes=row.q_nodes,
            omega_min=row.omega_min,
            omega_max=row.omega_max,
            omega_nodes=row.omega_nodes,
            eta=row.eta,
            ld_cutoff=row.ld_cutoff,
            ld_cutoff_mode=Symbol(row.ld_cutoff_mode),
            ld_threshold_mode=Symbol(row.ld_threshold_mode),
            p_num=8,
            t_num=4,
        ),
        allow_legacy_fd_fallback=false,
    )
    return Models.build_meson_thermo_contract_row(point)
end

function _assert_close(actual::Float64, expected::Float64; rtol::Float64=5e-3, atol::Float64=1e-8)
    @test isfinite(actual)
    @test isapprox(actual, expected; rtol=rtol, atol=atol)
end

@testset "Meson thermo fixedpoint regression" begin
    rows = _load_rows(BASELINE_PATH)

    for row in rows
        @testset "$(row.label)" begin
            res = _run_row(row)

            @test res.workflow == row.workflow
            @test replace(res.channel_set, "," => "|") == row.channel_set
            @test res.primary_channel == row.primary_channel
            @test res.secondary_channel == row.secondary_channel
            @test res.phase_shift_variant == row.phase_shift_variant
            @test res.thermo_derivation_mode == row.thermo_derivation_mode
            @test res.ld_cutoff_mode == row.ld_cutoff_mode
            @test res.ld_threshold_mode == row.ld_threshold_mode
            @test res.equilibrium_converged == row.equilibrium_converged

            _assert_close(res.P_meson, row.P_meson)
            _assert_close(res.P_meson_qp, row.P_meson_qp)
            _assert_close(res.P_meson_ld, row.P_meson_ld)
            _assert_close(res.P_quark_meanfield, row.P_quark_meanfield)
            _assert_close(res.P_total, row.P_total)
            _assert_close(res.P_primary, row.P_primary)
            _assert_close(res.P_secondary, row.P_secondary)
            _assert_close(res.P_primary_qp, row.P_primary_qp)
            _assert_close(res.P_primary_ld, row.P_primary_ld)
            _assert_close(res.P_secondary_qp, row.P_secondary_qp; atol=1e-7)
            _assert_close(res.P_secondary_ld, row.P_secondary_ld)

            _assert_close(res.entropy, row.entropy; rtol=1e-2)
            _assert_close(res.epsilon, row.epsilon; rtol=1e-2)
            _assert_close(res.trace_anomaly, row.trace_anomaly; rtol=1e-2)

            @test res.P_meson ≈ res.P_meson_qp + res.P_meson_ld rtol=1e-10 atol=1e-10
            @test res.P_primary ≈ res.P_primary_qp + res.P_primary_ld rtol=1e-10 atol=1e-10
            @test res.P_secondary ≈ res.P_secondary_qp + res.P_secondary_ld rtol=1e-10 atol=1e-10
            @test res.P_total ≈ res.P_quark_meanfield + res.P_meson rtol=1e-10 atol=1e-10
            @test res.ld_cutoff ≈ row.ld_cutoff rtol=0 atol=1e-12
        end
    end
end
