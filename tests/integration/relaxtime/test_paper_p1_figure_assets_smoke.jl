using Test
using JSON3

const REPO_ROOT_P1 = normpath(joinpath(@__DIR__, "..", "..", ".."))
const P1_SCRIPT_PATH = joinpath(REPO_ROOT_P1, "scripts", "relaxtime", "build_paper_p1_figure_assets.py")

function _python_cmd()
    py = Sys.which("python")
    py !== nothing && return py
    py3 = Sys.which("python3")
    py3 !== nothing && return py3
    return nothing
end

function _read_plain_csv(path::String)
    rows = Vector{Dict{String,String}}()
    header = String[]
    open(path, "r") do io
        for line in eachline(io)
            s = strip(line)
            isempty(s) && continue
            startswith(s, "#") && continue
            if isempty(header)
                header = [strip(x) for x in split(s, ',')]
                continue
            end
            vals = split(s, ',')
            row = Dict{String,String}()
            for (idx, key) in enumerate(header)
                row[key] = idx <= length(vals) ? strip(vals[idx]) : ""
            end
            push!(rows, row)
        end
    end
    return header, rows
end

@testset "P1 paper figure asset builder smoke" begin
    @test isfile(P1_SCRIPT_PATH)

    python = _python_cmd()
    if python === nothing
        @test true
        return
    end

    tmp = mktempdir()
    mott_csv = joinpath(tmp, "mott_grid.csv")
    isentropic_csv = joinpath(tmp, "isentropic.csv")
    phase_dir = joinpath(tmp, "phase")
    phase_ref_root = joinpath(tmp, "phase_reference")
    out_dir = joinpath(tmp, "assets")
    mkpath(phase_dir)
    mkpath(phase_ref_root)

    open(mott_csv, "w") do io
        println(io, "T_MeV,muB_MeV,xi,M_pi,M_K,Gamma_pi,Gamma_K,threshold_pi,threshold_K,gap_pi,gap_K")
        for (muB, tpi, tk) in ((0.0, 150.0, 155.0), (120.0, 152.0, 158.0))
            for T in (140.0, 160.0)
                gap_pi = T - tpi
                gap_k = T - tk
                println(io, join([T, muB, 0.0, 0.70, 0.80, 0.01, 0.02, 0.70 - gap_pi, 0.80 - gap_k, gap_pi, gap_k], ','))
            end
        end
    end
    write(joinpath(tmp, "mott_grid_combined_manifest.json"), JSON3.write(Dict(
        "schema_version" => "paper_p1_mott_combined_v1",
        "scan" => Dict("mott_phase" => Dict(
            "equilibrium_branch_mode" => "stable",
            "equilibrium_selector_policy" => "pressure_max_under_constraints",
            "equilibrium_selector_tiebreak" => "residual_norm_then_seed_index",
        )),
    )))

    open(isentropic_csv, "w") do io
        println(io, "path_family,path_profile,path_segment,path_point_index,path_order_key,path_label,sigma_target,T_MeV,muB_MeV,xi,M_pi,M_K,Gamma_pi,Gamma_K,threshold_pi,threshold_K,gap_pi,gap_K")
        println(io, "isentropic,default,fixed_sigma,0,140.0,sigma=30.000000,30.0,140.0,10.0,0.0,0.70,0.80,0.01,0.02,1.70,2.80,-1.0,-2.0")
        println(io, "isentropic,default,fixed_sigma,1,160.0,sigma=30.000000,30.0,160.0,30.0,0.0,0.70,0.80,0.01,0.02,-0.30,-1.20,1.0,2.0")
    end

    write(joinpath(phase_dir, "phase_summary.json"), JSON3.write(Dict(
        "xi" => 0.0,
        "cep" => Dict(
            "found" => true,
            "T_cep_MeV" => 120.0,
            "muq_cep_MeV" => 80.0,
            "muB_cep_MeV" => 240.0,
        ),
    )))
    write(joinpath(phase_dir, "first_order_boundary.csv"), "T_MeV,mu_transition_MeV,rho_hadron,rho_quark,area_residual,converged\n100.0,50.0,0.1,0.3,0.0,true\n")
    write(joinpath(phase_dir, "spinodal.csv"), "T_MeV,mu_spinodal_hadron_MeV,mu_spinodal_quark_MeV,rho_spinodal_hadron,rho_spinodal_quark\n100.0,45.0,55.0,0.1,0.3\n")
    write(joinpath(phase_dir, "crossover_line.csv"), "mu_MeV,T_crossover_MeV,rho,method,converged,derivative,variable\n40.0,130.0,0.2,peak,true,1.0,phi_u\n")
    write(joinpath(phase_ref_root, "boundary_demo.csv"), "xi,T_MeV,mu_transition_MeV,rho_hadron,rho_quark,curve_parameter,plot_order_key\n0.3,110.0,52.0,0.1,0.3,110.0,110.0\n")
    write(joinpath(phase_ref_root, "spinodals_demo.csv"), "xi,T_MeV,mu_spinodal_hadron_MeV,mu_spinodal_quark_MeV,rho_spinodal_hadron,rho_spinodal_quark,curve_parameter,plot_order_key\n0.3,110.0,47.0,57.0,0.1,0.3,110.0,110.0\n")
    write(joinpath(phase_ref_root, "crossover_demo.csv"), "xi,mu_MeV,T_crossover_MeV,rho,method,converged,derivative,variable,curve_parameter,plot_order_key\n0.3,42.0,132.0,0.2,peak,true,1.0,phi_u,42.0,42.0\n")
    write(joinpath(phase_ref_root, "cep_demo.csv"), "xi,T_CEP_MeV,muq_CEP_MeV,muB_CEP_MeV\n0.3,121.0,81.0,243.0\n")

    run(`$python $P1_SCRIPT_PATH --mott-grid-csv $mott_csv --isentropic-csv $isentropic_csv --phase-dir $phase_dir --out-dir $out_dir --skip-plots`)

    @test isfile(joinpath(out_dir, "mott_lines.csv"))
    @test isfile(joinpath(out_dir, "isentropic_trajectories.csv"))
    @test isfile(joinpath(out_dir, "isentropic_mott_crossings.csv"))
    @test isfile(joinpath(out_dir, "phase_overlay.csv"))
    @test isfile(joinpath(out_dir, "figure_manifest.json"))

    _, mott_rows = _read_plain_csv(joinpath(out_dir, "mott_lines.csv"))
    @test haskey(first(mott_rows), "bracket_kind")
    @test all(haskey(row, "bracket_mass_jump_inv_fm") for row in mott_rows)
    pi_mu0 = only(filter(r -> r["channel"] == "pi" && parse(Float64, r["muB_MeV"]) == 0.0, mott_rows))
    @test parse(Float64, pi_mu0["T_Mott_MeV"]) ≈ 150.0
    k_mu120 = only(filter(r -> r["channel"] == "K" && parse(Float64, r["muB_MeV"]) == 120.0, mott_rows))
    @test parse(Float64, k_mu120["T_Mott_MeV"]) ≈ 158.0

    _, crossing_rows = _read_plain_csv(joinpath(out_dir, "isentropic_mott_crossings.csv"))
    path_pi = only(filter(r -> r["channel"] == "pi", crossing_rows))
    @test parse(Float64, path_pi["T_MeV"]) ≈ 150.0
    @test parse(Float64, path_pi["muB_MeV"]) ≈ 20.0

    _, phase_rows = _read_plain_csv(joinpath(out_dir, "phase_overlay.csv"))
    first_order = only(filter(r -> r["kind"] == "first_order", phase_rows))
    @test parse(Float64, first_order["muB_MeV"]) ≈ 150.0
    @test haskey(first_order, "curve_parameter")
    @test haskey(first_order, "plot_order_key")
    @test parse(Float64, first_order["plot_order_key"]) ≈ 100.0
    @test any(r -> r["kind"] == "cep" && parse(Float64, r["muB_MeV"]) ≈ 240.0, phase_rows)

    out_dir_ref = joinpath(tmp, "assets_reference")
    run(`$python $P1_SCRIPT_PATH --mott-grid-csv $mott_csv --isentropic-csv $isentropic_csv --phase-reference-root $phase_ref_root --phase-reference-tag demo --out-dir $out_dir_ref --skip-plots`)
    _, phase_ref_rows = _read_plain_csv(joinpath(out_dir_ref, "phase_overlay.csv"))
    @test any(r -> r["kind"] == "first_order" && parse(Float64, r["xi"]) ≈ 0.3 && parse(Float64, r["muB_MeV"]) ≈ 156.0, phase_ref_rows)
    @test all(haskey(row, "plot_order_key") for row in phase_ref_rows)
    manifest_ref = JSON3.read(read(joinpath(out_dir_ref, "figure_manifest.json"), String))
    @test manifest_ref["inputs"]["phase_reference_tag"] == "demo"

    manifest = JSON3.read(read(joinpath(out_dir, "figure_manifest.json"), String))
    @test manifest["schema_version"] == "paper_p1_assets_v1"
    @test manifest["inputs"]["mott_source"]["equilibrium_branch_mode"] == "stable"
    @test manifest["inputs"]["mott_source"]["equilibrium_selector_policy"] == "pressure_max_under_constraints"
    @test manifest["counts"]["mott_line_points"] == 4
    @test manifest["counts"]["isentropic_crossings"] == 2

    has_matplotlib = success(`$python -c "import matplotlib"`)
    if has_matplotlib
        out_dir_plot = joinpath(tmp, "assets_with_plots")
        run(`$python $P1_SCRIPT_PATH --mott-grid-csv $mott_csv --isentropic-csv $isentropic_csv --phase-dir $phase_dir --out-dir $out_dir_plot --formats png`)
        @test isfile(joinpath(out_dir_plot, "figures", "p1_mott_phase_diagram.png"))
        @test isfile(joinpath(out_dir_plot, "figures", "p1_mott_phase_diagram_xi_0.png"))
        @test isfile(joinpath(out_dir_plot, "figures", "p1_isentropic_mott_paths.png"))
    end
end
