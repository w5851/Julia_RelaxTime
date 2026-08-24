#!/usr/bin/env julia

"""Solver-free RS/phase-reference consumer smoke for Issue #130.

This runner exercises the real gap-transport and phase-guided reference
resolution entrypoints for runtime, diagnostic, and explicit legacy modes. It
does not call an equilibrium solver or write any result/reference table.
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
using JSON3

const OUTPUT_PATH = isempty(ARGS) ?
    joinpath(PROJECT_ROOT, "docs", "analysis", "pnjl", "phase_reference",
        "issue130_rs_transport_runtime_parity_v1", "consumer_source_smoke.json") :
    normpath(abspath(ARGS[1]))

include(joinpath(PROJECT_ROOT, "scripts", "relaxtime", "run_phase_guided_transport_scan.jl"))

const CandidateCLI = Main.PhaseGuidedTransportScanCLI
const CandidatePlan = Main.PhaseGuidedTransportScanPlan
const PhaseEquilibrium = Main.GapTransportScanPhaseEquilibrium
const Adapter = Main.PhaseReferenceAdapter

function _dict_value(value)
    value isa Dict && return Dict(string(key) => _dict_value(item) for (key, item) in value)
    value isa NamedTuple && return Dict(string(key) => _dict_value(getproperty(value, key)) for key in keys(value))
    value isa AbstractVector && return [_dict_value(item) for item in value]
    value isa Symbol && return String(value)
    if value isa AbstractString || value isa Number || value isa Bool || value === nothing || value === missing
        return value
    end
    return string(value)
end

function _summary(source)
    summary = Adapter.source_summary(source)
    root_display = isempty(source.root) ? "" : replace(relpath(source.root, PROJECT_ROOT), '\\' => '/')
    return Dict(
        "source_kind" => String(Adapter.source_kind(source)),
        "source_layer" => String(Adapter.source_layer(source)),
        "runtime_enabled" => Adapter.source_runtime_enabled(source),
        "root" => root_display,
        "diagnostics" => _dict_value(summary),
    )
end

function _point_rows(source, mode::Symbol)
    mode_name = String(mode)
    phase_mode = mode
    boundary = PhaseEquilibrium.load_phase_boundary_data(
        0.0;
        phase_reference=source,
        phase_reference_mode=phase_mode,
    )
    crossover_T, crossover_xi = PhaseEquilibrium.interpolate_crossover_temperature(
        0.0,
        50.0;
        phase_reference=source,
        phase_reference_mode=phase_mode,
    )
    phase_kind, phase_T = CandidatePlan._phase_reference_for_mode_a(
        150.0,
        0.0;
        phase_reference=source,
        phase_reference_mode=phase_mode,
    )
    first_order_T, first_order_xi = CandidatePlan._interpolate_first_order_temperature(
        150.0,
        0.0;
        phase_reference=source,
        phase_reference_mode=phase_mode,
    )
    return Dict(
        "mode" => mode_name,
        "source_kind" => String(Adapter.source_kind(source)),
        "boundary_xi_used" => boundary.xi,
        "boundary_rows" => length(boundary.T_values),
        "boundary_first_T_MeV" => isempty(boundary.T_values) ? nothing : first(boundary.T_values),
        "boundary_first_muq_MeV" => isempty(boundary.mu_values) ? nothing : first(boundary.mu_values),
        "boundary_muB_CEP_MeV" => boundary.muB_CEP,
        "crossover_query_muq_MeV" => 50.0,
        "crossover_T_MeV" => crossover_T,
        "crossover_xi_used" => crossover_xi,
        "mode_a_query_muB_MeV" => 150.0,
        "mode_a_phase_kind" => String(phase_kind),
        "mode_a_phase_T_MeV" => phase_T,
        "first_order_T_MeV" => first_order_T,
        "first_order_xi_used" => first_order_xi,
        "solver_called" => false,
    )
end

function _source_for(mode::Symbol)
    args = String[
        "--mode", "fixed-muB-phase-scaled",
        "--case-name", "issue130_rs_parity_v1",
        "--xi-list", "0.0",
        "--muB-list", "150.0",
        "--alphaT-list", "1.0",
        "--dry-run",
        "--phase-reference-mode", String(mode),
    ]
    if mode === :diagnostic
        push!(args, "--phase-reference-root")
        push!(args, joinpath(PROJECT_ROOT, "data", "reference", "pnjl", "issue130_phase_reference_v1"))
    end
    opts = CandidateCLI.parse_args(args)
    return Main._load_runtime_phase_reference(opts)
end

function main()
    modes = (:runtime, :diagnostic, :legacy)
    sources = Dict{String,Any}()
    points = Any[]
    for mode in modes
        source = _source_for(mode)
        sources[String(mode)] = _summary(source)
        push!(points, _point_rows(source, mode))
    end

    payload = Dict(
        "schema_version" => "pnjl_issue130_rs_runtime_consumer_smoke_v1",
        "scope" => "solver_free_reference_resolution",
        "repo_head" => strip(read(`git -C $PROJECT_ROOT rev-parse HEAD`, String)),
        "candidate_root" => "data/reference/pnjl/issue130_phase_reference_v1",
        "legacy_root" => "data/reference/pnjl/{boundary.csv,cep.csv,crossover_dense.csv,spinodals.csv}",
        "query" => Dict("xi" => 0.0, "muB_MeV" => 150.0, "muq_MeV" => 50.0, "alpha_T" => 1.0),
        "solver_called" => false,
        "sources" => sources,
        "consumer_points" => points,
        "non_goals" => [
            "no equilibrium solver invocation",
            "no transport coefficient evaluation",
            "no data/reference write",
            "no legacy deletion or runtime fallback removal",
        ],
    )
    mkpath(dirname(OUTPUT_PATH))
    open(OUTPUT_PATH, "w") do io
        JSON3.pretty(io, payload)
        println(io)
    end
    println("Wrote solver-free consumer smoke: $(OUTPUT_PATH)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
