#!/usr/bin/env julia

if length(ARGS) < 2
    error("usage: julia --project=. scripts/perf/relaxtime/analyze_transport_point_trace_focus.jl <trace-a> <trace-b>")
end

const TRACE_A = normpath(ARGS[1])
const TRACE_B = normpath(ARGS[2])

const FOCUS_RE = r"TransportWorkflow|transport_coefficients|solve_transport_from_equilibrium|RelaxationTime|AFieldBuilder|ParameterTypes|ConfigLoader|GaussLegendre|FastGaussQuadrature|ForwardDiff|NLSolversBase|DifferentiationInterface|Models\.solve|Models\._solve_constraint_fixedmu|Models\.Conditions"

function load_focus(path::String)
    isfile(path) || error("missing trace file: $(path)")
    out = String[]
    for line in eachline(path)
        occursin(FOCUS_RE, line) || continue
        push!(out, strip(line))
    end
    return out
end

function bucket(line::String)
    if occursin(r"ConfigLoader|load_pnjl_config|Constants_PNJL", line)
        return :config_loader
    elseif occursin(r"GaussLegendre|FastGaussQuadrature", line)
        return :quadrature
    elseif occursin(r"TransportWorkflow|transport_coefficients|solve_transport_from_equilibrium|TransportIntegrationConfig", line)
        return :transport_api
    elseif occursin(r"RelaxTime|AFieldBuilder|ParameterTypes", line)
        return :relaxtime_runtime
    elseif occursin(r"ForwardDiff|DifferentiationInterface|NLSolversBase|Models\._solve_constraint_fixedmu|Models\.Conditions|Models\.solve", line)
        return :solver_ad
    end
    return :other
end

function summarize(lines::Vector{String})
    counts = Dict{Symbol, Int}()
    for line in lines
        key = bucket(line)
        counts[key] = get(counts, key, 0) + 1
    end
    return counts
end

focus_a = load_focus(TRACE_A)
focus_b = load_focus(TRACE_B)
set_a = Set(focus_a)
set_b = Set(focus_b)

common = sort!(collect(intersect(set_a, set_b)))
only_a = sort!(collect(setdiff(set_a, set_b)))
only_b = sort!(collect(setdiff(set_b, set_a)))

println("[transport-point-focus] A=$(TRACE_A)")
println("[transport-point-focus] B=$(TRACE_B)")
println("  focus_a=$(length(focus_a)) unique_a=$(length(set_a))")
println("  focus_b=$(length(focus_b)) unique_b=$(length(set_b))")
println("  common=$(length(common)) only_a=$(length(only_a)) only_b=$(length(only_b))")

println("\n[buckets.common]")
for key in sort!(collect(keys(summarize(common))); by=string)
    println("  $(key) = $(summarize(common)[key])")
end

println("\n[buckets.only_a]")
for key in sort!(collect(keys(summarize(only_a))); by=string)
    println("  $(key) = $(summarize(only_a)[key])")
end

println("\n[buckets.only_b]")
for key in sort!(collect(keys(summarize(only_b))); by=string)
    println("  $(key) = $(summarize(only_b)[key])")
end

println("\n[common.lines]")
for line in common
    println(line)
end

if !isempty(only_a)
    println("\n[only_a.lines]")
    for line in only_a
        println(line)
    end
end

if !isempty(only_b)
    println("\n[only_b.lines]")
    for line in only_b
        println(line)
    end
end
