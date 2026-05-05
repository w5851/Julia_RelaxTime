#!/usr/bin/env julia

if length(ARGS) < 2
    error("usage: julia --project=. scripts/perf/relaxtime/analyze_solver_ad_trace_focus.jl <trace-a> <trace-b>")
end

const TRACE_A = normpath(ARGS[1])
const TRACE_B = normpath(ARGS[2])

const FOCUS_RE = r"ForwardDiff|DifferentiationInterface|NLSolversBase|typeof\(Models\.solve\)|Models\._solve_constraint_fixedmu|Models\.Conditions"

function load_focus(path::String)
    isfile(path) || error("missing trace file: $(path)")
    raw_lines = String[]
    for line in eachline(path)
        occursin(FOCUS_RE, line) || continue
        push!(raw_lines, strip(line))
    end
    return raw_lines
end

function bucket(line::String)
    if occursin(r"DifferentiationInterface\.prepare_jacobian|ForwardDiffTwoArgJacobianPrep|JacobianConfig", line)
        return :prepare_jacobian
    elseif occursin(r"NLSolversBase|OnceDifferentiable", line)
        return :once_differentiable
    elseif occursin(r"ForwardDiff\.Tag", line)
        return :forwarddiff_tag
    elseif occursin(r"ForwardDiff\.Dual", line)
        return :forwarddiff_dual
    elseif occursin(r"ForwardDiff\.Partials", line)
        return :forwarddiff_partials
    elseif occursin(r"typeof\(Models\.solve\)|Models\._solve_constraint_fixedmu|Core\.kwcall", line)
        return :solve_kwcall
    elseif occursin(r"Models\.Conditions|_gap_conditions_with_model|residual!#140", line)
        return :conditions_residual
    elseif occursin(r"ForwardDiff", line)
        return :forwarddiff_misc
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

function tag_family(line::String)
    if occursin(r"Models\.Conditions\.var\"#_gap_conditions_with_model", line)
        return :conditions_gap
    elseif occursin(r"Models\.var\"#residual!#140", line)
        return :residual_140
    elseif occursin(r"Models\.var\"#129#130", line)
        return :solve_callback_129_130
    elseif occursin(r"Models\.var\"#256#257", line)
        return :thermo_probe_256_257
    elseif occursin(r"Models\.var\"#258#259", line)
        return :thermo_probe_258_259
    elseif occursin(r"ForwardDiff\.Tag", line)
        return :other_forwarddiff_tag
    end
    return :non_tag
end

function summarize_tag_families(lines::Vector{String})
    counts = Dict{Symbol, Int}()
    for line in lines
        family = tag_family(line)
        family == :non_tag && continue
        counts[family] = get(counts, family, 0) + 1
    end
    return counts
end

function print_bucket_summary(title::String, lines::Vector{String})
    counts = summarize(lines)
    println("\n[$(title)]")
    for key in sort!(collect(keys(counts)); by=string)
        println("  $(key) = $(counts[key])")
    end
end

function print_tag_family_summary(title::String, lines::Vector{String})
    counts = summarize_tag_families(lines)
    println("\n[$(title)]")
    for key in sort!(collect(keys(counts)); by=string)
        println("  $(key) = $(counts[key])")
    end
end

focus_a = load_focus(TRACE_A)
focus_b = load_focus(TRACE_B)
set_a = Set(focus_a)
set_b = Set(focus_b)

common = sort!(collect(intersect(set_a, set_b)))
only_a = sort!(collect(setdiff(set_a, set_b)))
only_b = sort!(collect(setdiff(set_b, set_a)))

println("[solver-ad-focus] A=$(TRACE_A)")
println("[solver-ad-focus] B=$(TRACE_B)")
println("  focus_a=$(length(focus_a)) unique_a=$(length(set_a))")
println("  focus_b=$(length(focus_b)) unique_b=$(length(set_b))")
println("  common=$(length(common)) only_a=$(length(only_a)) only_b=$(length(only_b))")

print_bucket_summary("buckets.common", common)
print_bucket_summary("buckets.only_a", only_a)
print_bucket_summary("buckets.only_b", only_b)
print_tag_family_summary("tag_families.common", common)
print_tag_family_summary("tag_families.only_a", only_a)
print_tag_family_summary("tag_families.only_b", only_b)

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
