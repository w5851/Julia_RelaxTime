#!/usr/bin/env julia
# Analyze cross-directory dependencies from Mermaid .mmd

using Dates

root = pwd()
indoc = joinpath(root, "docs", "architecture", "dependencies.mmd")
outdoc = joinpath(root, "docs", "architecture", "dependency_review.md")
strict = get(ENV, "DEPS_STRICT", "0") in ("1", "true", "TRUE", "yes", "YES")

if !isfile(indoc)
    error("Missing dependencies.mmd at: $indoc")
end

lines = readlines(indoc)

# parse nodes: id[label]
node_label = Dict{String,String}()
node_re = r"^\s*([A-Za-z0-9_]+)\[(.+)\]\s*$"
edge_re = r"^\s*([A-Za-z0-9_]+)\s*-->\s*([A-Za-z0-9_]+)\s*$"

for line in lines
    m = match(node_re, line)
    if m !== nothing
        node_label[m.captures[1]] = m.captures[2]
    end
end

function group_of(label::String)
    if occursin("/", label)
        return split(label, "/")[1]
    end
    return "root"
end

function is_file_label(label::String)
    return occursin("/", label)
end

# build edges
edges = Vector{Tuple{String,String,String,String}}()
for line in lines
    m = match(edge_re, line)
    if m === nothing
        continue
    end
    a = m.captures[1]
    b = m.captures[2]
    if !haskey(node_label, a) || !haskey(node_label, b)
        continue
    end
    la = node_label[a]
    lb = node_label[b]
    if !(is_file_label(la) && is_file_label(lb))
        continue
    end
    ga = group_of(la)
    gb = group_of(lb)
    push!(edges, (ga, gb, la, lb))
end

# allowed matrix
allowed = Dict(
    "root" => Set(["root"]),
    "utils" => Set(["root", "utils"]),
    "integration" => Set(["root", "utils", "integration"]),
    "simulation" => Set(["root", "utils", "integration", "simulation"]),
    "pnjl" => Set(["root", "utils", "integration", "pnjl", "relaxtime"]),
    "relaxtime" => Set(["root", "utils", "integration", "relaxtime"])
)

function is_exception(from_label::String, to_group::String)
    # allow pnjl/workflows -> relaxtime
    return occursin("pnjl/workflows/", from_label) && to_group == "relaxtime"
end

# summarize
cross = [(ga, gb, la, lb) for (ga, gb, la, lb) in edges if ga != gb]

# count edges by group pair
pair_counts = Dict{Tuple{String,String}, Int}()
for (ga, gb, _la, _lb) in cross
    pair_counts[(ga, gb)] = get(pair_counts, (ga, gb), 0) + 1
end

# violations
violations = Vector{Tuple{String,String,String,String}}()
for (ga, gb, la, lb) in cross
    if !(gb in get(allowed, ga, Set{String}()))
        if !is_exception(la, gb)
            push!(violations, (ga, gb, la, lb))
        end
    end
end

# render report
open(outdoc, "w") do io
    println(io, "# 依赖审计报告")
    println(io, "生成时间：", Dates.now())
    println(io, "\n来源：docs/architecture/dependencies.mmd\n")

    println(io, "## 跨目录依赖清单（汇总）\n")
    if isempty(pair_counts)
        println(io, "- 无跨目录依赖\n")
    else
        for (pair, cnt) in sort(collect(pair_counts); by=x->x[1])
            println(io, "- ", pair[1], " -> ", pair[2], ": ", cnt)
        end
        println(io)
    end

    println(io, "## 跨目录依赖明细\n")
    if isempty(cross)
        println(io, "- 无\n")
    else
        for (ga, gb, la, lb) in sort(cross)
            println(io, "- ", la, " -> ", lb, " (", ga, " -> ", gb, ")")
        end
        println(io)
    end

    println(io, "## 违规点（基于依赖矩阵）\n")
    if isempty(violations)
        println(io, "- 未发现违规\n")
    else
        for (ga, gb, la, lb) in sort(violations)
            println(io, "- ", la, " -> ", lb, " (", ga, " -> ", gb, ")")
        end
        println(io)
    end

    println(io, "## 调整建议\n")
    if isempty(violations)
        println(io, "- 当前依赖符合矩阵，建议继续保持。\n")
    else
        println(io, "- 若出现“底层依赖上层”，优先考虑下沉公共逻辑到 `utils/` 或 `integration/`。")
        println(io, "- 若属于流程编排层，可集中在 `workflows/` 并作为例外记录。")
        println(io, "- 如果仅为常量/类型共享，考虑抽离到 `src/Constants_PNJL.jl` 或新增 `src/core/`。")
        println(io)
    end
end

println("Wrote dependency review to: ", outdoc)

if strict && !isempty(violations)
    println("Dependency violations found in strict mode. Failing.")
    exit(1)
end
