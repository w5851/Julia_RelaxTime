#!/usr/bin/env julia
# Audit current source includes; generated Mermaid files are presentation only.
module DependencyAudit

using Dates

module StaticGraph
include(joinpath(@__DIR__, "gen_deps.jl"))
end

# Root compatibility shims, constants, types and config form the foundation.
const ALLOWED = Dict(
    "base" => Set(["base"]),
    "utils" => Set(["base", "utils"]),
    "integration" => Set(["base", "utils", "integration"]),
    "models" => Set(["base", "utils", "integration", "models"]),
    "relaxtime" => Set(["base", "utils", "integration", "relaxtime"]),
    "simulation" => Set(["base", "utils", "integration", "simulation"]),
)

# Existing bridges are visible debt, not permission for new src -> scripts edges.
const SCRIPT_BRIDGES = Dict(
    ("src/models/workflow_apps/ChargedGBUResearchWorkflow.jl",
     "scripts/analysis/relaxtime/causal_gbu_infinite_qgate.jl") => "charged GBU 计算内核迁移待独立验证",
    ("src/models/workflow_apps/ChargedGBUResearchWorkflow.jl",
     "scripts/relaxtime/workflow/charged_gbu_plot.jl") => "现有 charged GBU 绘图桥接",
    ("src/models/workflow_engine/adapters/RelaxtimeOrchestratorAdapter.jl",
     "scripts/relaxtime/config/WorkflowConfig.jl") => "现有 orchestrator 配置适配",
    ("src/models/workflow_engine/adapters/RelaxtimeOrchestratorAdapter.jl",
     "scripts/relaxtime/config/WorkflowConfigAudit.jl") => "现有 orchestrator 配置审计适配",
    ("src/models/workflow_engine/adapters/RelaxtimeOrchestratorAdapter.jl",
     "scripts/relaxtime/workflow/cross_section_orchestrated.jl") => "现有 orchestrator 工作流适配",
)

function dependency_group(path::String)
    group = StaticGraph.top_level_group(path)
    return group in ("root", "config", "constants", "types") ? "base" : group
end

function allowed_edge(from::String, to::String)
    ga, gb = dependency_group(from), dependency_group(to)
    gb in get(ALLOWED, ga, Set{String}()) && return true
    startswith(from, "src/models/workflow_apps/") && gb == "relaxtime" && return true
    # HTTP startup enters Models through its facade only.
    startswith(from, "src/simulation/") && to == "src/models/Models.jl" && return true
    return false
end

function audit_dependencies(root::String)
    unresolved = String[]
    adj = StaticGraph.collect_dependency_graph(root; unresolved)
    edges = sort([(a, b) for (a, bs) in adj for b in bs
                  if endswith(a, ".jl") && endswith(b, ".jl")])
    cross = [(a, b) for (a, b) in edges if dependency_group(a) != dependency_group(b)]
    bridges = [(a, b) for (a, b) in edges if haskey(SCRIPT_BRIDGES, (a, b))]
    violations = [(a, b) for (a, b) in edges
                  if !allowed_edge(a, b) && !haskey(SCRIPT_BRIDGES, (a, b))]
    return (; edges, cross, bridges, violations, unresolved=sort(unique(unresolved)))
end

function write_review(io::IO, result)
    println(io, "# 依赖审计报告\n")
    println(io, "生成时间：", Dates.now(), "\n")
    println(io, "来源：当前 src/ 源码；直接调用 gen_deps 的静态解析器，不读取旧 dependencies.mmd。\n")
    println(io, "只审计可定位的 include 文件边；相对导入名、宏展开、运行时参数和被加载脚本内部的依赖不在完整覆盖范围。")
    println(io, "解析 include 边：", length(result.edges), "；跨组边：", length(result.cross),
            "；已知脚本桥接：", length(result.bridges), "；违规：", length(result.violations),
            "；未解析或缺失：", length(result.unresolved), "。\n")
    println(io, "## 跨组依赖\n")
    for (a, b) in result.cross
        println(io, "- ", a, " → ", b)
    end
    isempty(result.cross) && println(io, "- 无。")
    println(io, "\n## 已知脚本桥接（后续独立迁移）\n")
    for edge in result.bridges
        println(io, "- ", edge[1], " → ", edge[2], "：", SCRIPT_BRIDGES[edge])
    end
    isempty(result.bridges) && println(io, "- 无。")
    println(io, "\n## 未解析或缺失的 include\n")
    for item in result.unresolved
        println(io, "- ", item)
    end
    isempty(result.unresolved) && println(io, "- 无。")
    println(io, "\n## 已解析 include 的矩阵违规\n")
    for (a, b) in result.violations
        println(io, "- ", a, " → ", b)
    end
    isempty(result.violations) && println(io, "- 未发现新增违规；这不是完整运行时依赖无环的证明。")
    println(io, "\n规则与边界见 [依赖规则](dependency_rules.md)。")
end

function main(; root::String=normpath(joinpath(@__DIR__, "..", "..")),
              strict::Bool=get(ENV, "DEPS_STRICT", "0") in ("1", "true", "TRUE", "yes", "YES"),
              io::IO=stdout)
    result = audit_dependencies(root)
    outdoc = joinpath(root, "docs", "architecture", "dependency_review.md")
    mkpath(dirname(outdoc))
    open(outdoc, "w") do output
        write_review(output, result)
    end
    println(io, "Wrote dependency review to: ", outdoc)
    println(io, "Resolved edges=$(length(result.edges)) violations=$(length(result.violations)) ",
            "known_bridges=$(length(result.bridges)) unresolved=$(length(result.unresolved))")
    !isempty(result.unresolved) && println(io, "Unresolved dynamic includes require manual review; coverage is partial.")
    source_errors = any(x -> occursin(": missing include target ", x) || occursin(": parse error:", x),
                        result.unresolved)
    return strict && (!isempty(result.violations) || source_errors) ? 1 : 0
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    exit(DependencyAudit.main())
end
