using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
const SOLVER_ROOT = joinpath(PROJECT_ROOT, "src", "models", "solver")

function _domain_for_path(path::String)
    rel = replace(relpath(path, SOLVER_ROOT), '\\' => '/')
    parts = split(rel, '/')
    if length(parts) >= 2 && parts[1] in ("api", "orchestrator", "spec", "runtime", "governance", "compat", "diagnostics", "config")
        return Symbol(parts[1])
    end
    return :root
end

function _solver_files()
    files = String[]
    for (root, _, names) in walkdir(SOLVER_ROOT)
        for name in names
            endswith(name, ".jl") || continue
            push!(files, joinpath(root, name))
        end
    end
    sort!(files)
    return files
end

function _occurs_identifier(source::String, ident::String)
    return occursin(Regex("\\b" * ident * "\\b"), source)
end

function _find_references(files::Vector{String}, symbols::Vector{String})
    refs = Dict{String, Vector{String}}()
    for file in files
        source = read(file, String)
        hits = String[]
        for sym in symbols
            _occurs_identifier(source, sym) || continue
            push!(hits, sym)
        end
        isempty(hits) || (refs[file] = hits)
    end
    return refs
end

const API_SYMBOLS = [
    "SolverResult",
    "SOLVER_CONTRACT_VERSION_V1",
    "coerce_solver_result",
    "solver_result_view",
    "solver_result_is_success",
]

const SPEC_SYMBOLS = [
    "ProblemSpec",
    "build_problem_spec",
    "ExtraConstraints",
]

const GOVERNANCE_SYMBOLS = [
    "build_governance_candidate",
    "execute_governance_selector",
    "normalize_selector_candidates",
    "evaluate_candidate_success",
    "execute_attempt_pool",
    "build_seed_pool",
    "classify_attempt_error",
    "normalize_error_message",
    "governance_quality_tag",
]

const RUNTIME_SYMBOLS = [
    "select_pressure_max_candidate",
    "select_residual_min_candidate",
    "default_hard_constraint_rules",
    "evaluate_hard_constraints",
]

const COMPAT_SYMBOLS = [
    "VarSchema",
    "SchemaRegistry",
    "register_schema!",
    "schema_for",
    "validate_schema",
    "named_to_vec",
    "vec_to_named",
    "build_pnjl_fixedmu_adapters",
    "build_pnjl_flavor_mu_adapters",
    "create_implicit_gap_solver",
    "create_flavor_mu_implicit_gap_solver",
    "create_pnjl_implicit_solver",
]

function _assert_no_refs(refs::Dict{String, Vector{String}}, predicate::Function; label::String)
    violations = String[]
    for (file, hits) in sort!(collect(refs); by=first)
        predicate(file) || continue
        push!(violations, "$(relpath(file, PROJECT_ROOT)): $(join(hits, ", "))")
    end
    @test isempty(violations)
    if !isempty(violations)
        @info label violations
    end
end

@testset "solver boundary rules R1-R5" begin
    files = _solver_files()
    api_refs = _find_references(files, API_SYMBOLS)
    governance_refs = _find_references(files, GOVERNANCE_SYMBOLS)
    runtime_refs = _find_references(files, RUNTIME_SYMBOLS)
    compat_refs = _find_references(files, COMPAT_SYMBOLS)

    @testset "R1 api -> orchestrator/spec only one-way" begin
        _assert_no_refs(api_refs, file -> begin
            d = _domain_for_path(file)
            d == :orchestrator || d == :spec
        end; label="R1 forbids orchestrator/spec -> api")
    end

    @testset "R2 runtime cannot depend on api" begin
        _assert_no_refs(api_refs, file -> _domain_for_path(file) == :runtime; label="R2 forbids runtime -> api")
    end

    @testset "R3 governance cannot depend on ConstraintSolver*" begin
        violations = String[]
        for file in files
            _domain_for_path(file) == :governance || continue
            src = read(file, String)
            bad = String[]
            occursin(r"\bConstraintSolver\w*\b", src) && push!(bad, "ConstraintSolver*")
            occursin(r"\b_solve_constraint_", src) && push!(bad, "_solve_constraint_*")
            isempty(bad) || push!(violations, "$(relpath(file, PROJECT_ROOT)): $(join(bad, ", "))")
        end
        @test isempty(violations)
        if !isempty(violations)
            @info "R3 violations" violations
        end
    end

    @testset "R4 compat only api/orchestrator may reference" begin
        _assert_no_refs(compat_refs, file -> begin
            d = _domain_for_path(file)
            d in (:spec, :governance, :runtime, :diagnostics, :config)
        end; label="R4 forbids compat refs from non-api/orchestrator domains")
    end

    @testset "R5 only orchestrator may bridge runtime + governance" begin
        violations = String[]
        for file in files
            d = _domain_for_path(file)
            d == :orchestrator && continue
            has_runtime = haskey(runtime_refs, file)
            has_governance = haskey(governance_refs, file)
            (has_runtime && has_governance) || continue
            push!(violations, "$(relpath(file, PROJECT_ROOT)): runtime=$(join(runtime_refs[file], ", ")) governance=$(join(governance_refs[file], ", "))")
        end
        @test isempty(violations)
        if !isempty(violations)
            @info "R5 violations" violations
        end
    end
end
