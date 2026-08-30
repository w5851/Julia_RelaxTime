using Test

const FORMULA_ROUTE_PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const FORMULA_ROUTE_CHECKER = joinpath(
    FORMULA_ROUTE_PROJECT_ROOT,
    "scripts",
    "dev",
    "check_formula_route_closure.jl",
)

if !isdefined(Main, :FormulaRouteClosure)
    include(FORMULA_ROUTE_CHECKER)
end

function _write_formula_route_fixture(root::AbstractString; production_authorized::Bool=false, marker::Bool=true)
    document_rel = "docs/reference/formula/fixture_route.md"
    test_rel = "tests/unit/config/fixture_formula_route_test.jl"
    document = joinpath(root, split(document_rel, '/')...)
    test_file = joinpath(root, split(test_rel, '/')...)
    mkpath(dirname(document))
    mkpath(dirname(test_file))
    markers = marker ? "route_id: fixture_route\n公式闭合链\n外部来源与项目转换\n公式 → 代码 → 测试映射\n生产边界与升格条件\n未决项与审查问题\n10.1234/example\n" : "route_id: fixture_route\n"
    write(document, markers)
    write(test_file, "# fixture\n")

    registry = """
schema_version = "v1"
allowed_statuses = ["draft", "candidate", "production_authorized", "deprecated"]
required_fields = ["id", "document", "status", "scope", "model_start", "approximations", "external_sources", "formula_tests", "production_authorized", "unresolved_items"]

[[route]]
id = "fixture_route"
document = "$(document_rel)"
status = "candidate"
scope = "fixture"
model_start = "fixture_model"
approximations = ["fixture approximation"]
external_sources = ["10.1234/example"]
formula_tests = ["$(test_rel)"]
production_authorized = $(production_authorized)
unresolved_items = ["fixture unresolved item"]
required_document_markers = ["route_id: fixture_route", "公式闭合链", "外部来源与项目转换", "公式 → 代码 → 测试映射", "生产边界与升格条件", "未决项与审查问题"]
"""
    registry_path = joinpath(root, "config", "governance", "formula_route_closure.toml")
    mkpath(dirname(registry_path))
    write(registry_path, registry)
    return registry_path
end

@testset "formula route closure registry" begin
    @test isfile(FORMULA_ROUTE_CHECKER)
    @test isfile(joinpath(FORMULA_ROUTE_PROJECT_ROOT, "config", "governance", "formula_route_closure.toml"))
    @test isempty(Main.FormulaRouteClosure.validate_registry(FORMULA_ROUTE_PROJECT_ROOT))

    mktempdir() do temp_root
        valid_registry = _write_formula_route_fixture(temp_root)
        @test isempty(Main.FormulaRouteClosure.validate_registry(
            temp_root;
            registry_rel=relpath(valid_registry, temp_root),
        ))

        invalid_registry = _write_formula_route_fixture(temp_root; production_authorized=true)
        invalid = Main.FormulaRouteClosure.validate_registry(
            temp_root;
            registry_rel=relpath(invalid_registry, temp_root),
        )
        @test any(item -> occursin("cannot authorize production", item), invalid)

        invalid_status_registry = _write_formula_route_fixture(temp_root)
        invalid_status = replace(
            read(invalid_status_registry, String),
            "status = \"candidate\"" => "status = \"not_a_route_state\"",
        )
        write(invalid_status_registry, invalid_status)
        status_violations = Main.FormulaRouteClosure.validate_registry(
            temp_root;
            registry_rel=relpath(invalid_status_registry, temp_root),
        )
        @test any(item -> occursin("status is not allowed", item), status_violations)

        escaped_document_registry = _write_formula_route_fixture(temp_root)
        escaped_document = replace(
            read(escaped_document_registry, String),
            "docs/reference/formula/fixture_route.md" => "docs/reference/outside_route.md",
        )
        write(escaped_document_registry, escaped_document)
        path_violations = Main.FormulaRouteClosure.validate_registry(
            temp_root;
            registry_rel=relpath(escaped_document_registry, temp_root),
        )
        @test any(item -> occursin("must stay under docs/reference/formula/", item), path_violations)

        missing_test_registry = _write_formula_route_fixture(temp_root)
        missing_test = replace(
            read(missing_test_registry, String),
            "tests/unit/config/fixture_formula_route_test.jl" =>
                "tests/unit/config/missing_formula_route_test.jl",
        )
        write(missing_test_registry, missing_test)
        missing_test_violations = Main.FormulaRouteClosure.validate_registry(
            temp_root;
            registry_rel=relpath(missing_test_registry, temp_root),
        )
        @test any(item -> occursin("formula_tests entry does not exist", item), missing_test_violations)

        missing_source_registry = _write_formula_route_fixture(temp_root)
        missing_source_document = replace(read(
            joinpath(temp_root, "docs", "reference", "formula", "fixture_route.md"),
            String,
        ), "10.1234/example" => "10.1234/omitted")
        write(
            joinpath(temp_root, "docs", "reference", "formula", "fixture_route.md"),
            missing_source_document,
        )
        source_violations = Main.FormulaRouteClosure.validate_registry(
            temp_root;
            registry_rel=relpath(missing_source_registry, temp_root),
        )
        @test any(item -> occursin("does not cite registered source", item), source_violations)

        missing_marker_registry = _write_formula_route_fixture(temp_root; marker=false)
        missing_marker = Main.FormulaRouteClosure.validate_registry(
            temp_root;
            registry_rel=relpath(missing_marker_registry, temp_root),
        )
        @test any(item -> occursin("missing required marker", item), missing_marker)
    end
end
