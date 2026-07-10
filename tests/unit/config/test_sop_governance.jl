using Dates
using Test

const SOP_CHECKER_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "scripts", "dev", "check_sop_governance.jl"))

if !isdefined(Main, :SopGovernance)
    include(SOP_CHECKER_PATH)
end

function _write_text(path::AbstractString, content::AbstractString)
    mkpath(dirname(path))
    write(path, content)
    return path
end

function _sop_content(; omit_section::Union{Nothing, String}=nothing)
    sections = [
        "## 1. 目的与适用范围",
        "## 2. 非适用范围",
        "## 3. 权威入口",
        "## 4. 物理口径、单位与参数约束",
        "## 5. 输入配置及优先级",
        "## 6. 环境与版本冻结",
        "## 7. Smoke 预检",
        "## 8. 收敛性验证",
        "## 9. 正式计算命令",
        "## 10. 输出目录与产物合同",
        "## 11. Regression / Validation 验收",
        "## 12. 失败点、断点续算与重跑",
        "## 13. Diagnostic 与 Formal Production 的边界",
        "## 14. 后处理与作图",
        "## 15. 关联公式、API 和测试",
        "## 16. 最后验证记录",
    ]
    kept = omit_section === nothing ? sections : filter(!=(omit_section), sections)
    return "# Fixture SOP\n\n" * join((section * "\n\nfixture\n" for section in kept), "\n")
end

function _registry_text(; last_verified::String="2026-07-10", duplicate_authority::Bool=false)
    required = join(("  \"$(section)\"," for section in [
        "## 1. 目的与适用范围",
        "## 2. 非适用范围",
        "## 3. 权威入口",
        "## 4. 物理口径、单位与参数约束",
        "## 5. 输入配置及优先级",
        "## 6. 环境与版本冻结",
        "## 7. Smoke 预检",
        "## 8. 收敛性验证",
        "## 9. 正式计算命令",
        "## 10. 输出目录与产物合同",
        "## 11. Regression / Validation 验收",
        "## 12. 失败点、断点续算与重跑",
        "## 13. Diagnostic 与 Formal Production 的边界",
        "## 14. 后处理与作图",
        "## 15. 关联公式、API 和测试",
        "## 16. 最后验证记录",
    ]), "\n")

    second = duplicate_authority ? """

[[sop]]
id = "duplicate"
path = "docs/guides/sop/duplicate.md"
status = "active"
authoritative_for = ["fixture_workflow"]
entrypoints = ["scripts/run_fixture.jl"]
stable_entrypoints = ["scripts/run_fixture.jl"]
configs = []
review_cycle_days = 30
last_verified = "$(last_verified)"
verification_commands = ["julia --project=. scripts/run_fixture.jl --help"]
""" : ""

    return """
schema_version = "v1"
allowed_statuses = ["draft", "active", "deprecated"]
stable_entrypoint_index = "docs/guides/scripts/README.md"
required_sections = [
$(required)
]
forbidden_patterns = ["src/pnjl/"]

[[sop]]
id = "fixture"
path = "docs/guides/sop/fixture.md"
status = "active"
authoritative_for = ["fixture_workflow"]
entrypoints = ["scripts/run_fixture.jl"]
stable_entrypoints = ["scripts/run_fixture.jl"]
configs = []
review_cycle_days = 30
last_verified = "$(last_verified)"
verification_commands = ["julia --project=. scripts/run_fixture.jl --help"]
$(second)
"""
end

function _make_fixture(; omit_section=nothing, last_verified="2026-07-10", duplicate_authority=false)
    root = mktempdir()
    _write_text(joinpath(root, "scripts", "run_fixture.jl"), "println(\"fixture\")\n")
    _write_text(
        joinpath(root, "docs", "guides", "scripts", "README.md"),
        "scripts/run_fixture.jl\n",
    )
    _write_text(
        joinpath(root, "docs", "guides", "sop", "fixture.md"),
        _sop_content(; omit_section=omit_section),
    )
    if duplicate_authority
        _write_text(
            joinpath(root, "docs", "guides", "sop", "duplicate.md"),
            _sop_content(),
        )
    end
    _write_text(
        joinpath(root, "config", "governance", "docs_authority_map.toml"),
        _registry_text(; last_verified=last_verified, duplicate_authority=duplicate_authority),
    )
    return root
end

@testset "SOP governance registry" begin
    @test isfile(SOP_CHECKER_PATH)

    valid_root = _make_fixture()
    @test isempty(Main.SopGovernance.validate_registry(
        valid_root;
        current_date=Date(2026, 7, 10),
    ))

    missing_section = "## 8. 收敛性验证"
    incomplete_root = _make_fixture(; omit_section=missing_section)
    incomplete = Main.SopGovernance.validate_registry(
        incomplete_root;
        current_date=Date(2026, 7, 10),
    )
    @test any(item -> occursin("missing required section", item), incomplete)
    @test any(item -> occursin(missing_section, item), incomplete)

    overdue_root = _make_fixture(; last_verified="2026-01-01")
    overdue = Main.SopGovernance.validate_registry(
        overdue_root;
        current_date=Date(2026, 7, 10),
    )
    @test any(item -> occursin("review is overdue", item), overdue)

    duplicate_root = _make_fixture(; duplicate_authority=true)
    duplicate = Main.SopGovernance.validate_registry(
        duplicate_root;
        current_date=Date(2026, 7, 10),
    )
    @test any(item -> occursin("is claimed by both", item), duplicate)
end
