#!/usr/bin/env julia

"""
文档一致性检查（轻量门禁）

检查目标：
1) active/guides/api 主导航文档中，严格禁止误导性旧主线路径。
2) archived 历史文档允许保留旧路径，但禁止把旧路径表述成“当前/主线/统一入口”。
3) 当 README 标注“修复中”时，guides 不应出现绝对化状态词。
4) 当 README 标注“已验证可用”时，guides 不应残留“修复中”表述。

用法：
  julia --project=. scripts/dev/check_docs_consistency.jl
"""

function project_root()
    return normpath(joinpath(@__DIR__, "..", ".."))
end

function relpath_display(path::String, root::String)
    rp = relpath(path, root)
    return replace(rp, '\\' => '/')
end

function collect_markdown_files(dir::String)
    files = String[]
    if !isdir(dir)
        return files
    end
    for (root, _, names) in walkdir(dir)
        for name in names
            endswith(lowercase(name), ".md") || continue
            push!(files, joinpath(root, name))
        end
    end
    sort!(files)
    return files
end

function main()
    root = project_root()
    readme = joinpath(root, "README.md")
    docs_dir = joinpath(root, "docs")
    guides_dir = joinpath(docs_dir, "guides")
    api_dir = joinpath(docs_dir, "api")
    active_dir = joinpath(docs_dir, "dev", "active")
    archived_dir = joinpath(docs_dir, "dev", "archived")

    guide_files = collect_markdown_files(guides_dir)
    api_files = collect_markdown_files(api_dir)
    active_files = collect_markdown_files(active_dir)
    archived_files = collect_markdown_files(archived_dir)

    strict_targets = vcat(guide_files, api_files, active_files)

    violations = String[]

    legacy_checks_common = [
        ("test_unit/", "历史测试目录，应改为 tests/unit/"),
        (r"\bjulia\s+server\.jl\b", "历史服务启动命令，应改为 scripts/server/server*.jl"),
        (".\\start.bat", "历史根目录启动脚本，应改为 .\\scripts\\server\\start.bat"),
        ("doc/domain-knowledge/", "历史文档目录，应改为 docs/reference/domain-knowledge/"),
    ]

    strict_legacy_checks = [
        (r"src/pnjl/PNJL\.jl", "旧主线入口路径，应改为 Models + src/models/entrypoints.jl"),
        (r"src/pnjl/", "旧主线路径，应避免在现行主导航中作为入口引用"),
    ]

    misleading_entry_pattern = r"((当前|主线|统一|现行)入口(为|是).*(src/pnjl|PNJL\.jl))|((src/pnjl|PNJL\.jl).*(当前|主线|统一|现行)入口(为|是))"

    for file in strict_targets
        isfile(file) || continue
        content = read(file, String)
        display = relpath_display(file, root)

        for (pattern, hint) in legacy_checks_common
            matched = pattern isa Regex ? occursin(pattern, content) : occursin(pattern, content)
            if matched
                push!(violations, "$(display): 命中不一致模式 => $(hint)")
            end
        end
        for (pattern, hint) in strict_legacy_checks
            matched = pattern isa Regex ? occursin(pattern, content) : occursin(pattern, content)
            if matched
                push!(violations, "$(display): 命中不一致模式 => $(hint)")
            end
        end
    end

    for file in archived_files
        content = read(file, String)
        if occursin(misleading_entry_pattern, content)
            display = relpath_display(file, root)
            push!(violations, "$(display): archived 文档把旧路径表述为当前/主线入口，存在误导风险")
        end
    end

    if isfile(readme)
        readme_content = read(readme, String)
        if occursin(r"src/pnjl/PNJL\.jl", readme_content) || occursin(r"src/pnjl/", readme_content)
            push!(violations, "README.md: 命中不一致模式 => 旧主线路径，应改为 Models + src/models/entrypoints.jl")
        end
    end

    readme_content = isfile(readme) ? read(readme, String) : ""
    readme_marks_wip = occursin("修复中", readme_content)
    readme_marks_validated = occursin("截面/弛豫时间链路（已验证可用）", readme_content) ||
                             occursin("截面/弛豫时间链路已验证可用", readme_content)

    if readme_marks_wip
        absolute_state = r"系统完全就绪|已知问题\s*[:：]\s*无"
        for file in guide_files
            content = read(file, String)
            if occursin(absolute_state, content)
                display = relpath_display(file, root)
                push!(violations, "$(display): README 标注部分链路“修复中”，此文件却出现绝对化状态词")
            end
        end
    end

    if readme_marks_validated
        stale_wip = r"截面/弛豫时间.*修复中|截面/弛豫时间.*仍在修复"
        for file in guide_files
            content = read(file, String)
            if occursin(stale_wip, content)
                display = relpath_display(file, root)
                push!(violations, "$(display): README 已标注“已验证可用”，此文件仍保留“修复中”表述")
            end
        end
    end

    if !isempty(violations)
        println("[docs-consistency] FAILED: 发现 $(length(violations)) 项不一致")
        for item in violations
            println(" - " * item)
        end
        exit(1)
    end

    println("[docs-consistency] OK: 未发现历史路径或状态口径冲突")
end

main()
