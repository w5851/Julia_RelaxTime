#!/usr/bin/env julia

"""
文档一致性检查（轻量门禁）

检查目标：
1) docs/guides/**/*.md 与 README.md 中是否含历史路径模式。
2) 当 README 标注“修复中”时，guides 不应出现绝对化状态词。

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
    guides_dir = joinpath(root, "docs", "guides")

    guide_files = collect_markdown_files(guides_dir)
    targets = guide_files

    violations = String[]

    legacy_checks = [
        ("test_unit/", "历史测试目录，应改为 tests/unit/"),
        (r"\bjulia\s+server\.jl\b", "历史服务启动命令，应改为 scripts/server/server*.jl"),
        (".\\start.bat", "历史根目录启动脚本，应改为 .\\scripts\\server\\start.bat"),
        ("doc/domain-knowledge/", "历史文档目录，应改为 docs/reference/domain-knowledge/"),
    ]

    for file in targets
        content = read(file, String)
        display = relpath_display(file, root)

        for (pattern, hint) in legacy_checks
            matched = pattern isa Regex ? occursin(pattern, content) : occursin(pattern, content)
            if matched
                push!(violations, "$(display): 命中不一致模式 => $(hint)")
            end
        end
    end

    readme_content = isfile(readme) ? read(readme, String) : ""
    readme_marks_wip = occursin("修复中", readme_content)

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
