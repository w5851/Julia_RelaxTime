#!/usr/bin/env julia

using Dates

const ROOT = pwd()
const ARCHIVED_DIR = joinpath(ROOT, "docs", "dev", "archived")

function has_valid_frontmatter(content::String)
    startswith(content, "---\n") || return false
    required_fields = ["title:", "archived:", "original:", "archived_date:"]
    for field in required_fields
        occursin(field, content) || return false
    end
    return true
end

function extract_title(content::String, filename::String)
    for line in split(content, '\n')
        m = match(r"^#\s+(.+)$", strip(line))
        if m !== nothing
            return strip(m.captures[1])
        end
    end
    base = replace(filename, r"\.md$" => "")
    base = replace(base, r"^\d{4}[-_]\d{2}[-_]\d{2}_" => "")
    return replace(base, "_" => " ")
end

function infer_archived_date(filename::String)
    m1 = match(r"^(\d{4})-(\d{2})-(\d{2})_", filename)
    if m1 !== nothing
        return "$(m1.captures[1])-$(m1.captures[2])-$(m1.captures[3])"
    end
    m2 = match(r"^(\d{4})_(\d{2})_(\d{2})_", filename)
    if m2 !== nothing
        return "$(m2.captures[1])-$(m2.captures[2])-$(m2.captures[3])"
    end
    return Dates.format(today(), "yyyy-mm-dd")
end

function infer_original_path(filename::String)
    base = replace(filename, r"^\d{4}[-_]\d{2}[-_]\d{2}_" => "")
    return "docs/dev/active/$(base)"
end

function wrap_with_frontmatter(content::String, filename::String)
    title = extract_title(content, filename)
    archived_date = infer_archived_date(filename)
    original = infer_original_path(filename)
    metadata = """
---
title: $title
archived: true
original: $original
archived_date: $archived_date
---

"""
    return metadata * "\n以下为原始内容（保留，以便审阅与历史参考）：\n\n---\n\n" * content
end

function main()
    isdir(ARCHIVED_DIR) || error("Archived directory not found: $(ARCHIVED_DIR)")
    files = sort([f for f in readdir(ARCHIVED_DIR) if endswith(f, ".md")])

    fixed = String[]
    skipped = String[]

    for file in files
        path = joinpath(ARCHIVED_DIR, file)
        content = read(path, String)
        if has_valid_frontmatter(content)
            push!(skipped, file)
            continue
        end

        new_content = wrap_with_frontmatter(content, file)
        open(path, "w") do io
            write(io, new_content)
        end
        push!(fixed, file)
    end

    println("[backfill-archived-frontmatter] fixed=$(length(fixed)), skipped=$(length(skipped))")
    if !isempty(fixed)
        println("Fixed files:")
        for file in fixed
            println("  - $(file)")
        end
    end
end

main()
