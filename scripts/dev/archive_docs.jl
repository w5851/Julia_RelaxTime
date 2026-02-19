#!/usr/bin/env julia
# Archive development documents from docs/dev/active to docs/dev/archived
# Adds metadata header and renames with date prefix

using Dates
using Printf

const ROOT = pwd()
const ACTIVE_DIR = joinpath(ROOT, "docs", "dev", "active")
const ARCHIVED_DIR = joinpath(ROOT, "docs", "dev", "archived")

# Ensure archived directory exists
mkpath(ARCHIVED_DIR)

"""
    extract_title(content::String) -> String

Extract title from markdown content (first # heading or filename-based fallback)
"""
function extract_title(content::String, filename::String)
    # Try to find first # heading
    for line in split(content, '\n')
        m = match(r"^#\s+(.+)$", strip(line))
        if m !== nothing
            return strip(m.captures[1])
        end
    end
    # Fallback: use filename without extension and date prefix
    base = replace(filename, r"^\d{4}[-_]\d{2}[-_]\d{2}[-_]" => "")
    base = replace(base, r"\.md$" => "")
    return replace(base, "_" => " ")
end

"""
    generate_metadata(title::String, original_path::String, archive_date::String) -> String

Generate YAML frontmatter for archived document
"""
function generate_metadata(title::String, original_path::String, archive_date::String)
    return """---
title: $title
archived: true
original: $original_path
archived_date: $archive_date
---

"""
end

"""
    archive_file(filepath::String; date::String=today_str(), dry_run::Bool=false) -> Bool

Archive a single file from active to archived directory
"""
function archive_file(filepath::String; date::String=Dates.format(today(), "yyyy-mm-dd"), dry_run::Bool=false)
    if !isfile(filepath)
        @error "File not found: $filepath"
        return false
    end
    
    # Read original content
    content = read(filepath, String)
    filename = basename(filepath)
    
    # Extract title
    title = String(extract_title(content, filename))
    
    # Generate relative path for metadata
    rel_path = "docs/dev/active/$filename"
    
    # Generate metadata
    metadata = generate_metadata(title, rel_path, date)
    
    # Generate archived filename and normalize to YYYY-MM-DD_...
    if !occursin(r"^\d{4}[-_]\d{2}[-_]\d{2}[-_]", filename)
        archived_filename = "$(date)_$filename"
    else
        m = match(r"^(\d{4})[-_](\d{2})[-_](\d{2})[-_](.+)$", filename)
        if m === nothing
            archived_filename = filename
        else
            yyyy, mm, dd, rest = m.captures
            archived_filename = "$(yyyy)-$(mm)-$(dd)_$(rest)"
        end
    end
    
    archived_path = joinpath(ARCHIVED_DIR, archived_filename)
    
    # Check if archived file already exists
    if isfile(archived_path) && !dry_run
        print("File already exists: $archived_path. Overwrite? (y/N): ")
        response = readline()
        if lowercase(strip(response)) != "y"
            println("Skipped: $filename")
            return false
        end
    end
    
    # Prepare archived content
    archived_content = metadata * "\n以下为原始内容（保留，以便审阅与历史参考）：\n\n---\n\n" * content
    
    if dry_run
        println("\n[DRY RUN] Would archive:")
        println("  From: $filepath")
        println("  To:   $archived_path")
        println("  Title: $title")
        return true
    end
    
    # Write archived file
    open(archived_path, "w") do io
        write(io, archived_content)
    end
    
    # Delete original file
    rm(filepath)
    
    println("✓ Archived: $filename -> $archived_filename")
    return true
end

"""
    list_active_files() -> Vector{String}

List all markdown files in active directory
"""
function list_active_files()
    if !isdir(ACTIVE_DIR)
        return String[]
    end
    files = String[]
    for f in readdir(ACTIVE_DIR)
        if endswith(f, ".md")
            push!(files, joinpath(ACTIVE_DIR, f))
        end
    end
    return sort(files)
end

"""
    check_archived_format(filepath::String) -> Bool

Check if an archived file has proper metadata format
"""
function check_archived_format(filepath::String)
    content = read(filepath, String)
    
    # Check for YAML frontmatter
    if match(r"^---\r?\n", content) === nothing
        return false
    end
    
    # Check for required fields
    required_fields = ["title:", "archived:", "original:", "archived_date:"]
    for field in required_fields
        if !occursin(field, content)
            return false
        end
    end
    
    return true
end

"""
    validate_archived_files()

Validate all files in archived directory
"""
function validate_archived_files()
    if !isdir(ARCHIVED_DIR)
        println("Archived directory not found: $ARCHIVED_DIR")
        return
    end
    
    files = [f for f in readdir(ARCHIVED_DIR) if endswith(f, ".md")]
    
    if isempty(files)
        println("No archived files found.")
        return
    end
    
    println("Validating $(length(files)) archived files...\n")
    
    valid_count = 0
    invalid_files = String[]
    
    for f in files
        filepath = joinpath(ARCHIVED_DIR, f)
        if check_archived_format(filepath)
            valid_count += 1
        else
            push!(invalid_files, f)
        end
    end
    
    println("Valid: $valid_count / $(length(files))")
    
    if !isempty(invalid_files)
        println("\nInvalid files (missing metadata):")
        for f in invalid_files
            println("  - $f")
        end
    end
end

"""
    interactive_mode()

Interactive mode for archiving files
"""
function interactive_mode()
    files = list_active_files()
    
    if isempty(files)
        println("No markdown files found in $ACTIVE_DIR")
        return
    end
    
    println("Active files:")
    for (i, f) in enumerate(files)
        println("  $i. $(basename(f))")
    end
    
    print("\nEnter file numbers to archive (comma-separated, or 'all'): ")
    input = readline()
    
    if lowercase(strip(input)) == "all"
        indices = 1:length(files)
    else
        try
            indices = [parse(Int, strip(s)) for s in split(input, ',')]
        catch
            println("Invalid input")
            return
        end
    end
    
    print("Archive date (YYYY-MM-DD, default: today): ")
    date_input = strip(readline())
    archive_date = isempty(date_input) ? Dates.format(today(), "yyyy-mm-dd") : date_input
    
    println("\nArchiving $(length(indices)) file(s)...\n")
    
    success_count = 0
    for i in indices
        if i < 1 || i > length(files)
            println("Skipping invalid index: $i")
            continue
        end
        if archive_file(files[i]; date=archive_date)
            success_count += 1
        end
    end
    
    println("\n✓ Successfully archived $success_count file(s)")
end

"""
    print_usage()

Print usage information
"""
function print_usage()
    println("""
Usage: julia scripts/dev/archive_docs.jl [OPTIONS] [FILES...]

Archive development documents from docs/dev/active to docs/dev/archived.

Options:
  -h, --help              Show this help message
  -i, --interactive       Interactive mode (select files to archive)
  -c, --check             Validate archived files format
  -d, --date DATE         Archive date (YYYY-MM-DD, default: today)
  --dry-run               Show what would be done without making changes
  -b, --batch FILES...    Batch archive specified files

Examples:
  # Interactive mode
  julia scripts/dev/archive_docs.jl -i

  # Batch archive specific files
  julia scripts/dev/archive_docs.jl file1.md file2.md

  # Batch archive with custom date
  julia scripts/dev/archive_docs.jl -d 2026-01-15 file1.md

  # Validate archived files
  julia scripts/dev/archive_docs.jl -c

  # Dry run
  julia scripts/dev/archive_docs.jl --dry-run file1.md
""")
end

# Main execution
function main()
    args = ARGS
    
    if isempty(args) || "-h" in args || "--help" in args
        print_usage()
        return
    end
    
    if "-c" in args || "--check" in args
        validate_archived_files()
        return
    end
    
    if "-i" in args || "--interactive" in args
        interactive_mode()
        return
    end
    
    # Parse options
    archive_date = Dates.format(today(), "yyyy-mm-dd")
    dry_run = false
    files_to_archive = String[]
    
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "-d" || arg == "--date"
            if i + 1 <= length(args)
                archive_date = args[i + 1]
                i += 2
            else
                println("Error: --date requires a value")
                return
            end
        elseif arg == "--dry-run"
            dry_run = true
            i += 1
        elseif arg == "-b" || arg == "--batch"
            # Remaining args are files
            files_to_archive = args[i+1:end]
            break
        elseif !startswith(arg, "-")
            push!(files_to_archive, arg)
            i += 1
        else
            println("Unknown option: $arg")
            print_usage()
            return
        end
    end
    
    if isempty(files_to_archive)
        println("No files specified. Use -i for interactive mode or specify files.")
        print_usage()
        return
    end
    
    # Archive specified files
    success_count = 0
    for filename in files_to_archive
        # Handle both full paths and just filenames
        filepath = if isfile(filename)
            filename
        elseif isfile(joinpath(ACTIVE_DIR, filename))
            joinpath(ACTIVE_DIR, filename)
        else
            println("File not found: $filename")
            continue
        end
        
        if archive_file(filepath; date=archive_date, dry_run=dry_run)
            success_count += 1
        end
    end
    
    if !dry_run
        println("\n✓ Successfully archived $success_count file(s)")
    end
end

# Run main if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
