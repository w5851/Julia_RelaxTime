#!/usr/bin/env julia

using Dates

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))

function _read_csv_sections(path::String)
    lines = readlines(path)
    meta = String[]
    header = ""
    rows = String[]
    for line in lines
        s = strip(line)
        isempty(s) && continue
        if startswith(s, "#") && isempty(header)
            push!(meta, line)
        elseif isempty(header)
            header = line
        else
            push!(rows, line)
        end
    end
    isempty(header) && throw(ArgumentError("missing header: $path"))
    return meta, header, rows
end

function _write_merged_csv(out_path::String, meta::Vector{String}, header::String, rows::Vector{String}, source_tags::Vector{String})
    mkpath(dirname(out_path))
    open(out_path, "w") do io
        for m in meta
            println(io, m)
        end
        println(io, "# merged_from=", join(source_tags, ","))
        println(io, "# merged_at_utc=", Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SSZ"))
        println(io, header)
        for r in rows
            println(io, r)
        end
    end
end

function _copy_prefixed_figures(out_dir::String, source_dirs::Vector{String})
    for sub in ("fig_gamma", "fig_mass")
        mkpath(joinpath(out_dir, sub))
    end
    for src in source_dirs
        tag = basename(src)
        for sub in ("fig_gamma", "fig_mass")
            fig_dir = joinpath(src, sub)
            isdir(fig_dir) || continue
            for f in readdir(fig_dir)
                src_path = joinpath(fig_dir, f)
                dst_path = joinpath(out_dir, sub, tag * "__" * f)
                cp(src_path, dst_path; force=true)
            end
        end
    end
end

function _write_readme(path::String, source_dirs::Vector{String})
    open(path, "w") do io
        println(io, "# mixed_window_globalmin_xi-0p3_muB0_merged")
        println(io)
        println(io, "Merged temporary experiment outputs for eta/eta_prime window diagnostics.")
        println(io)
        println(io, "## Scope")
        println(io, "- xi = -0.3")
        println(io, "- muB = 0")
        println(io, "- T window = 220..280 MeV")
        println(io)
        println(io, "## Sources")
        for d in source_dirs
            println(io, "- ", basename(d))
        end
        println(io)
        println(io, "## Files")
        println(io, "- window_eta_etap.csv: merged pointwise outputs (Gamma/mass, residual, selected method)")
        println(io, "- window_jump_summary.csv: jump-count summary under jump-threshold=0.25")
        println(io, "- fig_gamma/, fig_mass/: copied figures with source-prefix in filenames")
        println(io)
        println(io, "## Purpose")
        println(io, "This dataset compares optimizer-dependent jump behavior for eta/eta_prime curves and supports follow-up method selection.")
    end
end

function main()
    root = joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "mott_phase")
    src1 = joinpath(root, "mixed_window_globalmin_xi-0p3_muB0")
    src2 = joinpath(root, "mixed_window_globalmin_xi-0p3_muB0_next2")
    out = joinpath(root, "mixed_window_globalmin_xi-0p3_muB0_merged")

    sources = [src1, src2]
    all(isdir, sources) || throw(ArgumentError("source directories missing under $root"))

    meta1, header1, rows1 = _read_csv_sections(joinpath(src1, "window_eta_etap.csv"))
    _meta2, header2, rows2 = _read_csv_sections(joinpath(src2, "window_eta_etap.csv"))
    header1 == header2 || throw(ArgumentError("window_eta_etap header mismatch"))

    _write_merged_csv(
        joinpath(out, "window_eta_etap.csv"),
        meta1,
        header1,
        vcat(rows1, rows2),
        basename.(sources),
    )

    _jmeta1, jheader1, jrows1 = _read_csv_sections(joinpath(src1, "window_jump_summary.csv"))
    _jmeta2, jheader2, jrows2 = _read_csv_sections(joinpath(src2, "window_jump_summary.csv"))
    jheader1 == jheader2 || throw(ArgumentError("window_jump_summary header mismatch"))

    open(joinpath(out, "window_jump_summary.csv"), "w") do io
        println(io, jheader1)
        for r in vcat(jrows1, jrows2)
            println(io, r)
        end
    end

    _copy_prefixed_figures(out, sources)
    _write_readme(joinpath(out, "README.md"), sources)

    println("Merged outputs written to: ", out)
end

main()
