# Retained-state input and provenance for controlled threshold algorithm studies.
if !isdefined(Main, :RelaxTime)
    include(joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "RelaxTime.jl"))
end

module ThresholdUpgradeCommon

using CSV, JSON, SHA, Dates
using Main: RelaxTime
using Main.Constants_PNJL: G_fm2, K_fm5

export ROOT, AS, CS, retained_row, retained_state, rate_config, write_json,
       write_rows, start_run, finish_run, raw_cache, evaluate_cache,
       local_integral, direct_integral, reconstruct_tau

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const AS = RelaxTime.AverageScatteringRate
const CS = RelaxTime.TotalCrossSection

function write_json(path, value)
    mkpath(dirname(path))
    open(path * ".writing", "w") do io
        JSON.print(io, value, 2)
        println(io)
    end
    mv(path * ".writing", path; force=true)
end

function write_rows(path, rows)
    isempty(rows) && return
    mkpath(dirname(path))
    CSV.write(path, rows)
end

function run_source_identity(root, paths; package_manifest="")
    if isempty(package_manifest)
        ispath(joinpath(root, ".git")) || throw(ArgumentError(
            "Frozen source requires THRESHOLD_SOURCE_MANIFEST pointing to its package_manifest.json"))
        return (branch=strip(read(Cmd(["git", "-C", root, "branch", "--show-current"]), String)),
            head=strip(read(Cmd(["git", "-C", root, "rev-parse", "HEAD"]), String)),
            source_archive=nothing)
    end

    isfile(package_manifest) || throw(ArgumentError("Missing source package manifest: $package_manifest"))
    manifest_bytes = read(package_manifest)
    package = JSON.parse(String(copy(manifest_bytes)))
    files = get(package, "source_files", nothing)
    files isa AbstractDict && !isempty(files) ||
        throw(ArgumentError("Source package manifest requires a nonempty source_files map"))
    archive = joinpath(dirname(abspath(package_manifest)), "source.zip")
    isfile(archive) || throw(ArgumentError("Missing frozen source archive: $archive"))
    archive_hash = bytes2hex(sha256(read(archive)))
    archive_hash == get(package, "zip_sha256", nothing) ||
        throw(ArgumentError("Frozen source archive SHA256 mismatch"))
    get(package, "uncommitted_source", true) && get(package, "calculation_commit", nothing) !== nothing &&
        throw(ArgumentError("Uncommitted source cannot declare a calculation_commit"))
    for (relative, expected) in files
        relative = replace(relative, '\\'=>'/')
        (isabspath(relative) || ".." in split(relative, '/')) &&
            throw(ArgumentError("Source package path must stay inside the source root: $relative"))
        path = joinpath(root, relative)
        isfile(path) || throw(ArgumentError("Missing frozen source file: $relative"))
        bytes2hex(sha256(read(path))) == expected ||
            throw(ArgumentError("Frozen source file SHA256 mismatch: $relative"))
    end
    for path in paths
        relative = replace(relpath(path, root), '\\'=>'/')
        haskey(files, relative) ||
            throw(ArgumentError("Run source is absent from frozen package manifest: $relative"))
    end
    return (branch=nothing, head=nothing, source_archive=(
        path=archive, sha256=archive_hash, package_manifest=abspath(package_manifest),
        package_manifest_sha256=bytes2hex(sha256(manifest_bytes)),
        package_manifest_snapshot="source_package_manifest.json", verified_files=length(files),
        head_base=get(package, "head_base", nothing), source_branch=get(package, "source_branch", nothing),
        calculation_commit=get(package, "calculation_commit", nothing),
        uncommitted_source=get(package, "uncommitted_source", true)))
end

function start_run(outdir, inputdir, script, parameters; extra_sources=String[],
    package_manifest=get(ENV, "THRESHOLD_SOURCE_MANIFEST", ""))
    isdir(outdir) && !isempty(readdir(outdir)) &&
        error("Use a new empty output directory: $outdir")
    mkpath(outdir)
    inputs = NamedTuple[]
    for name in ("review_points.csv", "channel_plan.csv")
        src = joinpath(inputdir, name)
        dst = joinpath(outdir, "inputs", name)
        mkpath(dirname(dst))
        cp(src, dst)
        push!(inputs, (path=abspath(src), sha256=bytes2hex(sha256(read(src))),
                       snapshot=relpath(dst, outdir)))
    end
    paths = String[]
    for base in ("src", "config")
        for (dir, _, files) in walkdir(joinpath(ROOT, base))
            for name in files
                endswith(name, ".jl") || endswith(name, ".toml") || continue
                push!(paths, joinpath(dir, name))
            end
        end
    end
    append!(paths, [joinpath(ROOT, "Project.toml"), joinpath(ROOT, "Manifest.toml"),
                    abspath(script), abspath(@__FILE__)])
    append!(paths, extra_sources)
    identity = run_source_identity(ROOT, paths; package_manifest=package_manifest)
    isempty(package_manifest) || cp(package_manifest, joinpath(outdir, "source_package_manifest.json"))
    sources = NamedTuple[]
    for src in unique(paths)
        relative = relpath(src, ROOT)
        startswith(relative, "..") && error("Source snapshot must be inside this worktree")
        dst = joinpath(outdir, "source", relative)
        mkpath(dirname(dst))
        cp(src, dst)
        push!(sources, (path=relative, sha256=bytes2hex(sha256(read(src)))))
    end
    write_json(joinpath(outdir, "manifest.json"), (;
        status="analysis_only", created_utc=string(now(UTC)), julia=string(VERSION),
        root=ROOT, identity...,
        input_files=inputs, source_files=sources, parameters=parameters,
        runtime_sigma_cache_policies=String.(AS.SIGMA_CACHE_POLICIES),
        numerical_defaults_promoted_by_run=false))
end

finish_run(outdir; kwargs...) = write_json(joinpath(outdir, "status.json"),
    (; complete=true, finished_utc=string(now(UTC)), kwargs...))

function retained_row(inputdir, xi)
    rows = [r for r in CSV.File(joinpath(inputdir, "review_points.csv"))
        if String(r.review_mode) == "mode_a" && String(r.plot_panel) == "muB0.0" &&
           String(r.plot_series) == "alpha1.0" && abs(Float64(r.xi)-xi) < 1e-10]
    length(rows) == 1 || error("Expected one retained muB=0 alphaT=1 row at xi=$xi")
    return only(rows)
end

function retained_state(row)
    qp0 = (m=(u=Float64(row.m_u), d=Float64(row.m_d), s=Float64(row.m_s)),
           μ=(u=Float64(row.muq_fm), d=Float64(row.muq_fm), s=Float64(row.muq_fm)))
    tp = (T=Float64(row.T_fm), Φ=Float64(row.Phi),
          Φbar=Float64(row.Phibar), ξ=Float64(row.xi))
    av = RelaxTime.AFieldBuilder.build_A_triplet(qp0, tp;
        p_nodes=16, p_max=20.0, cos_nodes=4, use_aniso=true)
    # The revised production convention uses the same full-RS A values for
    # polarization and effective couplings. No gap equation is re-solved.
    au = av.u
    ass = av.s
    gu = RelaxTime.EffectiveCouplings.calculate_G_from_A(au, qp0.m.u)
    gs = RelaxTime.EffectiveCouplings.calculate_G_from_A(ass, qp0.m.s)
    kc = RelaxTime.EffectiveCouplings.calculate_effective_couplings(G_fm2, K_fm5, gu, gs)
    return (quark_params=merge(qp0, (A=av,)), thermo_params=tp, K_coeffs=kc)
end

rate_config() = (
    p_nodes=128, angle_nodes=20, phi_nodes=36, n_sigma_points=24,
    sigma_grid_n=560, sigma_cutoff=Main.Constants_PNJL.Λ_inv_fm,
    sigma_cache_policy=:validated_anchored, threshold_subtraction=true,
    asym_window=0.6, asym_fit_min_points=8, asym_extra_points=10,
    interpolation_mode=:linear)

function raw_cache(process, grid, values, stock; subtract=false)
    cache = AS.CrossSectionCache(process, copy(grid), copy(values))
    if subtract
        cache.asym_enabled = stock.asym_enabled
        cache.asym_requested = true
        cache.asym_s0 = stock.asym_s0
        cache.asym_A = stock.asym_A
        cache.asym_taper_start = stock.asym_taper_start
        cache.asym_taper_end = stock.asym_taper_end
        for i in eachindex(grid)
            cache.sigma_vals[i] = AS._regularized_sigma_value(cache, grid[i], values[i])
        end
    end
    return cache
end

function evaluate_cache(cache, state, s)
    AS._get_sigma_core(cache, s, state.quark_params, state.thermo_params,
        state.K_coeffs; n_points=24, interpolation_mode=:linear)
end

function local_integral(cache, state; xmax=0.6, order=12)
    s0 = 4state.quark_params.m.u^2
    edges = sort!(unique(vcat([0.0, xmax],
        [s-s0 for s in cache.s_vals if 0 < s-s0 < xmax])))
    z, w = RelaxTime.GaussLegendre.gauleg(-1.0, 1.0, order)
    total = 0.0
    for (a,b) in zip(edges[1:end-1], edges[2:end])
        for (zi,wi) in zip(z,w)
            x = (a+b)/2 + (b-a)*zi/2
            total += (b-a)*wi/2*x*evaluate_cache(cache, state, s0+x)
        end
    end
    return total
end

function direct_integral(process, state; xmax=0.6, order=128, t_nodes=24)
    s0 = 4state.quark_params.m.u^2
    q, w = RelaxTime.GaussLegendre.gauleg(0.0, sqrt(xmax), order)
    sum(wi*2qi^3*CS.total_cross_section(process, s0+qi^2,
        state.quark_params, state.thermo_params, state.K_coeffs; n_points=t_nodes)
        for (qi,wi) in zip(q,w))
end

function reconstruct_tau(rows)
    result = NamedTuple[]
    for xi in unique(r.xi for r in rows), method in unique(r.method for r in rows)
        selected = filter(r -> r.xi == xi && r.method == method, rows)
        length(selected) == 2 || continue
        old = first(selected).source_total
        all(r -> isapprox(r.source_total, old; rtol=1e-12), selected) ||
            error("Different retained total inverse times in a two-channel pair")
        new = old + sum(r.source_contribution*(r.rate/r.source_rate-1) for r in selected)
        push!(result, (xi=xi, method=method, inverse_tau=new, tau=1/new,
                       source_tau=1/old, replaced_channels=2))
    end
    return result
end

end
