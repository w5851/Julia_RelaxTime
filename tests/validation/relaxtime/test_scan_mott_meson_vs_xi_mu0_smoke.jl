using Test

if !isdefined(Main, :validation_targets_path)
    include(joinpath(@__DIR__, "..", "common", "data_paths.jl"))
end

const _SCAN_SCRIPT = joinpath(
    @__DIR__,
    "..", "..", "..",
    "scripts", "analysis", "scan_mott_meson_vs_xi_mu0.jl",
)

@testset "scan_mott_meson_vs_xi_mu0 smoke" begin
    @test isfile(_SCAN_SCRIPT)

    out = joinpath(
        dirname(_SCAN_SCRIPT),
        "..", "..", "data", "outputs", "results", "relaxtime", "scan",
        "mott_meson_vs_xi_mu0_smoke.csv",
    )
    out = normpath(out)

    if isfile(out)
        rm(out)
    end

    cmd = `julia --project=. $(abspath(_SCAN_SCRIPT)) --smoke --output $(out)`
    run(cmd)

    @test isfile(out)
    lines = readlines(out)
    @test length(lines) >= 2
    header = nothing
    for ln in lines
        s = strip(ln)
        isempty(s) && continue
        startswith(s, '#') && continue
        header = s
        break
    end
    @test header !== nothing
    @test occursin("xi", header)
    @test occursin("meson", header)
    @test occursin("T_Mott_MeV", header)
end
