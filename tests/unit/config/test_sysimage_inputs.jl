using Test

if !isdefined(Main, :SysimageInputs)
    include(joinpath(@__DIR__, "..", "..", "..", "scripts", "dev", "sysimage_inputs.jl"))
end

@testset "sysimage freshness follows build inputs" begin
    mktempdir() do root
        for rel in SysimageInputs.REQUIRED_FILES
            path = joinpath(root, rel)
            mkpath(dirname(path))
            write(path, "fixture\n")
        end
        mkpath(joinpath(root, "src"))
        write(joinpath(root, "src", "fixture.jl"), "answer() = 42\n")
        original = SysimageInputs.fingerprint(root)
        @test length(original) == 64
        @test SysimageInputs.fingerprint(root) == original
        mkpath(joinpath(root, "docs"))
        write(joinpath(root, "docs", "README.md"), "new documentation\n")
        mkpath(joinpath(root, "config", "governance"))
        write(joinpath(root, "config", "governance", "task_tracks.toml"), "updated_at = '2026-10-02'\n")
        @test SysimageInputs.fingerprint(root) == original

        # Every build input category affects freshness, even without a commit.
        for rel in ("src/fixture.jl", "src/added.jl", "config/models/pnjl/default.toml",
                    "config/physics/constants.toml", "Project.toml", "Manifest.toml",
                    "scripts/dev/precompile_workload.jl", "scripts/dev/build_sysimage.jl",
                    "scripts/dev/sysimage_inputs.jl")
            path = joinpath(root, rel)
            mkpath(dirname(path))
            previous = isfile(path) ? read(path) : nothing
            write(path, "changed\n")
            @test SysimageInputs.fingerprint(root) != original
            previous === nothing ? rm(path) : write(path, previous)
            @test SysimageInputs.fingerprint(root) == original
        end
        mv(joinpath(root, "src", "fixture.jl"), joinpath(root, "src", "renamed.jl"))
        @test SysimageInputs.fingerprint(root) != original
        rm(joinpath(root, "src", "renamed.jl"))
        @test SysimageInputs.fingerprint(root) != original
        rm(joinpath(root, "Project.toml"))
        @test_throws ArgumentError SysimageInputs.fingerprint(root)
    end
end
