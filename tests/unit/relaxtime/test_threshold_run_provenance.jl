using Test, JSON, SHA
if !isdefined(Main, :ThresholdUpgradeCommon)
    include(joinpath(@__DIR__, "..", "..", "..", "scripts", "analysis", "relaxtime", "ThresholdUpgradeCommon.jl"))
end

@testset "Frozen analysis source identity" begin
    mktempdir() do directory
        root = joinpath(directory, "project")
        mkpath(joinpath(root, "src"))
        source = joinpath(root, "src", "solver.jl")
        write(source, "identity_fixture = 1\n")
        write(joinpath(root, "Project.toml"), "[deps]\n")
        archive = joinpath(directory, "source.zip")
        # Archive extraction is the packager's responsibility; the entrypoint
        # binds retained archive bytes and the extracted source file hashes.
        write(archive, "retained archive bytes")
        package = Dict{String,Any}(
            "zip_sha256"=>bytes2hex(sha256(read(archive))),
            "source_files"=>Dict("src/solver.jl"=>bytes2hex(sha256(read(source))),
                "Project.toml"=>bytes2hex(sha256(read(joinpath(root, "Project.toml"))))),
            "head_base"=>repeat("a", 40), "source_branch"=>"fixture-branch",
            "calculation_commit"=>nothing, "uncommitted_source"=>true)
        manifest = joinpath(directory, "package_manifest.json")
        save_manifest() = write(manifest, JSON.json(package))
        save_manifest()
        identify() = ThresholdUpgradeCommon.run_source_identity(root, [source]; package_manifest=manifest)
        identity = identify()
        @test identity.head === nothing
        @test identity.branch === nothing
        @test identity.source_archive.calculation_commit === nothing
        @test identity.source_archive.head_base == package["head_base"]
        @test identity.source_archive.sha256 == package["zip_sha256"]
        @test identity.source_archive.package_manifest_sha256 == bytes2hex(sha256(read(manifest)))
        @test identity.source_archive.verified_files == 2
        @test_throws ArgumentError ThresholdUpgradeCommon.run_source_identity(root, [source])

        write(source, "identity_fixture = 2\n")
        @test_throws ArgumentError identify()
        write(source, "identity_fixture = 1\n")
        write(archive, "changed archive bytes")
        @test_throws ArgumentError identify()
        write(archive, "retained archive bytes")

        delete!(package["source_files"], "src/solver.jl")
        save_manifest()
        @test_throws ArgumentError identify()
        package["source_files"]["src/solver.jl"] = bytes2hex(sha256(read(source)))
        package["source_files"]["src/missing.jl"] = repeat("0", 64)
        save_manifest()
        @test_throws ArgumentError identify()
        delete!(package["source_files"], "src/missing.jl")
        package["calculation_commit"] = package["head_base"]
        save_manifest()
        @test_throws ArgumentError identify()
        package["calculation_commit"] = nothing
        package["source_files"]["../outside.jl"] = repeat("0", 64)
        save_manifest()
        @test_throws ArgumentError identify()
    end
end
