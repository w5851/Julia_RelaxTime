using Test

if !isdefined(Main, :ArchiveDocs)
    include(joinpath(@__DIR__, "..", "..", "..", "scripts", "dev", "archive_docs.jl"))
end
const AD = Main.ArchiveDocs

@testset "archive previews and help are read-only" begin
    mktempdir() do root
        @test AD.main(["--help"]; root, io=devnull) == 0
        @test AD.main(["--check"]; root, io=devnull) == 0
        @test isempty(readdir(root))
        active = joinpath(root, "docs", "dev", "active")
        mkpath(active)
        source = joinpath(active, "example.md")
        original = "# 测试: \"标题\" # 符号\r\n\r\n原文\r\n"
        write(source, original)
        @test AD.main(["--dry-run", "--date", "2026-10-07", "example.md"]; root, io=devnull) == 0
        @test read(source, String) == original
        @test !ispath(joinpath(root, "docs", "dev", "archived"))
        @test AD.main(["--date", "2026-10-07", "example.md"]; root, io=devnull) == 0
        target = joinpath(root, "docs", "dev", "archived", "2026-10-07_example.md")
        @test !ispath(source)
        @test endswith(read(target, String), original)
        fields = AD.metadata_fields(read(target, String))
        @test fields["title"] == "测试: \"标题\" # 符号"
        @test fields["original"] == "docs/dev/active/example.md"
        @test fields["task_status"] == "completed"
        @test AD.main(["--check", basename(target)]; root, io=devnull) == 0
        write(joinpath(dirname(target), "bad.md"), "---\n---\narchived: true\ntitle: x\noriginal: x\narchived_date: x\n")
        @test AD.main(["--check", basename(target)]; root, io=devnull) == 0
        @test AD.main(["--check"]; root, io=devnull) == 1
    end
end

@testset "archive rejects unsafe paths, conflicts and incomplete batches" begin
    mktempdir() do root
        active = joinpath(root, "docs", "dev", "active")
        archived = joinpath(root, "docs", "dev", "archived")
        mkpath(active)
        source = joinpath(active, "task.md")
        write(source, "# Task\n")
        outside = joinpath(root, "outside.md")
        write(outside, "must remain\n")
        @test AD.main([outside]; root, io=devnull) == 1
        @test read(outside, String) == "must remain\n"
        @test AD.main(["task.md", "missing.md"]; root, io=devnull) == 1
        @test AD.main(["task.md", "task.md"]; root, io=devnull) == 1
        if Sys.iswindows()
            @test AD.main(["task.md", "TASK.MD"]; root, io=devnull) == 1
        end
        @test AD.main(["--date", "2026-02-30", "task.md"]; root, io=devnull) == 1
        @test isfile(source)
        @test !ispath(archived)
        plan = AD.prepare_archive("task.md"; root, date="2026-10-07")
        write(source, "# Changed after preflight\n")
        @test_throws ArgumentError AD.execute_archive(plan)
        @test !ispath(plan.target)
        mkpath(archived)
        write(plan.target, "existing archive")
        @test AD.main(["--dry-run", "--date", "2026-10-07", "task.md"]; root, io=devnull) == 1
        @test AD.main(["--date", "2026-10-07", "task.md"]; root, io=devnull) == 1
        @test read(plan.target, String) == "existing archive"
        @test isfile(source)
    end
end

@testset "archive rejects directory redirection after preview" begin
    mktempdir() do root
        active = joinpath(root, "docs", "dev", "active")
        archived = joinpath(root, "docs", "dev", "archived")
        external = joinpath(root, "external")
        mkpath(active)
        mkpath(external)
        source = joinpath(active, "task.md")
        write(source, "# Task\n")
        plan = AD.prepare_archive("task.md"; root, date="2026-10-07")
        symlink(external, archived; dir_target=true)
        @test_throws ArgumentError AD.execute_archive(plan)
        @test AD.main(["--dry-run", "task.md"]; root, io=devnull) == 1
        @test read(source, String) == "# Task\n"
        @test isempty(readdir(external))
    end
end

@testset "archive terminal reasons and scalar validation" begin
    mktempdir() do root
        active = joinpath(root, "docs", "dev", "active")
        mkpath(active)
        write(joinpath(active, "task.md"), "# Task\n")
        @test AD.main(["--status", "cancelled", "task.md"]; root, io=devnull) == 1
        @test AD.main(["--status", "cancelled", "--reason", "Scope withdrawn", "--date", "2026-10-07", "task.md"]; root, io=devnull) == 0
        path = joinpath(root, "docs", "dev", "archived", "2026-10-07_task.md")
        valid = read(path, String)
        @test AD.metadata_fields(valid)["archive_reason"] == "Scope withdrawn"
        for invalid in (
            replace(valid, "archived: true" => "archived: \"true\""),
            replace(valid, "archived: true" => "archived: false"),
            replace(valid, "archived: true" => "archived: true\narchived: true"),
            replace(valid, "2026-10-07" => "2026-02-30"),
            replace(valid, "title: \"Task\"" => "title: unquoted: invalid"),
            replace(valid, "title: \"Task\"" => "title: null"),
            replace(valid, "title: \"Task\"" => "title: 123"),
            replace(valid, "title: \"Task\"" => "title: true"),
        )
            write(path, invalid)
            @test !AD.check_archived_format(path)
            @test AD.main(["--check", basename(path)]; root, io=devnull) == 1
        end
    end
    script = joinpath(@__DIR__, "..", "..", "..", "scripts", "dev", "archive_docs.jl")
    project = normpath(joinpath(@__DIR__, "..", "..", ".."))
    process = run(pipeline(ignorestatus(`$(Base.julia_cmd()) --startup-file=no --project=$project $script --unknown`); stdout=devnull, stderr=devnull))
    @test process.exitcode == 1
end
