module PNJLBenchmarkMetadata

using SHA
using LinearAlgebra

function file_digest(path::String)
    return isfile(path) ? bytes2hex(sha256(read(path))) : nothing
end

function config_digest(root::String)
    files = String[]
    for subdir in ("config/models/pnjl", "config/physics")
        for (dir, _, names) in walkdir(joinpath(root, subdir))
            append!(files, [joinpath(dir, name) for name in names if endswith(name, ".toml")])
        end
    end
    records = [replace(relpath(path, root), '\\' => '/') * ":" * file_digest(path) for path in sort(files)]
    return bytes2hex(sha256(join(records, "\n")))
end

function metadata(root::String)
    git(args...) = readchomp(Cmd(["git", "-C", root, args...]))
    source_paths = ["src", "config/models", "config/physics", "Project.toml", "Manifest.toml",
                    "benchmark/pnjl", "scripts/perf/pnjl"]
    dirty = !isempty(git("status", "--porcelain", "--untracked-files=normal", "--", source_paths...))
    workload_paths = ["benchmark/pnjl/single_point_solver_perf.jl",
                      "scripts/perf/pnjl/scan_perf.jl", "benchmark/pnjl/benchmark_metadata.jl"]
    workload_hash = bytes2hex(sha256(join([file_digest(joinpath(root, path)) for path in workload_paths], "\n")))
    cpu_info = Sys.cpu_info()
    return (
        provenance=(
            commit=git("rev-parse", "HEAD"),
            branch=get(ENV, "GITHUB_HEAD_REF", "") == "" ?
                   get(ENV, "GITHUB_REF_NAME", git("branch", "--show-current")) : ENV["GITHUB_HEAD_REF"],
            source_dirty=dirty,
            repository=get(ENV, "GITHUB_REPOSITORY", ""),
            run_id=get(ENV, "GITHUB_RUN_ID", ""),
            run_attempt=get(ENV, "GITHUB_RUN_ATTEMPT", ""),
            event=get(ENV, "GITHUB_EVENT_NAME", "local"),
        ),
        environment=(
            julia_version=string(VERSION),
            julia_threads=Threads.nthreads(),
            blas_threads=BLAS.get_num_threads(),
            blas_config=string(BLAS.get_config()),
            os=string(Sys.KERNEL),
            arch=string(Sys.ARCH),
            cpu_name=Sys.CPU_NAME,
            cpu_model=isempty(cpu_info) ? "unknown" : first(cpu_info).model,
            cpu_threads=Sys.CPU_THREADS,
            runner_os=get(ENV, "RUNNER_OS", "local"),
            runner_arch=get(ENV, "RUNNER_ARCH", string(Sys.ARCH)),
            runner_image=get(ENV, "ImageOS", "local"),
            runner_image_version=get(ENV, "ImageVersion", ""),
            pnjl_profile=get(ENV, "PNJL_PARAM_PROFILE", "default"),
            physics_profile=get(ENV, "PHYSICS_PARAM_PROFILE", "default"),
            project_sha256=file_digest(joinpath(root, "Project.toml")),
            manifest_sha256=file_digest(joinpath(root, "Manifest.toml")),
            config_sha256=config_digest(root),
            workload_sha256=workload_hash,
        ),
    )
end

end # module
