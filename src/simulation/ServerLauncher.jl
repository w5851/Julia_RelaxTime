module ServerLauncher

using HTTP

include(joinpath(@__DIR__, "FullServerApp.jl"))
using .FullServerApp
include(joinpath(@__DIR__, "ServerWarmup.jl"))
using .ServerWarmup

export run_full_server, parse_port, server_runtime_policy, runtime_policy_env
export list_server_warmup_steps, resolve_server_warmup_profile
export run_server_warmup, run_server_warmup_from_env

const DEFAULT_PORT = 8080
const DEFAULT_DEPLOY_PROFILE = "localhost"

function _normalize_profile(profile::AbstractString)
    normalized = lowercase(strip(String(profile)))
    if normalized in ("localhost", "staging", "remote")
        return normalized
    end
    @warn "Unknown deploy profile; fallback to localhost" profile=profile
    return DEFAULT_DEPLOY_PROFILE
end

function _policy_for_profile(profile::String)
    if profile == "localhost"
        return (
            profile="localhost",
            cors_allow_origins="*",
            scan_max_running=2,
            scan_max_pending=32,
            scan_job_timeout_seconds=0,
            warmup_profile="none",
            warmup_strict=false,
        )
    elseif profile == "staging"
        return (
            profile="staging",
            cors_allow_origins="https://staging.jrt.local",
            scan_max_running=2,
            scan_max_pending=64,
            scan_job_timeout_seconds=180,
            warmup_profile="point",
            warmup_strict=false,
        )
    else
        return (
            profile="remote",
            cors_allow_origins="https://api.jrt.example.com",
            scan_max_running=4,
            scan_max_pending=128,
            scan_job_timeout_seconds=300,
            warmup_profile="service_core",
            warmup_strict=false,
        )
    end
end

function server_runtime_policy(profile::AbstractString=get(ENV, "JRT_DEPLOY_PROFILE", DEFAULT_DEPLOY_PROFILE))
    normalized = _normalize_profile(profile)
    return _policy_for_profile(normalized)
end

function runtime_policy_env(profile::AbstractString=get(ENV, "JRT_DEPLOY_PROFILE", DEFAULT_DEPLOY_PROFILE))
    policy = server_runtime_policy(profile)
    return Dict(
        "JRT_DEPLOY_PROFILE" => policy.profile,
        "JRT_CORS_ALLOW_ORIGINS" => policy.cors_allow_origins,
        "PNJL_SCAN_MAX_RUNNING" => string(policy.scan_max_running),
        "PNJL_SCAN_MAX_PENDING" => string(policy.scan_max_pending),
        "PNJL_SCAN_JOB_TIMEOUT_SECONDS" => string(policy.scan_job_timeout_seconds),
        "JRT_SERVER_WARMUP_PROFILE" => policy.warmup_profile,
        "JRT_SERVER_WARMUP_STRICT" => policy.warmup_strict ? "1" : "0",
    )
end

function parse_port(args::Vector{String})
    port = DEFAULT_PORT
    if length(args) >= 1
        try
            port = parse(Int, args[1])
            if port < 1024 || port > 65535
                @warn "端口号应在1024-65535之间，使用默认端口8080"
                port = DEFAULT_PORT
            end
        catch
            @warn "无效的端口号，使用默认端口8080"
        end
    end
    return port
end

function print_banner(port::Int, policy)
    println("\n" * "="^60)
    println("🚀 散射计算服务器启动中...")
    println("="^60)
    println("📍 服务地址: http://localhost:$port")
    println("📡 API端点:")
    println("   POST http://localhost:$port/compute")
    println("   GET  http://localhost:$port/health")
    println("   GET  http://localhost:$port/api/modules")
    println("   POST http://localhost:$port/api/modules/pnjl-gap/run")
    println("   POST http://localhost:$port/api/modules/pnjl-scan/jobs")
    println("   GET  http://localhost:$port/api/modules/pnjl-scan/jobs/{job_id}")
    println("   GET  http://localhost:$port/api/modules/pnjl-scan/jobs/{job_id}/result")
    println("\n⚙️  部署策略:")
    println("   profile: $(policy.profile)")
    println("   cors_allow_origins: $(policy.cors_allow_origins)")
    println("   scan_max_running: $(policy.scan_max_running)")
    println("   scan_max_pending: $(policy.scan_max_pending)")
    println("   scan_job_timeout_seconds: $(policy.scan_job_timeout_seconds)")
    println("   warmup_profile: $(policy.warmup_profile)")
    println("   warmup_strict: $(policy.warmup_strict)")
    println("\n📁 静态文件:")
    println("   http://localhost:$port/")
    println("   http://localhost:$port/index.html")
    println("   http://localhost:$port/simple_test.html")
    println("\n💡 提示:")
    println("   • 在浏览器中打开: http://localhost:$port")
    println("   • 按 Ctrl+C 停止服务器")
    println("="^60 * "\n")
end

function run_full_server(repo_root::String, args::Vector{String}=String[])
    port = parse_port(args)
    policy = server_runtime_policy(get(ENV, "JRT_DEPLOY_PROFILE", DEFAULT_DEPLOY_PROFILE))
    env_map = runtime_policy_env(policy.profile)
    for (k, v) in env_map
        haskey(ENV, k) || (ENV[k] = v)
    end
    app = FullServerApp.build_app(repo_root)
    warmup_profile = resolve_server_warmup_profile(get(ENV, "JRT_SERVER_WARMUP_PROFILE", policy.warmup_profile))
    warmup_strict = get(ENV, "JRT_SERVER_WARMUP_STRICT", policy.warmup_strict ? "1" : "0")

    print_banner(port, policy)

    if warmup_profile !== :none
        @info "Running server warmup" warmup_profile=warmup_profile strict=warmup_strict
        run_server_warmup_from_env()
    end

    try
        HTTP.serve(app, "0.0.0.0", port; verbose=false)
    catch e
        if e isa Base.IOError || e isa ArgumentError
            @error "无法启动服务器: 端口 $port 可能已被占用" exception=e
            println("\n💡 尝试使用其他端口:")
            println("   julia scripts/server/server_full.jl 8081")
            rethrow(e)
        else
            rethrow(e)
        end
    end
end

end
