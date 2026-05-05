module JRTLauncherApp

export julia_main, real_main, parse_launcher_args

using Pkg

const APP_ROOT = normpath(joinpath(@__DIR__, ".."))
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const MODELS_PATH = joinpath(REPO_ROOT, "src", "models", "Models.jl")
const PHASE_SUPPORT_PATH = joinpath(REPO_ROOT, "scripts", "pnjl", "phase_cli_support.jl")

struct LauncherCommand
    name::Symbol
    passthrough::Vector{String}
end

function _usage(io::IO=stdout)
    println(io, "JRT launcher app PoC")
    println(io, "")
    println(io, "用法:")
    println(io, "  jrt-launcher help")
    println(io, "  jrt-launcher phase [phase-cli-args...]")
    println(io, "")
    println(io, "说明:")
    println(io, "  - 本 PoC 当前只承载 phase 子命令。")
    println(io, "  - phase 参数保持与 scripts/pnjl/calculate_phase_structure.jl 一致。")
    println(io, "  - 当前 bundle 仍面向仓库内本地验证，不宣称已脱离源码仓库独立分发。")
end

function parse_launcher_args(args::Vector{String})
    isempty(args) && return LauncherCommand(:help, String[])

    cmd = lowercase(strip(first(args)))
    if cmd in ("help", "-h", "--help")
        return LauncherCommand(:help, collect(String.(args[2:end])))
    elseif cmd == "phase"
        return LauncherCommand(:phase, collect(String.(args[2:end])))
    end

    throw(ArgumentError("unsupported launcher command: $(first(args))"))
end

function _ensure_phase_runtime!()
    Pkg.activate(REPO_ROOT)
    if !isdefined(Main, :Models)
        Base.include(Main, MODELS_PATH)
    end
    if !isdefined(Main, :PhaseCliSupport)
        Base.include(Main, PHASE_SUPPORT_PATH)
    end
    return nothing
end

function _invoke_phase_main(args::Vector{String})
    phase_cli = getfield(Main, :PhaseCliSupport)
    models = getfield(Main, :Models)
    phase_main = getfield(phase_cli, :main)
    return phase_main(models, REPO_ROOT, args)
end

function _run_phase(args::Vector{String})
    _ensure_phase_runtime!()
    return Base.invokelatest(_invoke_phase_main, args)
end

function real_main(args::Vector{String}=collect(String.(ARGS)))
    cmd = parse_launcher_args(args)
    if cmd.name === :help
        _usage()
        return 0
    elseif cmd.name === :phase
        return _run_phase(cmd.passthrough)
    end
    throw(ArgumentError("unsupported launcher command: $(cmd.name)"))
end

function julia_main()::Cint
    try
        return Cint(real_main())
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
end

end # module JRTLauncherApp
