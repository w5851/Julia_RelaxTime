# A startup pipe keeps the wrapper from running until its process scope is known.
# The payload is always passed as argv, never interpolated into the starter code.
@static if Sys.iswindows()
    const _ScriptTaskScope = Ptr{Cvoid}
    const _SCRIPT_TASK_STARTER = String[
        joinpath(Sys.BINDIR, Base.julia_exename()), "--startup-file=no", "--threads=1", "--compile=min", "-e",
        "readline(stdin) == \"go\" || exit(125); exit(run(pipeline(ignorestatus(Cmd(ARGS)), stdout=stdout, stderr=stderr)).exitcode)",
        "--",
    ]

    struct _JobBasicAccounting
        user_time::Int64
        kernel_time::Int64
        period_user_time::Int64
        period_kernel_time::Int64
        page_faults::UInt32
        total_processes::UInt32
        active_processes::UInt32
        terminated_processes::UInt32
    end

    _close_script_task_scope(scope::_ScriptTaskScope) =
        ccall((:CloseHandle, "kernel32"), stdcall, Cint, (Ptr{Cvoid},), scope)

    function _open_script_task_scope(proc::Base.Process)
        scope = ccall((:CreateJobObjectW, "kernel32"), stdcall, Ptr{Cvoid},
            (Ptr{Cvoid}, Ptr{UInt16}), C_NULL, C_NULL)
        Base.windowserror("CreateJobObjectW", scope == C_NULL)
        try
            # PROCESS_SET_QUOTA | PROCESS_TERMINATE
            handle = ccall((:OpenProcess, "kernel32"), stdcall, Ptr{Cvoid},
                (UInt32, Cint, UInt32), 0x0101, 0, getpid(proc))
            Base.windowserror("OpenProcess", handle == C_NULL)
            try
                ok = ccall((:AssignProcessToJobObject, "kernel32"), stdcall, Cint,
                    (Ptr{Cvoid}, Ptr{Cvoid}), scope, handle)
                Base.windowserror("AssignProcessToJobObject", ok == 0)
            finally
                _close_script_task_scope(handle)
            end
        catch
            _close_script_task_scope(scope)
            rethrow()
        end
        return scope
    end

    function _terminate_script_task_scope(scope::_ScriptTaskScope)
        ok = ccall((:TerminateJobObject, "kernel32"), stdcall, Cint,
            (Ptr{Cvoid}, UInt32), scope, 1)
        Base.windowserror("TerminateJobObject", ok == 0)
        return nothing
    end

    function _script_task_scope_active(scope::_ScriptTaskScope)
        info = Ref{_JobBasicAccounting}()
        ok = ccall((:QueryInformationJobObject, "kernel32"), stdcall, Cint,
            (Ptr{Cvoid}, Cint, Ref{_JobBasicAccounting}, UInt32, Ptr{UInt32}),
            scope, 1, info, sizeof(_JobBasicAccounting), C_NULL) # JobObjectBasicAccountingInformation
        Base.windowserror("QueryInformationJobObject", ok == 0)
        return info[].active_processes > 0
    end
else
    const _ScriptTaskScope = Cint
    const _SCRIPT_TASK_STARTER = String["sh", "-c",
        "IFS= read -r gate && [ \"\$gate\" = go ] || exit 125; exec \"\$@\"", "script-task"]
    _open_script_task_scope(proc::Base.Process) = getpid(proc)  # libuv detach creates a new session
    _close_script_task_scope(::_ScriptTaskScope) = nothing

    function _signal_script_task_scope(group::_ScriptTaskScope, signal)
        rc = ccall(:kill, Cint, (Cint, Cint), -group, signal)
        rc == 0 && return true
        errno = Libc.errno()
        errno == Libc.ESRCH && return false
        throw(SystemError("kill process group", errno))
    end

    _terminate_script_task_scope(group::_ScriptTaskScope) = _signal_script_task_scope(group, Base.SIGKILL)

    function _script_task_scope_active(group::_ScriptTaskScope)
        _signal_script_task_scope(group, 0) || return false
        # kill(group, 0) also sees defunct children awaiting reaping. They cannot write output.
        rows = read(addenv(`ps -Ao pgid=,stat=`, "LC_ALL" => "C"), String)
        return any(eachsplit(rows, '\n')) do row
            fields = split(row)
            length(fields) == 2 && tryparse(Cint, fields[1]) == group && first(fields[2]) != 'Z'
        end
    end
end

struct _ScriptTaskProcess
    proc::Base.Process
    scope::_ScriptTaskScope
end

function _wait_script_task_process(tree::_ScriptTaskProcess)
    wait(tree.proc)
    while _script_task_scope_active(tree.scope)
        sleep(0.05)
    end
    return iszero(tree.proc.termsignal) ? tree.proc.exitcode : 128 + tree.proc.termsignal
end

function _start_script_task_process(parts::Vector{String}, out, err; dir=_PROJECT_ROOT)
    starter = Cmd(Cmd([_SCRIPT_TASK_STARTER; parts]); dir, detach=!Sys.iswindows())
    gate = Pipe()
    try
        proc = run(pipeline(starter; stdin=gate, stdout=out, stderr=err); wait=false)
        scope = nothing
        try
            close(gate.out)
            scope = _open_script_task_scope(proc)
            write(gate.in, "go\n")
            return _ScriptTaskProcess(proc, scope)
        catch
            if scope === nothing
                kill(proc, Base.SIGKILL)  # starter has not received permission to spawn
                wait(proc)
            else
                try
                    _terminate_script_task_scope(scope)
                    _wait_script_task_process(_ScriptTaskProcess(proc, scope))
                finally
                    _close_script_task_scope(scope)
                end
            end
            rethrow()
        end
    finally
        close(gate)
    end
end
