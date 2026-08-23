#!/usr/bin/env julia

"""
    magnetic_residual_context_probe.jl

Historical A/B-M1 entrypoint retained as a no-op marker. The old experiment
compared two finite-difference residual implementations; the magnetic
stationarity residual is now AD-only, so that comparison is no longer a valid
or maintained performance path. Its recorded result remains in the active
performance ledger.
"""

println("status=retired diagnostic=A/B-M1 finite-difference residual comparison")
println("reason=magnetic stationarity residual production path is ForwardDiff/AD-only")
println("solver_called=false production_written=false")
