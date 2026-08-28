# Issue #130 RS transport energy convergence gate

Verdict: `production-grade`

- source commit: `05be2c05186f8e12baf3097b68f8619e53d19711`
- GitHub Actions runs: `4`
- p104 -> p128 worst relative drift: `0.0092405457`
- tau old/new invariance worst relative drift: `0.00092242053`
- xi=0 transport old/new invariance worst relative drift: `8.3725102e-08`
- failed checks: `0`

This gate covers the old convergence audit's most sensitive mode-A and mode-B panels
at xi = -0.5, -0.1, 0.0, 0.35, 0.5.  It is a prerequisite for the full
high-resolution production rerun; the Action artifacts themselves remain diagnostic-only.
