# Stage-C Feasibility And Replay Evidence

This group contains four versioned Stage-C diagnostic packages:

- `cep_hybrid_stagec_offline_feasibility_v1/`: offline semantic feasibility replay.
- `cep_hybrid_stagec_certificate_feasibility_v2/`: Julia-core certificate feasibility replay.
- `cep_hybrid_stagec_extrema_guard_feasibility_v1/`: discrete-extrema guard feasibility candidate.
- `cep_hybrid_stagec_tolerance_replay_v1/`: tolerance-contract replay combining the base and revalidation runs.

The packages answer related but distinct questions and retain separate source runs, caps, verdicts, and evidence tables. They are grouped physically without merging their provenance or promotion boundaries.

These artifacts are solver-free diagnostic evidence. A `feasible_candidate` result authorizes at most the next focused production investigation; it does not promote a method to production, reference, or transport use.
