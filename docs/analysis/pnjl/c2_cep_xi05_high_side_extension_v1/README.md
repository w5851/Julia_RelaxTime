# C2 `xi=0.5` high-side extension v1

This versioned diagnostic plan extends the author-reviewed `xi=0.5` CEP
observations above `T=107.0625 MeV` by three equal steps of `0.0625 MeV`:

`107.125`, `107.1875`, and `107.25 MeV`.

The numerical job evaluates each temperature independently with the frozen
calculation SHA. The runner does not read oracle labels or author decisions to
select a point. The source artifact remains numerical evidence; the aggregate
collector is solver-free and produces plots for the next bracket decision.

This is diagnostic-only. It does not modify C0/C1/C2, the production CEP
workflow, Maxwell semantics, reference promotion, or transport.
