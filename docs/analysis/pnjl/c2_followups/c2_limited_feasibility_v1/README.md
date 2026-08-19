# C2 Limited Feasibility v1 Inputs

This directory freezes the 17 CEP failure brackets used by the Issue #130
limited-feasibility workflow. The table is an input contract, not a numerical
result and not a production/reference artifact.

- `cep_failures.csv` is copied from the C2 convergence audit and is immutable
  for this feasibility version.
- The first runnable scope is Stage-C density. CEP and crossover remain
  deferred scope plans until the density gate is reviewed.
- No equilibrium solver is called by the aggregate replay; numerical fine-pool
  jobs are the only solver-calling part of the workflow.

The calculation SHA, source run, postprocess SHA and artifact hashes for an
actual Actions run are recorded in the external evidence package, not in this
frozen input directory.
