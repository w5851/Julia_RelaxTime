---
name: literature-reproduction-spike
description: Run an isolated Julia_RelaxTime literature-reproduction spike for a specified paper figure, formula, table, threshold, or numerical claim. Use a tempN sandbox, audit stated and unstated conventions, and maintain an operational status; issue the scientific verdict aligned or insufficient-information only after the verdict gate is complete.
---

# Literature Reproduction Spike

## Outcome model

Track an operational status throughout the work:

- `in_progress`: implementation or evidence collection is incomplete.
- `blocked`: a concrete dependency, access issue, or technical blocker prevents completion.
- `verdict_ready`: the independent implementation and ambiguity audit satisfy the verdict gate.

Only `verdict_ready` work may end with one of these scientific verdicts:

- `与文献对齐`
- `文献信息不足以复现`

Do not convert implementation failure, unavailable source material, or an unfinished audit into either scientific verdict.

## Isolation rules

- Work under a dedicated `tempN/<paper-or-topic>_reproduction_spike/` sandbox.
- Keep `src/`, `tests/`, `docs/`, and `config/` read-only unless the user separately authorizes mainline implementation.
- Do not call an existing high-level project implementation to generate the target result.
- Reimplement the target formula, solver path, and post-processing from the cited literature.
- Separate `explicit in paper`, `directly derived`, and `unstated choice` evidence.
- Never tune hidden parameters merely to make a curve look similar.

## Sandbox contract

Include:

- `README.md` with target, scope, non-goals, and acceptance criteria;
- `literature_facts.md` with page/equation/figure provenance;
- minimal reproduction and plotting scripts;
- `output/` containing CSV, figures, diagnostics, and configuration summaries;
- a status/verdict note with evidence and remaining blockers.

Use `paper-figure-digitize` when curve extraction, coordinate calibration, or paper-background overlays are required, and keep those artifacts inside the sandbox.

## Workflow

1. Lock the paper, target figure/table/equation, numerical feature, and error criterion.
2. Build the literature fact table before coding.
3. Implement the smallest independent loop and retain paper-symbol-to-code-variable mapping.
4. Run the most literal stated convention first and save intermediate diagnostics.
5. If results disagree, parameterize only unstated choices capable of changing topology, thresholds, scale, or branch selection.
6. Exclude false causes such as units, signs, grids, integration precision, initial values, phase unwrapping, or root selection.
7. Compare the target, first pass, and justified variants with quantitative error measures where possible.
8. Update operational status; evaluate the verdict gate only when the implementation and audit are complete.

## Verdict gate

Use `与文献对齐` only when:

- the key qualitative features align;
- important locations, thresholds, peaks, jumps, or shapes agree within a stated error;
- the necessary convention is determined by the paper or required upstream references;
- the result does not depend on hidden tuning.

Use `文献信息不足以复现` only when:

- the independent implementation is complete;
- key unstated conventions were tested;
- multiple reasonable conventions produce materially different results;
- the paper and required upstream references cannot identify the actual convention.

Otherwise retain `in_progress` or `blocked` and report the missing work or dependency.

## Final report

Report the operational status, scientific verdict if available, sandbox path, key files, tested variants, stated facts, result-sensitive unstated choices, quantitative comparison, and residual blockers. If the verdict is `文献信息不足以复现`, list the specific missing conventions rather than imposing an arbitrary minimum count.

After a successful verdict, enter mainline implementation only through a separate authorized task using the relevant Julia, regression, API-documentation, or paper-writing skills.
