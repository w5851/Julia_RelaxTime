# Literature To Implementation Protocol

Status: project routing protocol.

`Julia_RelaxTime` uses literature review to support computational decisions: model formulas, numerical methods, parameter choices, validation targets, regression baselines, documentation, and paper handoff notes.

The main citation and BibTeX library remains in `D:\Desktop\paper\bib`.

## Project Agents

- `relax-literature-search-strategist`: targeted searches for formulas, parameterizations, algorithms, validation points, and reproducibility signals.
- `relax-method-reviewer`: method assumptions, equation correctness, units, reproducibility, implementation risk, and validation needs.
- `relax-evidence-synthesizer`: method comparison, engineering decision paths, validation plans, and documentation implications.
- `relax-gap-analyst`: gaps between literature and implementation, turned into experiments, regression tests, docs tasks, or paper handoffs.

## Boundary With paper

- Do not edit `D:\Desktop\paper\bib` from this project.
- Do not create a competing master bibliography in this repository.
- When citation metadata, duplicate cleanup, or citekey normalization is needed, hand off to the `paper-citation-curator` in `D:\Desktop\paper`.
- Keep local literature notes focused on implementation relevance, not complete bibliography governance.

## Evidence-To-Code Rules

- Literature evidence alone is not enough to change source code, tests, or baselines.
- Every implementation recommendation should state assumptions, equations, units, parameter definitions, validation targets, and regression risk.
- Use repository test layering: unit, integration, regression, validation, and benchmark as appropriate.
- Do not loosen tolerances to match literature without explaining numerical and physical meaning.
- Treat archived development docs as historical evidence unless explicitly requested for provenance.

## Expected Handoff

For a literature-backed implementation question:

1. Define the code/model question.
2. Search for targeted evidence.
3. Review methods for reproducibility and assumptions.
4. Synthesize options into a recommended implementation path.
5. Define validation and regression coverage.
6. Hand citation cleanup or manuscript positioning back to `D:\Desktop\paper` when needed.
