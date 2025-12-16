# Debug scripts

This folder holds one-off diagnostics, cross-language comparisons, and exploratory scripts.

Convention:
- `scripts/relaxtime/`: “production” scripts that generate data/figures/results.
- `scripts/debug/relaxtime/`: debug/compare/find scripts that help validate implementation details.

These scripts typically:
- locate the Julia project root by searching for `Project.toml` upward from `@__DIR__`, so they can be moved without breaking;
- may depend on external artifacts (e.g. C++ debug output files).
