# Test Entrypoints

`test/` is the unified test entry layer for this repository.

- `test/runtests.jl`: standard `Pkg.test` / `TestItemRunner` entrypoint.
- `test/unit.jl`: thin wrapper for `tests/unit/runtests.jl`.
- `test/integration.jl`: thin wrapper for `tests/integration/runtests.jl`.
- `test/regression.jl`: thin wrapper for `tests/regression/runtests.jl`.
- `test/validation.jl`: thin wrapper for `tests/validation/runtests.jl`.

The `tests/` tree remains the implementation layer that stores the actual test
families, fixtures, and data.