# Baseline Test Results Summary

This document summarizes the baseline test results from the existing test suite before the parameter struct migration.

## Test Execution Date
Generated during task 1 setup of parameter struct migration.

## Baseline Test Run

Command executed:
```bash
julia --project=. -e 'include("tests/unit/runtests.jl")'
```

## Results Summary

### Overall Results
- **Profile**: Smoke (default)
- **Total Tests**: 310
- **Passed**: 309
- **Failed**: 0
- **Errored**: 1
- **Broken**: 0
- **Total Time**: ~1m29.4s

### Test Breakdown

#### 1. PNJL solve(FixedMu) random physical smoke
- **Tests**: 154
- **Status**: ✅ All passed
- **Time**: 6.7s

#### 2. GaussLegendre Module
- **Tests**: 38
- **Status**: ✅ All passed
- **Time**: 2.6s

#### 3. CauchyPV Module
- **Tests**: 48
- **Status**: ✅ All passed
- **Time**: 1.8s

#### 4. B0_correction Module
- **Tests**: 70 (69 passed, 1 errored)
- **Status**: ⚠️ One error (pre-existing)
- **Time**: 21.5s
- **Error Details**: 
  - Test: "精度测试" (Precision test)
  - Issue: `MethodError: no method matching B0_correction(...; rtol::Float64)`
  - The method doesn't support the `rtol` keyword argument
  - **Note**: This is a pre-existing issue, not related to struct migration

### Test Categories

1. **基本功能测试** (Basic functionality tests): 2 passed
2. **ξ = 0 时应该返回零** (Should return zero when ξ = 0): 2 passed
3. **不同 ξ 值的系统性测试** (Systematic tests for different ξ values): 11 passed
4. **k = 0 情况下的 ξ 依赖性** (ξ dependency for k = 0): 4 passed
5. **不同物理参数下的 ξ 依赖性** (ξ dependency for different physical parameters): 3 passed
6. **精度测试** (Precision tests): 1 errored
7. **ξ 线性近似验证** (ξ linear approximation verification): 1 passed
8. **B0_correction 相对于 B0 的大小比较** (Size comparison of B0_correction relative to B0): 10 passed
9. **不同 k 值下 B0_correction 相对于 B0 的比较** (Comparison of B0_correction relative to B0 for different k values): 18 passed
10. **B0_correction/B0 比值的 ξ 依赖性** (ξ dependency of B0_correction/B0 ratio): 18 passed

## Backward Compatibility Target

The parameter struct migration must maintain these results:
- All 309 passing tests must continue to pass
- The 1 pre-existing error should remain unchanged (not fixed or broken further)
- Numerical results should be identical within floating-point tolerance (rtol=1e-12, atol=1e-14)

## Verification Strategy

After completing the migration:
1. Run the same test command: `julia --project=. -e 'include("tests/unit/runtests.jl")'`
2. Verify that all 309 tests still pass
3. Verify that the same 1 error still occurs (B0_correction precision test)
4. Compare numerical outputs to ensure no changes

## Notes

- The baseline was captured with the smoke test profile (default)
- Full test suite can be run with `UNIT_PROFILE=full` environment variable
- Performance tests are excluded by default (can be included with `UNIT_INCLUDE_PERF=1`)
- Some tests are in the default skip list (see `tests/unit/runtests.jl` for details)

## Files

- Full baseline output: `baseline_results.txt`
- This summary: `BASELINE_SUMMARY.md`
