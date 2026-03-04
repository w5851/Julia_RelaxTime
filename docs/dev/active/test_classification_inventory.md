# 测试文件分类清单

生成日期：2026-03-04
依据：`docs/dev/active/2026-03-03_项目架构待优化点评估与治理计划.md` §10

---

## 分类定义

| 标签 | 说明 | 目标位置 |
|------|------|---------|
| **UNIT** | 纯单元测试，FIRST 原则 | `tests/unit/<subsystem>/` |
| **INTEGRATION** | 多模块协作/工作流测试 | `tests/integration/<subsystem>/` |
| **VALIDATION** | 基线回归验证 | `tests/validation/<subsystem>/` |
| **BENCHMARK** | BenchmarkTools 性能基准 | `benchmark/<subsystem>/` |
| **ANALYSIS** | 探索/诊断脚本 | `scripts/analysis/` |
| **PERF** | profiling/timing 探针 | `scripts/perf/` |

---

## 1. tests/unit/relaxtime/ (31 .jl)

| # | 文件 | 分类 | 迁移动作 |
|---|------|------|---------|
| 1 | `test_aniso_A_switch_smoke.jl` | UNIT | 保留（重命名去掉 _smoke 可选） |
| 2 | `test_average_scattering_rate.jl` | UNIT | 保留 |
| 3 | `test_b0_correction.jl` | UNIT | 保留 |
| 4 | `test_bulk_viscosity_derivatives.jl` | UNIT | 保留 |
| 5 | `test_default_lambda_cutoff.jl` | UNIT | 保留 |
| 6 | `test_differential_cross_section.jl` | UNIT | 保留 |
| 7 | `test_effective_couplings.jl` | UNIT | 保留 |
| 8 | `test_frame_transformations.jl` | UNIT | 保留 |
| 9 | `test_hybrid_strategy.jl` | UNIT | 保留 |
| 10 | `test_meson_mass_mott_transition.jl` | UNIT | 保留 |
| 11 | `test_meson_mass_workflow_smoke.jl` | INTEGRATION | → `tests/integration/relaxtime/` |
| 12 | `test_meson_propagator.jl` | UNIT | 保留 |
| 13 | `test_new_scattering_processes.jl` | UNIT | 保留 |
| 14 | `test_oneloopintegrals.jl` | UNIT | 保留 |
| 15 | `test_oneloopintegrals_aniso.jl` | UNIT | 保留 |
| 16 | `test_particle_symbols.jl` | UNIT | 保留 |
| 17 | `test_polarization_cache.jl` | UNIT | 保留 |
| 18 | `test_relaxation_time.jl` | UNIT | 保留 |
| 19 | `test_total_cross_section_model_ready_fixture_smoke.jl` | VALIDATION | → `tests/validation/relaxtime/` |
| 20 | `test_transport_coefficients.jl` | UNIT | 保留 |
| 21 | `test_transport_legacy_models_bridge_smoke.jl` | INTEGRATION | → `tests/integration/relaxtime/` |
| 22 | `test_transport_workflow.jl` | INTEGRATION | → `tests/integration/relaxtime/` |
| 23 | `test_transport_workflow_smoke.jl` | INTEGRATION | → `tests/integration/relaxtime/` |
| 24 | `test_transport_workflow_solver_backend_switch_smoke.jl` | INTEGRATION | → `tests/integration/relaxtime/` |
| 25 | `test_transport_workflow_toml_prefer_energy_aniso_smoke.jl` | INTEGRATION | → `tests/integration/relaxtime/` |
| 26 | `test_workflow_paramtypes_contract_fields_smoke.jl` | INTEGRATION | → `tests/integration/relaxtime/` |
| 27 | `test_workflow_paramtypes_depwarn_smoke.jl` | INTEGRATION | → `tests/integration/relaxtime/` |
| 28 | `test_workflow_paramtypes_equivalence_smoke.jl` | INTEGRATION | → `tests/integration/relaxtime/` |
| 29 | `test_workflow_paramtypes_mixedmode_smoke.jl` | INTEGRATION | → `tests/integration/relaxtime/` |
| 30 | `test_workflow_paramtypes_validation_smoke.jl` | INTEGRATION | → `tests/integration/relaxtime/` |

## 2. tests/unit/models/ (29 .jl)

| # | 文件 | 分类 | 迁移动作 |
|---|------|------|---------|
| 1 | `test_gap_residual_generic_smoke.jl` | INTEGRATION | → `tests/integration/models/` |
| 2 | `test_models_derivatives_dual_smoke.jl` | INTEGRATION | → `tests/integration/models/` |
| 3 | `test_models_dispatch_interface_smoke.jl` | INTEGRATION | → `tests/integration/models/` |
| 4 | `test_models_implicitdiff_njl2_smoke.jl` | INTEGRATION | → `tests/integration/models/` |
| 5 | `test_models_implicitdiff_smoke.jl` | INTEGRATION | → `tests/integration/models/` |
| 6 | `test_models_legacy_adapter_smoke.jl` | INTEGRATION | → `tests/integration/models/` |
| 7 | `test_models_legacy_njl_smoke.jl` | INTEGRATION | → `tests/integration/models/` |
| 8 | `test_models_native_solver_phase1_smoke.jl` | INTEGRATION | → `tests/integration/models/` |
| 9 | `test_models_phase0_smoke.jl` | INTEGRATION | → `tests/integration/models/` |
| 10 | `test_models_unified_entrypoints_smoke.jl` | INTEGRATION | → `tests/integration/models/` |
| 11 | `test_njl2_model_factory.jl` | UNIT | 保留 |
| 12 | `test_njl2_omega.jl` | UNIT | 保留 |
| 13 | `test_njl_model_factory.jl` | UNIT | 保留 |
| 14 | `test_njl_omega.jl` | UNIT | 保留 |
| 15 | `test_phase_artifacts_promotion_smoke.jl` | INTEGRATION | → `tests/integration/models/` |
| 16 | `test_phase_cli_direct_smoke.jl` | INTEGRATION | → `tests/integration/models/` |
| 17 | `test_phase_core_algorithms_smoke.jl` | INTEGRATION | → `tests/integration/models/` |
| 18 | `test_phase_legacy_path_detach_smoke.jl` | INTEGRATION | → `tests/integration/models/` |
| 19 | `test_phase_pipeline_smoke.jl` | INTEGRATION | → `tests/integration/models/` |
| 20 | `test_pnjl_gk_polyakov_injection_smoke.jl` | INTEGRATION | → `tests/integration/models/` |
| 21 | `test_pnjl_integrals_forwarddiff_smoke.jl` | INTEGRATION | → `tests/integration/models/` |
| 22 | `test_pnjl_lambda_injection_smoke.jl` | INTEGRATION | → `tests/integration/models/` |
| 23 | `test_pnjl_models_integrals_smoke.jl` | INTEGRATION | → `tests/integration/models/` |
| 24 | `test_pnjl_params_injection_smoke.jl` | INTEGRATION | → `tests/integration/models/` |
| 25 | `test_pnjl_solve_gap_backend_switch_smoke.jl` | INTEGRATION | → `tests/integration/models/` |
| 26 | `test_pnjl_solve_gap_generic_smoke.jl` | INTEGRATION | → `tests/integration/models/` |
| 27 | `test_pnjl_thermo_bridge_multipoint_smoke.jl` | INTEGRATION | → `tests/integration/models/` |
| 28 | `test_rpnjl_bridge_smoke.jl` | INTEGRATION | → `tests/integration/models/` |
| 29 | `test_rpnjl_model_factory_smoke.jl` | INTEGRATION | → `tests/integration/models/` |

## 3. tests/unit/pnjl/ (32 .jl)

| # | 文件 | 分类 | 迁移动作 |
|---|------|------|---------|
| 1 | `test_aniso_gap_solver.jl` | UNIT | 保留 |
| 2 | `test_bulk_viscosity.jl` | UNIT | 保留 |
| 3 | `test_constraint_fixedpoint_baseline_smoke.jl` | VALIDATION | → `tests/validation/pnjl/` |
| 4 | `test_constraint_modes_parity_smoke.jl` | INTEGRATION | → `tests/integration/pnjl/` |
| 5 | `test_core_integrals.jl` | UNIT | 保留 |
| 6 | `test_core_thermodynamics.jl` | UNIT | 保留 |
| 7 | `test_implicit_jacobian.jl` | UNIT | 保留 |
| 8 | `test_magnetic_coupling_GB.jl` | UNIT | 保留 |
| 9 | `test_magnetic_density.jl` | UNIT | 保留 |
| 10 | `test_magnetic_energy_levels.jl` | UNIT | 保留 |
| 11 | `test_magnetic_fixedpoint_baseline_smoke.jl` | VALIDATION | → `tests/validation/pnjl/` |
| 12 | `test_magnetic_nmax_convergence.jl` | ANALYSIS | → `scripts/analysis/pnjl/` |
| 13 | `test_magnetic_omega_components.jl` | UNIT | 保留 |
| 14 | `test_pnjl_quark_distribution_sum.jl` | UNIT | 保留 |
| 15 | `test_quark_distribution_aniso.jl` | UNIT | 保留 |
| 16 | `test_quark_distribution_antiderivative.jl` | UNIT | 保留 |
| 17 | `test_scans.jl` | INTEGRATION | → `tests/integration/pnjl/` |
| 18 | `test_scan_config_equivalence_smoke.jl` | INTEGRATION | → `tests/integration/pnjl/` |
| 19 | `test_scan_fixedpoint_baseline_smoke.jl` | VALIDATION | → `tests/validation/pnjl/` |
| 20 | `test_scan_input_contract_validation_smoke.jl` | INTEGRATION | → `tests/integration/pnjl/` |
| 21 | `test_scan_solver_boundary_error_smoke.jl` | INTEGRATION | → `tests/integration/pnjl/` |
| 22 | `test_solver_conditions.jl` | UNIT | 保留 |
| 23 | `test_solver_constraints_models_backend_smoke.jl` | INTEGRATION | → `tests/integration/pnjl/` |
| 24 | `test_solver_constraint_modes.jl` | UNIT | 保留 |
| 25 | `test_solver_implicit.jl` | UNIT | 保留 |
| 26 | `test_solver_random_physical_smoke.jl` | INTEGRATION | → `tests/integration/pnjl/` |
| 27 | `test_solver_seed_strategies.jl` | UNIT | 保留 |
| 28 | `test_thermo_derivatives.jl` | UNIT | 保留 |
| 29 | `test_tmu_scan_smoke.jl` | INTEGRATION | → `tests/integration/pnjl/` |
| 30 | `test_tmu_scan_solver_backend_models_smoke.jl` | INTEGRATION | → `tests/integration/pnjl/` |
| 31 | `test_trho_scan_smoke.jl` | INTEGRATION | → `tests/integration/pnjl/` |
| 32 | `test_trho_scan_solver_backend_models_smoke.jl` | INTEGRATION | → `tests/integration/pnjl/` |

## 4. tests/unit/njl/ (2 .jl)

| # | 文件 | 分类 | 迁移动作 |
|---|------|------|---------|
| 1 | `test_njl_core.jl` | UNIT | 保留 |
| 2 | `test_njl2_core.jl` | UNIT | 保留 |

## 5. tests/unit/config/ (3 .jl)

| # | 文件 | 分类 | 迁移动作 |
|---|------|------|---------|
| 1 | `test_config_profile_smoke.jl` | INTEGRATION | → `tests/integration/config/` |
| 2 | `test_pnjl_profile_dynamic_constants_smoke.jl` | INTEGRATION | → `tests/integration/config/` |
| 3 | `test_pnjl_rpnjl_critical_params_smoke.jl` | INTEGRATION | → `tests/integration/config/` |

## 6. tests/unit/integration/ (3 .jl) → 改名为 tests/unit/numerics/

| # | 文件 | 分类 | 迁移动作 |
|---|------|------|---------|
| 1 | `test_cauchypv.jl` | UNIT | 保留（目录改名） |
| 2 | `test_gausslegendre.jl` | UNIT | 保留（目录改名） |
| 3 | `test_momentum_mapping.jl` | UNIT | 保留（目录改名） |

## 7. tests/unit/struct_migration/ (22 .jl) → 拆散

| # | 文件 | 分类 | 迁移动作 |
|---|------|------|---------|
| 1 | `runtests.jl` | — | 删除（入口不再需要） |
| 2 | `baseline_results.txt` | — | 归档或删除 |
| 3 | `test_parameter_types.jl` | UNIT | → `tests/unit/types/` |
| 4 | `test_particle_symbols_struct.jl` | UNIT | → `tests/unit/relaxtime/` 合并 |
| 5 | `test_particle_symbols_edge_cases.jl` | UNIT | → `tests/unit/relaxtime/` 合并 |
| 6 | `test_ensure_quark_params_has_A_property.jl` | UNIT | → `tests/unit/relaxtime/` 合并 |
| 7 | `test_scattering_amplitude_minimal.jl` | UNIT | → `tests/unit/relaxtime/` 合并 |
| 8 | `test_scattering_amplitude_property.jl` | UNIT | → `tests/unit/relaxtime/` 合并 |
| 9 | `test_scattering_amplitude_edge_cases.jl` | UNIT | → `tests/unit/relaxtime/` 合并 |
| 10 | `test_differential_cross_section_property.jl` | UNIT | → `tests/unit/relaxtime/` 合并 |
| 11 | `test_differential_cross_section_edge_cases.jl` | UNIT | → `tests/unit/relaxtime/` 合并 |
| 12 | `test_total_propagator_property.jl` | UNIT | → `tests/unit/relaxtime/` 合并 |
| 13 | `test_total_propagator_edge_cases.jl` | UNIT | → `tests/unit/relaxtime/` 合并 |
| 14 | `test_total_cross_section_property.jl` | UNIT | → `tests/unit/relaxtime/` 合并 |
| 15 | `test_total_cross_section_edge_cases.jl` | UNIT | → `tests/unit/relaxtime/` 合并 |
| 16 | `test_average_scattering_rate_property.jl` | UNIT | → `tests/unit/relaxtime/` 合并 |
| 17 | `test_average_scattering_rate_edge_cases.jl` | UNIT | → `tests/unit/relaxtime/` 合并 |
| 18 | `test_relaxation_time_property.jl` | UNIT | → `tests/unit/relaxtime/` 合并 |
| 19 | `test_relaxation_time_edge_cases.jl` | UNIT | → `tests/unit/relaxtime/` 合并 |
| 20 | `test_utils.jl` | UNIT | → `tests/unit/relaxtime/` 合并 |
| 21 | `test_integration_full_chain.jl` | INTEGRATION | → `tests/integration/relaxtime/` |
| 22 | `test_mixed_usage_patterns.jl` | INTEGRATION | → `tests/integration/relaxtime/` |
| 23 | `test_conversion_roundtrip_property.jl` | INTEGRATION | → `tests/integration/relaxtime/` |

## 8. tests/perf/pnjl/ (3 .jl)

| # | 文件 | 分类 | 迁移动作 |
|---|------|------|---------|
| 1 | `scan_perf.jl` | PERF | → `scripts/perf/pnjl/` |
| 2 | `single_point_solver_perf.jl` | BENCHMARK | → `benchmark/pnjl/` |
| 3 | `test_bulk_derivative_coeffs_performance.jl` | BENCHMARK | → `benchmark/pnjl/` |

## 9. tests/perf/relaxtime/ (16 .jl + 3 .md)

| # | 文件 | 分类 | 迁移动作 |
|---|------|------|---------|
| 1 | `benchmark_adaptive_nodes.jl` | BENCHMARK | → `benchmark/relaxtime/` |
| 2 | `benchmark_adaptive_nodes_v2.jl` | BENCHMARK | → `benchmark/relaxtime/` |
| 3 | `benchmark_average_scattering_rate_N60_w0cdf_pchip.jl` | BENCHMARK | → `benchmark/relaxtime/` |
| 4 | `benchmark_bulk_viscosity_derivatives.jl` | BENCHMARK | → `benchmark/relaxtime/` |
| 5 | `benchmark_struct_vs_namedtuple.jl` | BENCHMARK | → `benchmark/relaxtime/` |
| 6 | `benchmark_total_cross_section_hotspot.jl` | BENCHMARK | → `benchmark/relaxtime/` |
| 7 | `benchmark_with_roots.jl` | BENCHMARK | → `benchmark/relaxtime/` |
| 8 | `comprehensive_benchmark.jl` | BENCHMARK | → `benchmark/relaxtime/` |
| 9 | `compare_quadgk_vs_hybrid.jl` | ANALYSIS | → `scripts/analysis/relaxtime/` |
| 10 | `performance_comparison.jl` | PERF | → `scripts/perf/relaxtime/` |
| 11 | `profile_hybrid_overhead.jl` | PERF | → `scripts/perf/relaxtime/` |
| 12 | `profile_hybrid_v2.jl` | PERF | → `scripts/perf/relaxtime/` |
| 13 | `test_average_scattering_rate_performance.jl` | BENCHMARK | → `benchmark/relaxtime/` |
| 14 | `test_total_cross_section_performance.jl` | BENCHMARK | → `benchmark/relaxtime/` |
| 15 | `test_total_propagator_performance.jl` | BENCHMARK | → `benchmark/relaxtime/` |
| 16 | `timing_analysis.jl` | PERF | → `scripts/perf/relaxtime/` |
| — | `README.md` | — | → `benchmark/relaxtime/README.md` |
| — | `test_differential_cross_section_performance.md` | — | → `benchmark/relaxtime/` |
| — | `test_total_cross_section_performance.md` | — | → `benchmark/relaxtime/` |

## 10. tests/analysis/ (50 .jl + md) → 整体迁至 scripts/analysis/

全部分类为 **ANALYSIS**，迁移到 `scripts/analysis/`，保留原子目录结构。

**根目录** (18 .jl + 3 .md):
`analyze_ss_scattering_issue.jl`, `check_mandelstam_thresholds.jl`, `deep_trace_warning.jl`,
`find_source_of_minus_5_88.jl`, `gauss_jacobi_vs_legendre.jl`, `test_cache_debug.jl`,
`test_cache_first_call_analysis.jl`, `test_cache_hit_rate.jl`, `test_cache_overhead_analysis.jl`,
`test_cache_simple.jl`, `test_cache_size_impact.jl`, `test_const_nodes.jl`,
`test_cosθ_half_precision.jl`, `test_gauss_exports.jl`, `test_sequential_calls.jl`,
`test_ss_recommended_params.jl`, `test_timing_methods_comparison.jl`, `trace_ss_warning.jl`

**convergence/** (3 .jl + 1 .md): `run_convergence.jl`, `test_A_convergence.jl`, `test_B0_convergence.jl`
**pnjl/** (2 .jl): `quadrature_node_scan.jl`, `theta_node_scan.jl`
**relaxtime/** (18 .jl): 全部分析脚本
**relaxtime_diagnostics/** (4 .jl): 诊断脚本
**relaxtime_validation/** (5 .jl): 验证脚本（宽松断言，非正式 CI 回归）
**struct_migration/** (4 .md): 文档，归档
**unit_summaries/** (12 .md): 测试摘要文档，归档

---

## 汇总统计

| 分类 | 文件数 | 说明 |
|------|--------|------|
| UNIT（保留在 tests/unit/） | 46 | 纯模块级测试 |
| INTEGRATION（→ tests/integration/） | 56 | 工作流/跨模块 smoke |
| VALIDATION（→ tests/validation/） | 4 | 基线回归 |
| BENCHMARK（→ benchmark/） | 14 | BenchmarkTools 基准 |
| ANALYSIS（→ scripts/analysis/） | 51 | 探索/诊断脚本 |
| PERF（→ scripts/perf/） | 8 | profiling/timing |
| 其他（runner/md/data） | 15 | runtests.jl/文档/数据 |
| **总计** | **194** | — |
