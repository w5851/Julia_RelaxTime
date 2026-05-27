using UUIDs: uuid4

const _SCRIPT_TASK_JOBS_LOCK = ReentrantLock()
const _SCRIPT_TASK_JOBS = Dict{String, Dict{String, Any}}()
const _SCRIPT_TASK_JOB_SEQ = Ref(0)
const _SCRIPT_TASK_MAX_RUNNING = 1
const _SCRIPT_TASK_MAX_PENDING = 16
const _SCRIPT_TASK_SEMAPHORE = Base.Semaphore(_SCRIPT_TASK_MAX_RUNNING)
const _SCRIPT_TASK_STATUS_QUEUED = "queued"
const _SCRIPT_TASK_STATUS_RUNNING = "running"
const _SCRIPT_TASK_STATUS_SUCCEEDED = "succeeded"
const _SCRIPT_TASK_STATUS_FAILED = "failed"
const _SCRIPT_TASK_STATUS_CANCELLED = "cancelled"
const _SCRIPT_TASK_TERMINAL = Set([
    _SCRIPT_TASK_STATUS_SUCCEEDED,
    _SCRIPT_TASK_STATUS_FAILED,
    _SCRIPT_TASK_STATUS_CANCELLED,
])
const _SCRIPT_TASK_JOB_ROOT = joinpath(_PROJECT_ROOT, "data", "outputs", "frontend_jobs")
const _SCRIPT_TASK_LOG_ROOT = joinpath(_PROJECT_ROOT, "data", "outputs", "logs", "frontend_jobs")

function _job_token_path(parts::AbstractString...)
    return joinpath("data", "outputs", "frontend_jobs", parts...)
end

function _script_task_catalog_items()
    return [
        Dict(
            "id" => "run-unified-scan",
            "script" => "scripts/models/run_unified_scan.jl",
            "name" => "统一 Models 扫描",
            "purpose" => "统一触发 Models 主链的 T-mu / T-rho 扫描，适合做前端可复现的小网格或点阵扫描。",
            "use_cases" => ["T-mu/T-rho 网格 smoke", "Models 主链扫描契约验证", "生成可追踪 CSV"],
            "key_params" => ["scan kind: tmu|trho", "model_kind", "T_values", "mu_values/rho_values", "xi_values", "output_path"],
            "default_preset" => "smoke",
            "estimated_time" => Dict("smoke" => "约 1-3 分钟（冷启动更久）", "canonical" => "数分钟到十几分钟", "custom" => "取决于网格规模"),
            "local_recommendation" => "建议本机运行 smoke；canonical 需确认网格规模和输出路径。",
            "outputs" => ["扫描 CSV", "stdout/stderr 日志"],
            "manifest_candidates" => String[],
            "presets" => Dict(
                "smoke" => Dict(
                    "label" => "smoke",
                    "heavy" => false,
                    "args" => ["scan", "tmu", "--model_kind=PNJL", "--T_values=150", "--mu_values=0", "--xi_values=0.0", "--output_path={job_output}/unified_scan_smoke.csv", "--overwrite=true"],
                ),
                "canonical" => Dict(
                    "label" => "canonical",
                    "heavy" => true,
                    "args" => ["scan", "tmu", "--model_kind=PNJL", "--T_values=150,160", "--mu_values=0,100", "--xi_values=0.0", "--output_path={job_output}/unified_scan_canonical.csv", "--overwrite=true"],
                ),
                "custom" => Dict("label" => "custom", "heavy" => true, "args" => String[]),
            ),
        ),
        Dict(
            "id" => "run-conserved-charge-susceptibilities",
            "script" => "scripts/pnjl/run_conserved_charge_susceptibilities.jl",
            "name" => "守恒荷易感性/累积量",
            "purpose" => "计算 PNJL 守恒荷广义磁化率、累积量与 baryon Ssigma/kappa_sigma2。",
            "use_cases" => ["单点易感性 smoke", "沿 T 或 muB 做轻量扫描", "导出 CSV 供论文图表或回归检查"],
            "key_params" => ["T", "muB/muQ/muS", "scan axis", "scan_min/max/step", "xi", "p_num/t_num", "output"],
            "default_preset" => "smoke",
            "estimated_time" => Dict("smoke" => "约 1-3 分钟（AD 冷启动可能更久）", "canonical" => "数分钟", "custom" => "取决于扫描轴和阶数成本"),
            "local_recommendation" => "本机适合 smoke 和小范围扫描；大范围扫描建议先确认节点数。",
            "outputs" => ["CSV 结果", "stdout/stderr 日志"],
            "manifest_candidates" => String[],
            "presets" => Dict(
                "smoke" => Dict(
                    "label" => "smoke",
                    "heavy" => false,
                    "args" => ["--T=150", "--muB=0", "--scan=none", "--xi=0.0", "--p_num=8", "--t_num=4", "--output={job_output}/susceptibilities_smoke.csv"],
                ),
                "canonical" => Dict(
                    "label" => "canonical",
                    "heavy" => true,
                    "args" => ["--scan=muB", "--scan_min=0", "--scan_max=200", "--scan_step=100", "--xi=0.0", "--p_num=24", "--t_num=8", "--output={job_output}/susceptibilities_muB_scan.csv"],
                ),
                "custom" => Dict("label" => "custom", "heavy" => true, "args" => String[]),
            ),
        ),
        Dict(
            "id" => "run-relaxtime-orchestrator",
            "script" => "scripts/relaxtime/run_relaxtime_orchestrator.jl",
            "name" => "Relaxtime 编排工作流",
            "purpose" => "按配置编排 transport 或 cross-section 工作流，生成带 run_id 的输出目录。",
            "use_cases" => ["验证 workflow config", "运行 cross-section 子集", "生成 transport 编排产物"],
            "key_params" => ["subcommand: transport|cross-section", "config", "outdir", "processes", "resume/overwrite"],
            "default_preset" => "smoke",
            "estimated_time" => Dict("smoke" => "数分钟，取决于截面过程数", "canonical" => "较长，建议独占本机", "custom" => "取决于配置"),
            "local_recommendation" => "smoke 可本机试跑；transport/canonical 应确认 CPU 占用和输出目录。",
            "outputs" => ["orchestrated 输出目录", "run_id 摘要", "stdout/stderr 日志"],
            "manifest_candidates" => ["run_manifest.json", "effective_config.json"],
            "presets" => Dict(
                "smoke" => Dict(
                    "label" => "smoke",
                    "heavy" => false,
                    "args" => ["cross-section", "--outdir", "{job_output}/orchestrator_cross_section_smoke", "--overwrite", "--processes", "uubar_to_uubar"],
                ),
                "canonical" => Dict(
                    "label" => "canonical",
                    "heavy" => true,
                    "args" => ["transport", "--outdir", "{job_output}/orchestrator_transport_canonical", "--overwrite"],
                ),
                "custom" => Dict("label" => "custom", "heavy" => true, "args" => String[]),
            ),
        ),
        Dict(
            "id" => "run-gap-transport-scan",
            "script" => "scripts/relaxtime/run_gap_transport_scan.jl",
            "name" => "Gap + Transport 扫描",
            "purpose" => "批量串联 PNJL 平衡求解、tau 与 RTA 输运系数，输出扫描 CSV 与 provenance sidecars。",
            "use_cases" => ["输运单点/小网格 smoke", "xi/T/muB 扫描", "生成 gap_transport_scan.csv"],
            "key_params" => ["tmin/tmax/tstep", "mubmin/mubmax/mubstep", "xi-list", "mode", "compute_bulk", "output", "provenance-dir"],
            "default_preset" => "smoke",
            "estimated_time" => Dict("smoke" => "数分钟", "canonical" => "几十分钟或更久", "custom" => "取决于网格和 tau 设置"),
            "local_recommendation" => "本机默认只跑 smoke；canonical 需确认节点数、bulk 和网格规模。",
            "outputs" => ["gap_transport CSV", "failed_points CSV", "run_manifest/effective_config", "stdout/stderr 日志"],
            "manifest_candidates" => ["run_manifest.json", "effective_config.json"],
            "presets" => Dict(
                "smoke" => Dict(
                    "label" => "smoke",
                    "heavy" => false,
                    "args" => ["--output", "{job_output}/gap_transport_smoke.csv", "--failed-points-output", "{job_output}/gap_transport_failed_points.csv", "--provenance-dir", "{job_output}", "--tmin", "150", "--tmax", "150", "--tstep", "10", "--mubmin", "0", "--mubmax", "0", "--mubstep", "100", "--xi-list", "0.0", "--mode", "finite_15", "--no-compute-bulk", "--p-num", "8", "--t-num", "4", "--tau-p-nodes", "8", "--tau-angle-nodes", "2", "--tau-phi-nodes", "4", "--tau-n-sigma", "3", "--sigma-grid-n", "8", "--overwrite"],
                ),
                "canonical" => Dict(
                    "label" => "canonical",
                    "heavy" => true,
                    "args" => ["--output", "{job_output}/gap_transport_canonical.csv", "--failed-points-output", "{job_output}/gap_transport_failed_points.csv", "--provenance-dir", "{job_output}", "--tmin", "120", "--tmax", "200", "--tstep", "20", "--mubmin", "0", "--mubmax", "400", "--mubstep", "200", "--xi-list", "0.0,0.2", "--mode", "finite_15", "--overwrite"],
                ),
                "custom" => Dict("label" => "custom", "heavy" => true, "args" => String[]),
            ),
        ),
        Dict(
            "id" => "run-phase-guided-transport-scan",
            "script" => "scripts/relaxtime/run_phase_guided_transport_scan.jl",
            "name" => "Phase-guided Transport 扫描",
            "purpose" => "沿相结构参考生成 phase-guided transport 采样计划并可执行扫描。",
            "use_cases" => ["dry-run 生成 sampling plan", "canonical case 扫描", "为后续 plot wrapper 准备 case_dir"],
            "key_params" => ["mode", "case-name", "outdir", "dry-run", "compute_bulk", "overwrite/resume"],
            "default_preset" => "smoke",
            "estimated_time" => Dict("smoke" => "约几十秒到数分钟（dry-run）", "canonical" => "较长，取决于 mode 和采样数", "custom" => "取决于采样计划"),
            "local_recommendation" => "smoke 仅 dry-run，适合先理解任务；实际扫描需确认。",
            "outputs" => ["sampling_plan.csv", "README.md", "effective_config.json", "run_manifest.json", "phase_guided_transport_scan.csv（非 dry-run）"],
            "manifest_candidates" => ["run_manifest.json", "effective_config.json"],
            "presets" => Dict(
                "smoke" => Dict(
                    "label" => "smoke",
                    "heavy" => false,
                    "args" => ["--mode", "fixed-T-sparse-muB", "--case-name", "frontend_smoke", "--outdir", "{job_output}/phase_guided_smoke", "--dry-run", "--overwrite"],
                ),
                "canonical" => Dict(
                    "label" => "canonical",
                    "heavy" => true,
                    "args" => ["--mode", "fixed-muB-phase-scaled", "--case-name", "frontend_canonical", "--outdir", "{job_output}/phase_guided_canonical", "--overwrite"],
                ),
                "custom" => Dict("label" => "custom", "heavy" => true, "args" => String[]),
            ),
        ),
        Dict(
            "id" => "run-phase-guided-transport-plots",
            "script" => "scripts/relaxtime/run_phase_guided_transport_plots.jl",
            "name" => "Phase-guided 图表生成",
            "purpose" => "对已有 phase-guided case_dir 生成 plot_manifest 与 PNG 图层，并回写 README。",
            "use_cases" => ["复用扫描 case_dir 做可视化", "生成 plot review 包", "检查 phase-guided 输出是否可画"],
            "key_params" => ["case-dir", "csv", "fig-dir", "python", "overwrite"],
            "default_preset" => "smoke",
            "estimated_time" => Dict("smoke" => "需要已有 case_dir；通常数十秒到数分钟", "canonical" => "取决于图数", "custom" => "取决于输入 CSV"),
            "local_recommendation" => "必须先有有效 case_dir；前端默认只展示和要求 custom 确认。",
            "outputs" => ["PNG figures", "plot_manifest.json", "README plot-review 追加段落"],
            "manifest_candidates" => ["plot_manifest.json"],
            "presets" => Dict(
                "smoke" => Dict("label" => "smoke", "heavy" => false, "args" => ["--case-dir", "{job_output}/missing_case_dir"]),
                "canonical" => Dict("label" => "canonical", "heavy" => true, "args" => ["--case-dir", "{job_output}/missing_case_dir", "--overwrite"]),
                "custom" => Dict("label" => "custom", "heavy" => true, "args" => String[]),
            ),
        ),
        Dict(
            "id" => "run-manual-relaxation-workflow",
            "script" => "scripts/relaxtime/run_manual_relaxation_scan_workflow.jl",
            "name" => "手动 Relaxation 扫描组合",
            "purpose" => "生成 cross_section、temperature_scan_muB0_xi0 和 fixed_temperature_xi_scan_muB0 三类产物。",
            "use_cases" => ["快速生成 cross_section 预览", "温度扫描/固定温度 xi 扫描", "生成 figures/relaxtime 组合产物"],
            "key_params" => ["sections", "mode", "tau nodes", "sigma-grid-n", "no-plots", "base-output-dir"],
            "default_preset" => "smoke",
            "estimated_time" => Dict("smoke" => "数分钟，默认只跑 cross_section 且不画图", "canonical" => "较长，会跑全套默认 sections", "custom" => "取决于 sections 和节点数"),
            "local_recommendation" => "建议先跑 cross_section smoke；全套 canonical 需确认。",
            "outputs" => ["CSV 结果", "figures（若开启 plots）", "run_manifest.json/effective_config.json"],
            "manifest_candidates" => ["run_manifest.json", "effective_config.json"],
            "presets" => Dict(
                "smoke" => Dict(
                    "label" => "smoke",
                    "heavy" => false,
                    "args" => ["--sections", "cross_section", "--no-plots", "--tau-p-nodes", "8", "--tau-angle-nodes", "2", "--tau-phi-nodes", "4", "--tau-n-sigma", "3", "--sigma-grid-n", "8", "--base-output-dir", "{job_output}/manual_relaxation_smoke", "--overwrite"],
                ),
                "canonical" => Dict(
                    "label" => "canonical",
                    "heavy" => true,
                    "args" => ["--sections", "all", "--base-output-dir", "{job_output}/manual_relaxation_canonical", "--overwrite"],
                ),
                "custom" => Dict("label" => "custom", "heavy" => true, "args" => String[]),
            ),
        ),
    ]
end

function _extended_script_task_catalog_items()
    return [
        Dict(
            "id" => "run-gap-meson-mass-scan",
            "script" => "scripts/relaxtime/run_gap_meson_mass_scan.jl",
            "name" => "Gap + Meson Mass/Mott 扫描",
            "purpose" => "串联 PNJL 平衡求解、介子质量/宽度与 Mott threshold gap，输出后续 Mott 定位可消费的 CSV。",
            "use_cases" => ["Mott 相关 smoke 网格", "生成 gap_meson_mass_scan.csv", "为 Mott phase 后处理准备输入"],
            "key_params" => ["tmin/tmax/tstep", "mubmin/mubmax/mubstep", "xi-list", "p-num/t-num/max-iter", "output"],
            "default_preset" => "smoke",
            "estimated_time" => Dict("smoke" => "数分钟", "canonical" => "十几分钟或更久", "custom" => "取决于 T/muB/xi 网格"),
            "local_recommendation" => "建议先跑单点 smoke；canonical 会扫描多个温度和 muB 点，应确认耗时。",
            "outputs" => ["介子质量/宽度/Mott gap CSV", "stdout/stderr 日志"],
            "manifest_candidates" => String[],
            "presets" => Dict(
                "smoke" => Dict(
                    "label" => "smoke",
                    "heavy" => false,
                    "args" => ["--output", "{job_output}/gap_meson_mass_smoke.csv", "--tmin", "150", "--tmax", "150", "--tstep", "10", "--mubmin", "0", "--mubmax", "0", "--mubstep", "100", "--xi-list", "0.0", "--p-num", "8", "--t-num", "4", "--max-iter", "20", "--overwrite"],
                ),
                "canonical" => Dict(
                    "label" => "canonical",
                    "heavy" => true,
                    "args" => ["--output", "{job_output}/gap_meson_mass_canonical.csv", "--tmin", "120", "--tmax", "220", "--tstep", "20", "--mubmin", "0", "--mubmax", "400", "--mubstep", "200", "--xi-list", "0.0,0.2", "--overwrite"],
                ),
                "custom" => Dict("label" => "custom", "heavy" => true, "args" => String[]),
            ),
        ),
        Dict(
            "id" => "run-mott-phase-scan",
            "script" => "scripts/relaxtime/run_mott_phase_scan.jl",
            "name" => "Mott Phase 扫描",
            "purpose" => "固定 muB 口径下扫描介子 Mott phase 资产，生成 scan CSV、effective_config 与 run_manifest。",
            "use_cases" => ["Mott phase smoke", "canonical xi/T 温度扫描", "为 derived CSV 与 plot modes 准备输入"],
            "key_params" => ["outdir", "output", "tmin/tmax/tstep", "xi-list", "include-mixed", "p-num/t-num/max-iter"],
            "default_preset" => "smoke",
            "estimated_time" => Dict("smoke" => "数分钟", "canonical" => "较长，取决于 T 与 xi 点数", "custom" => "取决于温度/xi 网格"),
            "local_recommendation" => "smoke 只跑极小温度窗口；canonical 建议空闲时运行。",
            "outputs" => ["mott_phase_scan.csv", "effective_config.json", "run_manifest.json", "stdout/stderr 日志"],
            "manifest_candidates" => ["run_manifest.json", "effective_config.json"],
            "presets" => Dict(
                "smoke" => Dict(
                    "label" => "smoke",
                    "heavy" => false,
                    "args" => ["--outdir", "{job_output}/mott_phase_smoke", "--output", "mott_phase_scan.csv", "--tmin", "140", "--tmax", "144", "--tstep", "4", "--xi-list", "0.0", "--p-num", "8", "--t-num", "4", "--max-iter", "20", "--overwrite"],
                ),
                "canonical" => Dict(
                    "label" => "canonical",
                    "heavy" => true,
                    "args" => ["--outdir", "{job_output}/mott_phase_canonical", "--output", "mott_phase_scan.csv", "--tmin", "120", "--tmax", "260", "--tstep", "5", "--xi-list", "-0.3,0.0,0.3", "--overwrite"],
                ),
                "custom" => Dict("label" => "custom", "heavy" => true, "args" => String[]),
            ),
        ),
        Dict(
            "id" => "run-relaxtime-scan-dispatcher",
            "script" => "scripts/relaxtime/run_scan.jl",
            "name" => "Relaxtime 子命令分发器",
            "purpose" => "统一转发 gap-transport、phase-guided-transport、tau-vs-t 与 manual-workflow 子命令。",
            "use_cases" => ["用统一入口调用 gap-transport smoke", "保留脚本迁移期兼容入口", "对比直接脚本与分发入口行为"],
            "key_params" => ["subcommand", "pass-through args", "output/outdir"],
            "default_preset" => "smoke",
            "estimated_time" => Dict("smoke" => "数分钟", "canonical" => "取决于子命令", "custom" => "完全由子命令参数决定"),
            "local_recommendation" => "更推荐前端直接选具体任务；该入口适合验证分发层。",
            "outputs" => ["转发目标脚本的输出", "stdout/stderr 日志"],
            "manifest_candidates" => ["run_manifest.json", "effective_config.json"],
            "presets" => Dict(
                "smoke" => Dict(
                    "label" => "smoke",
                    "heavy" => false,
                    "args" => ["phase-guided-transport", "--mode", "fixed-T-sparse-muB", "--case-name", "frontend_dispatch_smoke", "--outdir", "{job_output}/dispatcher_phase_guided_smoke", "--dry-run", "--overwrite"],
                ),
                "canonical" => Dict(
                    "label" => "canonical",
                    "heavy" => true,
                    "args" => ["gap-transport", "--output", "{job_output}/dispatcher_gap_transport_canonical.csv", "--failed-points-output", "{job_output}/dispatcher_gap_transport_failed_points.csv", "--provenance-dir", "{job_output}", "--tmin", "120", "--tmax", "200", "--tstep", "20", "--mubmin", "0", "--mubmax", "400", "--mubstep", "200", "--xi-list", "0.0,0.2", "--mode", "finite_15", "--overwrite"],
                ),
                "custom" => Dict("label" => "custom", "heavy" => true, "args" => String[]),
            ),
        ),
        Dict(
            "id" => "run-aniso-phase-template",
            "script" => "scripts/pnjl/run_aniso_phase_template.jl",
            "name" => "各向异性 PNJL 相图模板",
            "purpose" => "一键执行各向异性 PNJL T-mu 扫描、相结构计算，并可选生成相图。",
            "use_cases" => ["各向异性相图 smoke", "生成 experiment_templates/aniso_phase 资产", "跳过绘图只验证扫描/相结构"],
            "key_params" => ["profile=smoke|full", "xi-values", "tag", "skip-plot", "python"],
            "default_preset" => "smoke",
            "estimated_time" => Dict("smoke" => "数分钟到十几分钟", "canonical" => "较长", "custom" => "取决于 profile/xi 数量"),
            "local_recommendation" => "前端 smoke 默认跳过绘图；full profile 需要明确确认。",
            "outputs" => ["scan CSV", "phase CSV/summary", "manifest.txt", "可选 figures"],
            "manifest_candidates" => ["manifest.txt"],
            "presets" => Dict(
                "smoke" => Dict(
                    "label" => "smoke",
                    "heavy" => false,
                    "args" => ["--profile=smoke", "--xi-values=0.0", "--tag={job_id}", "--skip-plot"],
                ),
                "canonical" => Dict(
                    "label" => "canonical",
                    "heavy" => true,
                    "args" => ["--profile=full", "--xi-values=0.0,0.2,0.4", "--tag={job_id}", "--skip-plot"],
                ),
                "custom" => Dict("label" => "custom", "heavy" => true, "args" => String[]),
            ),
        ),
        Dict(
            "id" => "run-magnetic-point",
            "script" => "scripts/pnjl/run_magnetic_point.jl",
            "name" => "PNJL 磁场单点",
            "purpose" => "执行默认 T/mu/eB 下的磁场单点诊断，主要通过 stdout 查看 omega、pressure 与密度。",
            "use_cases" => ["磁场单点 smoke", "确认 magnetic 模块能加载", "查看默认点数值摘要"],
            "key_params" => ["当前脚本无 CLI 参数；使用脚本内默认 T/mu/eB"],
            "default_preset" => "smoke",
            "estimated_time" => Dict("smoke" => "通常较快"),
            "local_recommendation" => "该入口是专题诊断，不建议作为批量扫描替代。",
            "outputs" => ["stdout 诊断", "stderr 日志"],
            "manifest_candidates" => String[],
            "presets" => Dict(
                "smoke" => Dict("label" => "smoke", "heavy" => false, "args" => String[]),
            ),
        ),
        Dict(
            "id" => "run-magnetic-eb-scan",
            "script" => "scripts/pnjl/run_magnetic_eb_scan.jl",
            "name" => "PNJL 磁场 eB 扫描",
            "purpose" => "按脚本默认 eB 网格输出磁场热力学与收敛性 CSV。",
            "use_cases" => ["磁场 eB 轴 smoke", "生成 pnjl_magnetic_eb_scan.csv", "检查 nmax 收敛列"],
            "key_params" => ["当前脚本无 CLI 参数；使用脚本默认 T/mu/eB/points/output"],
            "default_preset" => "smoke",
            "estimated_time" => Dict("smoke" => "通常数分钟内", "canonical" => "同 smoke，脚本当前无 CLI 扩展"),
            "local_recommendation" => "当前输出路径由脚本默认决定；若要隔离到 job 目录，需后续给脚本补 CLI 参数。",
            "outputs" => ["pnjl_magnetic_eb_scan.csv", "stdout/stderr 日志"],
            "manifest_candidates" => String[],
            "presets" => Dict(
                "smoke" => Dict("label" => "smoke", "heavy" => false, "args" => String[]),
            ),
        ),
        Dict(
            "id" => "run-magnetic-stability-scan",
            "script" => "scripts/pnjl/run_magnetic_stability_scan.jl",
            "name" => "PNJL 磁场稳定性扫描",
            "purpose" => "扫描默认 T/mu/eB 组合并输出 pass/fail 与失败点 CSV，用于磁场专题稳定性检查。",
            "use_cases" => ["磁场稳定性 smoke", "生成 stability/failures CSV", "排查 nmax_not_converged"],
            "key_params" => ["当前脚本无 CLI 参数；使用脚本默认网格/output/failures_output"],
            "default_preset" => "smoke",
            "estimated_time" => Dict("smoke" => "数分钟到十几分钟", "custom" => "脚本当前不消费自定义 CLI 参数"),
            "local_recommendation" => "这是专题检查入口；运行前应确认默认输出可被覆盖。",
            "outputs" => ["stability scan CSV", "failures CSV", "stdout/stderr 日志"],
            "manifest_candidates" => String[],
            "presets" => Dict(
                "smoke" => Dict("label" => "smoke", "heavy" => false, "args" => String[]),
            ),
        ),
        Dict(
            "id" => "run-mott-phase-derived-csv",
            "script" => "scripts/relaxtime/run_mott_phase_derived_csv.jl",
            "name" => "Mott 派生 CSV",
            "purpose" => "读取 Mott phase scan CSV 并追加派生质量组合列，供 plot modes 使用。",
            "use_cases" => ["后处理已有 Mott CSV", "生成 M_u_plus_M_d/M_u_plus_M_s", "为 Mott 图像模式准备输入"],
            "key_params" => ["--in", "--out"],
            "default_preset" => "custom",
            "estimated_time" => Dict("custom" => "通常较快，取决于输入 CSV 行数"),
            "local_recommendation" => "需要已有 scan CSV；前端仅提供 custom 参数入口。",
            "outputs" => ["derived CSV", "stdout/stderr 日志"],
            "requires_existing_input" => true,
            "manifest_candidates" => String[],
            "presets" => Dict(
                "custom" => Dict("label" => "custom", "heavy" => true, "args" => String[]),
            ),
        ),
        Dict(
            "id" => "run-mott-phase-plot-modes",
            "script" => "scripts/relaxtime/run_mott_phase_plot_modes.jl",
            "name" => "Mott 图像模式生成",
            "purpose" => "读取 derived CSV 并生成 mode_a/mode_b/mode_ab 图像资产。",
            "use_cases" => ["后处理 Mott derived CSV", "生成 mode review figures", "快速检查 xi 分组图"],
            "key_params" => ["--in", "--out-dir", "--include-mixed"],
            "default_preset" => "custom",
            "estimated_time" => Dict("custom" => "通常数十秒到数分钟，取决于图数和 Python 环境"),
            "local_recommendation" => "需要已有 derived CSV 与可用 Python 绘图依赖；前端仅提供 custom 参数入口。",
            "outputs" => ["PNG figures under mode_a/mode_b/mode_ab", "stdout/stderr 日志"],
            "requires_existing_input" => true,
            "manifest_candidates" => String[],
            "presets" => Dict(
                "custom" => Dict("label" => "custom", "heavy" => true, "args" => String[]),
            ),
        ),
        Dict(
            "id" => "run-offline-transport-patch",
            "script" => "scripts/relaxtime/run_offline_transport_patch.jl",
            "name" => "离线 Transport 补点",
            "purpose" => "读取带 quality_flag 的 transport scan CSV，选择异常点并逐点重算输出 patch CSV，不自动回写原文件。",
            "use_cases" => ["修补质量异常点", "限制 max-points 做 Top-K 复算", "生成 gap_transport_scan_patch.csv"],
            "key_params" => ["--input", "--output", "--max-points", "--reason-filter", "tau/tr 节点参数", "--compute-bulk"],
            "default_preset" => "custom",
            "estimated_time" => Dict("custom" => "取决于异常点数量和 tau 参数；可能很长"),
            "local_recommendation" => "必须指定已有 scan CSV；建议 custom 里先设 --max-points 1 做试跑。",
            "outputs" => ["patch CSV", "stdout/stderr 日志"],
            "requires_existing_input" => true,
            "manifest_candidates" => String[],
            "presets" => Dict(
                "custom" => Dict("label" => "custom", "heavy" => true, "args" => String[]),
            ),
        ),
    ]
end

const SCRIPT_TASK_CATALOG = vcat(_script_task_catalog_items(), _extended_script_task_catalog_items())

function _script_task_by_id(task_id::AbstractString)
    for task in SCRIPT_TASK_CATALOG
        String(task["id"]) == String(task_id) && return task
    end
    return nothing
end

function _public_script_task_catalog()
    return Dict(
        "status" => "ok",
        "tasks" => SCRIPT_TASK_CATALOG,
        "template_tiers" => ["smoke", "canonical", "custom"],
        "default_allowed_tier" => "smoke",
        "heavy_requires_confirmation" => true,
        "non_default_frontend_policy" => Dict(
            "scripts/dev" => "read-only diagnostics; not default runnable tasks",
            "scripts/analysis" => "read-only diagnostics; not default runnable tasks",
            "scripts/debug" => "read-only diagnostics; not default runnable tasks",
            "scripts/perf" => "read-only diagnostics; not default runnable tasks",
        ),
    )
end

function handle_script_task_catalog(req::HTTP.Request)
    req.method == "GET" || return HTTP.Response(405, ["Content-Type" => "text/plain"], "Method Not Allowed")
    return _script_task_json_response(200, _public_script_task_catalog())
end

@inline function _script_task_json_response(status::Int, payload::AbstractDict)
    headers = [
        "Content-Type" => "application/json; charset=utf-8",
        "Access-Control-Allow-Origin" => "*",
    ]
    return HTTP.Response(status, headers, JSON3.write(_pnjl_json_safe_value(payload)))
end

@inline function _script_task_error_response(status::Int, code::String, message::String; extras::AbstractDict=Dict{String, Any}())
    return _script_task_json_response(status, _error_payload(code, message; extras=extras))
end

function _script_task_custom_args(raw)
    raw === nothing && return String[]
    raw isa AbstractVector || throw(ArgumentError("custom_args must be an array of strings"))
    length(raw) <= 80 || throw(ArgumentError("custom_args is too long"))
    args = String[]
    for item in raw
        item isa AbstractString || throw(ArgumentError("custom_args must contain strings only"))
        s = String(item)
        (contains(s, '\n') || contains(s, '\r')) && throw(ArgumentError("custom_args cannot contain newlines"))
        length(s) <= 300 || throw(ArgumentError("custom_args item too long"))
        push!(args, s)
    end
    return args
end

function _script_task_template_args(task::Dict{String, Any}, preset::String, job_id::String, custom_args::Vector{String})
    presets = task["presets"]::Dict
    haskey(presets, preset) || throw(ArgumentError("unknown preset: $(preset)"))
    preset_data = presets[preset]::Dict
    args = preset == "custom" ? custom_args : String.(preset_data["args"])
    output_dir = joinpath(_SCRIPT_TASK_JOB_ROOT, job_id)
    job_output = replace(relpath(output_dir, _PROJECT_ROOT), '\\' => '/')
    return [replace(String(arg), "{job_id}" => job_id, "{job_output}" => job_output) for arg in args]
end

function _script_task_resolve_request(body_dict::Dict{Symbol, Any}, job_id::String)
    task_id = String(get(body_dict, :task_id, get(body_dict, :id, "")))
    isempty(task_id) && throw(ArgumentError("task_id is required"))
    task = _script_task_by_id(task_id)
    task === nothing && throw(ArgumentError("unknown task_id: $(task_id)"))

    preset = lowercase(String(get(body_dict, :preset, get(task, "default_preset", "smoke"))))
    custom_args = _script_task_custom_args(get(body_dict, :custom_args, nothing))
    preset == "custom" && isempty(custom_args) && throw(ArgumentError("custom preset requires custom_args"))

    presets = task["presets"]::Dict
    haskey(presets, preset) || throw(ArgumentError("unknown preset: $(preset)"))
    preset_data = presets[preset]::Dict
    heavy = Bool(get(preset_data, "heavy", false))
    confirm_heavy = _to_bool(get(body_dict, :confirm_heavy, false), false)
    if heavy && !confirm_heavy
        throw(ArgumentError("preset $(preset) is heavy and requires confirm_heavy=true"))
    end

    args = _script_task_template_args(task, preset, job_id, custom_args)
    return (
        task=task,
        task_id=task_id,
        preset=preset,
        args=args,
        custom_args=custom_args,
        heavy=heavy,
    )
end

function _script_task_command_parts(task::Dict{String, Any}, args::Vector{String})
    script = String(task["script"])
    if Sys.iswindows()
        wrapper = joinpath("scripts", "dev", "run_with_sysimage.ps1")
        if isfile(joinpath(_PROJECT_ROOT, wrapper))
            return String["powershell", "-ExecutionPolicy", "Bypass", "-File", wrapper, "-MismatchPolicy", "fallback", script, args...], "run_with_sysimage.ps1"
        end
    else
        wrapper = joinpath("scripts", "dev", "run_with_sysimage.sh")
        if isfile(joinpath(_PROJECT_ROOT, wrapper))
            return String["sh", wrapper, "--mismatch-policy=fallback", script, args...], "run_with_sysimage.sh"
        end
    end
    julia_bin = joinpath(Sys.BINDIR, Base.julia_exename())
    return String[julia_bin, "--project=.", script, args...], "julia --project fallback"
end

function _script_task_command_preview(task::Dict{String, Any}, args::Vector{String})
    parts, wrapper_kind = _script_task_command_parts(task, args)
    return Dict(
        "wrapper" => wrapper_kind,
        "argv" => parts,
        "display" => join(map(x -> occursin(' ', x) ? "\"$(x)\"" : x, parts), " "),
    )
end

function _script_task_queue_snapshot()
    lock(_SCRIPT_TASK_JOBS_LOCK) do
        queued = count(job -> get(job, "status", "") == _SCRIPT_TASK_STATUS_QUEUED, values(_SCRIPT_TASK_JOBS))
        running = count(job -> get(job, "status", "") == _SCRIPT_TASK_STATUS_RUNNING, values(_SCRIPT_TASK_JOBS))
        return (queued=queued, running=running)
    end
end

function _script_task_queue_position(job_id::String)
    lock(_SCRIPT_TASK_JOBS_LOCK) do
        haskey(_SCRIPT_TASK_JOBS, job_id) || return nothing
        job = _SCRIPT_TASK_JOBS[job_id]
        get(job, "status", "") == _SCRIPT_TASK_STATUS_QUEUED || return nothing
        seq = Int(get(job, "seq", typemax(Int)))
        pos = 1
        for other in values(_SCRIPT_TASK_JOBS)
            if get(other, "status", "") == _SCRIPT_TASK_STATUS_QUEUED && Int(get(other, "seq", typemax(Int))) < seq
                pos += 1
            end
        end
        return pos
    end
end

function _with_script_task_job(f::Function, job_id::String)
    lock(_SCRIPT_TASK_JOBS_LOCK) do
        haskey(_SCRIPT_TASK_JOBS, job_id) || return nothing
        return f(_SCRIPT_TASK_JOBS[job_id])
    end
end

function _new_script_task_event(job::Dict{String, Any}, event_code::String; extras::AbstractDict=Dict{String, Any}())
    event = Dict{String, Any}(
        "job_id" => get(job, "job_id", nothing),
        "task_id" => get(job, "task_id", nothing),
        "state" => get(job, "status", nothing),
        "timestamp" => _timestamp_now(),
        "event_code" => event_code,
    )
    for (k, v) in pairs(extras)
        event[string(k)] = v
    end
    return event
end

function _append_script_task_event!(job::Dict{String, Any}, event_code::String; extras::AbstractDict=Dict{String, Any}())
    events = get!(job, "events") do
        Vector{Dict{String, Any}}()
    end
    event = _new_script_task_event(job, event_code; extras=extras)
    push!(events, event)
    @info "SCRIPT_TASK_JOB_EVENT" event
    return event
end

function _set_script_task_status!(job::Dict{String, Any}, status::String)
    job["status"] = status
    status == _SCRIPT_TASK_STATUS_RUNNING && get(job, "started_at", nothing) === nothing && (job["started_at"] = _timestamp_now())
    status in _SCRIPT_TASK_TERMINAL && (job["ended_at"] = _timestamp_now())
    return nothing
end

function _new_script_task_job(job_id::String, resolved)
    mkpath(_SCRIPT_TASK_JOB_ROOT)
    mkpath(_SCRIPT_TASK_LOG_ROOT)
    seq = lock(_SCRIPT_TASK_JOBS_LOCK) do
        _SCRIPT_TASK_JOB_SEQ[] += 1
        _SCRIPT_TASK_JOB_SEQ[]
    end
    output_dir = joinpath(_SCRIPT_TASK_JOB_ROOT, job_id)
    stdout_path = joinpath(_SCRIPT_TASK_LOG_ROOT, "$(job_id).out.log")
    stderr_path = joinpath(_SCRIPT_TASK_LOG_ROOT, "$(job_id).err.log")
    mkpath(output_dir)
    job = Dict{String, Any}(
        "job_id" => job_id,
        "seq" => seq,
        "kind" => "script-task",
        "task_id" => resolved.task_id,
        "task_name" => resolved.task["name"],
        "preset" => resolved.preset,
        "status" => _SCRIPT_TASK_STATUS_QUEUED,
        "created_at" => _timestamp_now(),
        "started_at" => nothing,
        "ended_at" => nothing,
        "progress" => Dict("total" => 1, "completed" => 0, "percent" => 0.0),
        "events" => Vector{Dict{String, Any}}(),
        "output_dir" => output_dir,
        "stdout_path" => stdout_path,
        "stderr_path" => stderr_path,
        "args" => resolved.args,
        "heavy" => resolved.heavy,
        "command" => _script_task_command_preview(resolved.task, resolved.args),
        "result" => nothing,
        "error" => nothing,
        "reason_code" => nothing,
        "process" => nothing,
    )
    _append_script_task_event!(job, "created")
    lock(_SCRIPT_TASK_JOBS_LOCK) do
        _SCRIPT_TASK_JOBS[job_id] = job
    end
    return job_id
end

function _tail_file(path::String, max_chars::Int=8000)
    isfile(path) || return ""
    data = read(path, String)
    length(data) <= max_chars && return data
    return last(data, max_chars)
end

function _script_task_artifacts(job::Dict{String, Any})
    output_dir = String(job["output_dir"])
    artifacts = Dict{String, Any}[]
    if isdir(output_dir)
        for (root, _dirs, files) in walkdir(output_dir)
            for file in files
                path = joinpath(root, file)
                push!(artifacts, Dict(
                    "path" => path,
                    "relative_path" => replace(relpath(path, _PROJECT_ROOT), '\\' => '/'),
                    "bytes" => filesize(path),
                    "mtime" => string(Dates.unix2datetime(mtime(path))),
                ))
            end
        end
        sort!(artifacts; by=item -> String(item["relative_path"]))
    end
    return artifacts
end

function _script_task_manifest_paths(job::Dict{String, Any})
    task = _script_task_by_id(String(job["task_id"]))
    task === nothing && return String[]
    candidates = String.(get(task, "manifest_candidates", String[]))
    isempty(candidates) && return String[]
    paths = String[]
    output_dir = String(job["output_dir"])
    isdir(output_dir) || return paths
    for (root, _dirs, files) in walkdir(output_dir)
        for file in files
            file in candidates && push!(paths, joinpath(root, file))
        end
    end
    sort!(paths)
    return paths
end

function _build_script_task_result(job::Dict{String, Any}, exit_code)
    artifacts = _script_task_artifacts(job)
    manifests = _script_task_manifest_paths(job)
    return Dict(
        "task_id" => job["task_id"],
        "task_name" => job["task_name"],
        "preset" => job["preset"],
        "exit_code" => exit_code,
        "output_dir" => job["output_dir"],
        "stdout_path" => job["stdout_path"],
        "stderr_path" => job["stderr_path"],
        "stdout_tail" => _tail_file(String(job["stdout_path"])),
        "stderr_tail" => _tail_file(String(job["stderr_path"])),
        "artifacts" => artifacts,
        "manifest_paths" => manifests,
        "command" => job["command"],
    )
end

function _launch_script_task_job(job_id::String)
    Threads.@spawn begin
        Base.acquire(_SCRIPT_TASK_SEMAPHORE)
        try
            job_snapshot = _with_script_task_job(job_id) do job
                String(get(job, "status", "")) == _SCRIPT_TASK_STATUS_QUEUED || return nothing
                _set_script_task_status!(job, _SCRIPT_TASK_STATUS_RUNNING)
                _append_script_task_event!(job, "started")
                return Dict(
                    "task_id" => job["task_id"],
                    "args" => copy(job["args"]),
                    "stdout_path" => job["stdout_path"],
                    "stderr_path" => job["stderr_path"],
                )
            end
            job_snapshot === nothing && return

            task = _script_task_by_id(String(job_snapshot["task_id"]))
            task === nothing && error("unknown task")
            args = String.(job_snapshot["args"])
            parts, _wrapper = _script_task_command_parts(task, args)
            cmd = Cmd(Cmd(parts); dir=_PROJECT_ROOT)

            exit_code = -1
            open(String(job_snapshot["stdout_path"]), "w") do out
                open(String(job_snapshot["stderr_path"]), "w") do err
                    proc = run(pipeline(cmd, stdout=out, stderr=err); wait=false)
                    _with_script_task_job(job_id) do job
                        job["process"] = proc
                    end
                    wait(proc)
                    exit_code = proc.exitcode
                end
            end

            _with_script_task_job(job_id) do job
                job["process"] = nothing
                if String(get(job, "status", "")) == _SCRIPT_TASK_STATUS_CANCELLED
                    job["result"] = _build_script_task_result(job, exit_code)
                    return
                end
                if exit_code == 0
                    _set_script_task_status!(job, _SCRIPT_TASK_STATUS_SUCCEEDED)
                    job["progress"] = Dict("total" => 1, "completed" => 1, "percent" => 100.0)
                    job["result"] = _build_script_task_result(job, exit_code)
                    _append_script_task_event!(job, "ended"; extras=Dict("outcome" => "succeeded", "exit_code" => exit_code))
                else
                    _set_script_task_status!(job, _SCRIPT_TASK_STATUS_FAILED)
                    job["error"] = "script task failed"
                    job["reason_code"] = "EXECUTION_ERROR"
                    job["result"] = _build_script_task_result(job, exit_code)
                    _append_script_task_event!(job, "ended"; extras=Dict("outcome" => "failed", "exit_code" => exit_code))
                end
            end
        catch e
            _with_script_task_job(job_id) do job
                if !(String(get(job, "status", "")) in _SCRIPT_TASK_TERMINAL)
                    _set_script_task_status!(job, _SCRIPT_TASK_STATUS_FAILED)
                    job["error"] = sprint(showerror, e)
                    job["reason_code"] = "LAUNCH_ERROR"
                    job["result"] = _build_script_task_result(job, -1)
                    _append_script_task_event!(job, "ended"; extras=Dict("outcome" => "failed", "error_type" => string(typeof(e))))
                end
            end
        finally
            Base.release(_SCRIPT_TASK_SEMAPHORE)
        end
    end
end

function handle_script_task_job_create(req::HTTP.Request)
    req.method == "POST" || return HTTP.Response(405, ["Content-Type" => "text/plain"], "Method Not Allowed")
    try
        body = isempty(req.body) ? Dict{Symbol,Any}() : JSON3.read(String(req.body))
        body_dict = body isa Dict ? body : _to_symbol_dict(body)
        snap = _script_task_queue_snapshot()
        if snap.queued >= _SCRIPT_TASK_MAX_PENDING
            return _script_task_error_response(429, "QUEUE_FULL", "script task queue is full")
        end

        job_id = string(uuid4())
        resolved = _script_task_resolve_request(body_dict, job_id)
        job_id = _new_script_task_job(job_id, resolved)
        _launch_script_task_job(job_id)
        pos = _script_task_queue_position(job_id)
        return _script_task_json_response(202, Dict(
            "status" => "accepted",
            "job_id" => job_id,
            "kind" => "script-task",
            "task_id" => resolved.task_id,
            "preset" => resolved.preset,
            "status_url" => "/api/modules/script-tasks/jobs/$(job_id)",
            "result_url" => "/api/modules/script-tasks/jobs/$(job_id)/result",
            "queue" => Dict("position" => pos, "max_running" => _SCRIPT_TASK_MAX_RUNNING, "max_pending" => _SCRIPT_TASK_MAX_PENDING),
            "command" => _with_script_task_job(job_id) do job
                job["command"]
            end,
        ))
    catch e
        return _script_task_error_response(400, "INVALID_REQUEST", sprint(showerror, e))
    end
end

function _script_task_status_payload(job_id::String)
    pos = _script_task_queue_position(job_id)
    snap = _script_task_queue_snapshot()
    return _with_script_task_job(job_id) do job
        Dict(
            "status" => "ok",
            "job_id" => job["job_id"],
            "kind" => "script-task",
            "task_id" => job["task_id"],
            "task_name" => job["task_name"],
            "preset" => job["preset"],
            "job_status" => job["status"],
            "created_at" => job["created_at"],
            "started_at" => job["started_at"],
            "ended_at" => job["ended_at"],
            "progress" => job["progress"],
            "error" => job["status"] == _SCRIPT_TASK_STATUS_FAILED ? job["error"] : nothing,
            "reason_code" => job["reason_code"],
            "queue" => Dict("position" => pos, "queued" => snap.queued, "running" => snap.running, "max_running" => _SCRIPT_TASK_MAX_RUNNING, "max_pending" => _SCRIPT_TASK_MAX_PENDING),
            "events" => job["events"],
            "command" => job["command"],
            "logs" => Dict(
                "stdout_path" => job["stdout_path"],
                "stderr_path" => job["stderr_path"],
                "stdout_tail" => _tail_file(String(job["stdout_path"])),
                "stderr_tail" => _tail_file(String(job["stderr_path"])),
            ),
        )
    end
end

function handle_script_task_job_status(job_id::String)
    payload = _script_task_status_payload(job_id)
    payload === nothing && return _script_task_error_response(404, "JOB_NOT_FOUND", "job not found"; extras=Dict("job_id" => job_id))
    return _script_task_json_response(200, payload)
end

function handle_script_task_job_cancel(job_id::String)
    payload = _with_script_task_job(job_id) do job
        status = String(get(job, "status", ""))
        status in _SCRIPT_TASK_TERMINAL && return _error_payload("JOB_NOT_CANCELLABLE", "job already in terminal state"; extras=Dict("job_id" => job_id, "job_status" => status))
        proc = get(job, "process", nothing)
        if proc !== nothing
            try
                kill(proc)
            catch err
                @warn "failed to kill script task process" job_id=job_id exception=(err, catch_backtrace())
            end
        end
        _set_script_task_status!(job, _SCRIPT_TASK_STATUS_CANCELLED)
        job["error"] = "script task cancelled"
        job["reason_code"] = "CANCELLED"
        _append_script_task_event!(job, "ended"; extras=Dict("outcome" => "cancelled"))
        return Dict("status" => "ok", "job_id" => job_id, "kind" => "script-task", "task_id" => job["task_id"], "job_status" => job["status"])
    end
    payload === nothing && return _script_task_error_response(404, "JOB_NOT_FOUND", "job not found"; extras=Dict("job_id" => job_id))
    if payload["status"] == "error"
        return _script_task_json_response(409, payload)
    end
    return _script_task_json_response(200, payload)
end

function handle_script_task_job_result(job_id::String)
    payload = _with_script_task_job(job_id) do job
        status = String(job["status"])
        if !(status in _SCRIPT_TASK_TERMINAL)
            return _error_payload("JOB_NOT_FINISHED", "job is not finished"; extras=Dict("job_id" => job_id, "job_status" => status))
        end
        result = job["result"] === nothing ? _build_script_task_result(job, status == _SCRIPT_TASK_STATUS_CANCELLED ? -1 : 0) : job["result"]
        return Dict(
            "status" => "ok",
            "job_id" => job_id,
            "kind" => "script-task",
            "task_id" => job["task_id"],
            "job_status" => status,
            "result" => result,
        )
    end
    payload === nothing && return _script_task_error_response(404, "JOB_NOT_FOUND", "job not found"; extras=Dict("job_id" => job_id))
    if payload["status"] == "error"
        return _script_task_json_response(409, payload)
    end
    return _script_task_json_response(200, payload)
end
