using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
Models.pnjl_module()
const PNJL = Models.pnjl_module()

function _read_single_data_line(path)
    lines = readlines(path)
    @test length(lines) == 2
    return lines[2]
end

@testset "ScanConfig equivalence smoke: TmuScan" begin
    tmp_dir = mktempdir()
    out_kw = joinpath(tmp_dir, "tmu_kw.csv")
    out_cfg = joinpath(tmp_dir, "tmu_cfg.csv")

    kw = (
        T_values=[90.0],
        mu_values=[10.0],
        xi_values=[0.0],
        overwrite=true,
        resume=false,
        p_num=10,
        t_num=4,
        iterations=60,
    )

    stats_kw = Models.run_tmu_scan(; kw..., output_path=out_kw)

    cfg = PNJL.TmuScanConfig(
        T_values=[90.0],
        mu_values=[10.0],
        xi_values=[0.0],
        output_path=out_cfg,
        overwrite=true,
        resume=false,
        p_num=10,
        t_num=4,
        nlsolve_kwargs=(iterations=60,),
    )
    stats_cfg = Models.run_tmu_scan(cfg)

    @test stats_kw.total == stats_cfg.total == 1
    @test stats_kw.success == stats_cfg.success

    line_kw = _read_single_data_line(out_kw)
    line_cfg = _read_single_data_line(out_cfg)
    @test split(line_kw, ',')[1:3] == split(line_cfg, ',')[1:3]
end

@testset "ScanConfig passthrough: semantic options" begin
    tmp_dir = mktempdir()
    out_tmu = joinpath(tmp_dir, "tmu_cfg_semantic.csv")
    out_trho = joinpath(tmp_dir, "trho_cfg_semantic.csv")

    cfg_tmu = PNJL.TmuScanConfig(
        T_values=[90.0],
        mu_values=[10.0],
        xi_values=[0.0],
        output_path=out_tmu,
        overwrite=true,
        resume=false,
        solver_backend=:auto,
        auto_pnjl_backend=:models,
        semantic_mode=:ground_state,
        p_num=10,
        t_num=4,
        nlsolve_kwargs=(iterations=60,),
    )
    stats_tmu = Models.run_tmu_scan(cfg_tmu)
    @test stats_tmu.total == 1
    @test stats_tmu.success == 1

    cfg_trho = PNJL.TrhoScanConfig(
        T_values=[90.0],
        rho_values=[0.2],
        xi_values=[0.0],
        output_path=out_trho,
        overwrite=true,
        resume=false,
        reverse_rho=false,
        solver_backend=:models,
        semantic_mode=:constrained_manifold,
        p_num=10,
        t_num=4,
        nlsolve_kwargs=(iterations=60,),
    )
    stats_trho = Models.run_trho_scan(cfg_trho)
    @test stats_trho.total == 1
    @test stats_trho.success == 1
end

@testset "ScanConfig equivalence smoke: TrhoScan" begin
    tmp_dir = mktempdir()
    out_kw = joinpath(tmp_dir, "trho_kw.csv")
    out_cfg = joinpath(tmp_dir, "trho_cfg.csv")

    kw = (
        T_values=[90.0],
        rho_values=[0.2],
        xi_values=[0.0],
        overwrite=true,
        resume=false,
        reverse_rho=false,
        p_num=10,
        t_num=4,
        iterations=60,
    )

    stats_kw = Models.run_trho_scan(; kw..., output_path=out_kw)

    cfg = PNJL.TrhoScanConfig(
        T_values=[90.0],
        rho_values=[0.2],
        xi_values=[0.0],
        output_path=out_cfg,
        overwrite=true,
        resume=false,
        reverse_rho=false,
        p_num=10,
        t_num=4,
        nlsolve_kwargs=(iterations=60,),
    )
    stats_cfg = Models.run_trho_scan(cfg)

    @test stats_kw.total == stats_cfg.total == 1
    @test stats_kw.success == stats_cfg.success

    line_kw = _read_single_data_line(out_kw)
    line_cfg = _read_single_data_line(out_cfg)
    @test split(line_kw, ',')[1:3] == split(line_cfg, ',')[1:3]
end
