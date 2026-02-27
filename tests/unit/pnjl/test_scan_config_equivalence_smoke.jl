using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
Models.legacy_pnjl_module()
const PNJL = Models.legacy_pnjl_module()

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

    stats_kw = run_tmu_scan(; kw..., output_path=out_kw)

    cfg = TmuScanConfig(
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
    stats_cfg = run_tmu_scan(cfg)

    @test stats_kw.total == stats_cfg.total == 1
    @test stats_kw.success == stats_cfg.success

    line_kw = _read_single_data_line(out_kw)
    line_cfg = _read_single_data_line(out_cfg)
    @test split(line_kw, ',')[1:3] == split(line_cfg, ',')[1:3]
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

    stats_kw = run_trho_scan(; kw..., output_path=out_kw)

    cfg = TrhoScanConfig(
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
    stats_cfg = run_trho_scan(cfg)

    @test stats_kw.total == stats_cfg.total == 1
    @test stats_kw.success == stats_cfg.success

    line_kw = _read_single_data_line(out_kw)
    line_cfg = _read_single_data_line(out_cfg)
    @test split(line_kw, ',')[1:3] == split(line_cfg, ',')[1:3]
end
