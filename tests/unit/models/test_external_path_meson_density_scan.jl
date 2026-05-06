using Test

const PROJECT_ROOT_EP = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT_EP, "src", "models", "Models.jl"))
end

const _EPMDS = Models.ExternalPathMesonDensityScan

@testset "ExternalPathMesonDensityScan" begin
    @testset "default output path exists as string" begin
        @test _EPMDS.DEFAULT_EXTERNAL_PATH_MESON_DENSITY_OUTPUT_PATH isa String
        @test !isempty(_EPMDS.DEFAULT_EXTERNAL_PATH_MESON_DENSITY_OUTPUT_PATH)
    end

    @testset "entrypoint exists" begin
        @test isdefined(_EPMDS, :run_external_path_meson_density_scan)
        @test _EPMDS.run_external_path_meson_density_scan isa Function
        @test Models.run_external_path_meson_density_scan isa Function
    end

    @testset "header keeps path metadata first" begin
        @test startswith(_EPMDS.HEADER, "path_source,path_case_id,path_line_style,path_point_index,path_order_key,")
        @test occursin("muB_MeV", _EPMDS.HEADER)
        @test occursin("kpi_ratio", _EPMDS.HEADER)
    end

    @testset "continuation resets on path group change" begin
        pt_a0 = (path_source="figA", path_case_id="caseA", path_line_style="solid")
        pt_a1 = (path_source="figA", path_case_id="caseA", path_line_style="solid")
        pt_b0 = (path_source="figB", path_case_id="caseB", path_line_style="dashed")

        @test _EPMDS._path_group_key(pt_a0) == _EPMDS._path_group_key(pt_a1)
        @test _EPMDS._path_group_key(pt_a0) != _EPMDS._path_group_key(pt_b0)
    end
end
