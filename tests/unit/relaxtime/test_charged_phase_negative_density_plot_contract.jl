using Test

const _CHARGED_PHASE_PLOT_SCRIPT = joinpath(
    @__DIR__, "..", "..", "..", "scripts", "analysis", "relaxtime",
    "plot_charged_phase_negative_density_fig2_like.jl",
)

@testset "charged phase negative-density plot contract" begin
    @test isfile(_CHARGED_PHASE_PLOT_SCRIPT)
    source = read(_CHARGED_PHASE_PLOT_SCRIPT, String)
    @test Meta.parseall(source) isa Expr
    @test occursin("phase_object=:inverse_propagator", source)
    @test occursin("phase_sign=-1", source)
    @test occursin("strict_phase_profile", source)
    @test occursin("ordered_pv_cut", source)
    @test occursin("ordered_retarded", source)
    @test occursin("raw_phase", source)
    @test occursin("unwrapped_phase", source)
    @test occursin("anchored_phase", source)
    @test occursin("fold_0_pi", source)
    @test occursin("phase_display_fold_0_pi", source)
    @test occursin("gbu_display_fold_0_pi", source)
    @test occursin("delta - sin(2*delta)/2", source)
    @test occursin("current_integrand", source)
    @test occursin("current_integrand_raw", source)
    @test occursin("gbu_integrand", source)
    @test occursin("q_shell_current_fm3_per_dq", source)
    @test occursin("q_shell_current_raw_fm3_per_dq", source)
    @test occursin("threshold_MeV", source)
    @test occursin("charged_phase_fig2_like.png", source)
    @test occursin("charged_phase_negative_density_attribution.png", source)
    @test occursin("diagnostic_only_not_production", source)
end

@testset "Fig2 shell algebra executes without a solver or plotting runtime" begin
    sandbox = Module(:ChargedPlotAlgebra)
    helpers = Set((:_phase_derivative,:_trapz,:_gbu_phase,:_gbu_phase_derivative,:_fold_0_pi,:_shell_summary))
    for expr in Meta.parseall(read(_CHARGED_PHASE_PLOT_SCRIPT,String)).args
        expr isa Expr || continue
        definition = expr.head === :macrocall ? expr.args[end] : expr
        definition isa Expr && definition.head in (:function, :(=)) || continue
        signature = definition.args[1]
        signature isa Expr && signature.head === :(::) && (signature=signature.args[1])
        signature isa Expr && signature.head === :call && signature.args[1] in helpers || continue
        Core.eval(sandbox,expr)
    end
    omega = [0.5,1.0,1.5,2.0]
    delta = [0.0,0.0,Float64(π),Float64(π)]
    profile = (raw_phase=delta,unwrapped_phase=delta,anchored_phase=delta)
    for T in (0.2,0.8,1.4)
        shell = sandbox._shell_summary(profile,omega,fill(2.0,4),fill(6.0,4),1.0,T)
        @test shell.current_shell_fm3_per_dq ≈ 1/π^2
        @test shell.gbu_shell_fm3_per_dq ≈ shell.current_shell_fm3_per_dq
        @test shell.current_shell_unwrapped_fm3_per_dq ≈ shell.current_shell_fm3_per_dq
    end
    shell0 = sandbox._shell_summary(profile,omega,fill(2.0,4),fill(6.0,4),0.0,0.8)
    @test shell0.current_shell_fm3_per_dq == shell0.gbu_shell_fm3_per_dq == 0
end

@testset "Frozen gap and historical comparison script contracts" begin
    directory = dirname(_CHARGED_PHASE_PLOT_SCRIPT)
    gap = read(joinpath(directory,"audit_charged_gap_spectrum.jl"),String)
    comparison = read(joinpath(directory,"compare_charged_phase_temp7.jl"),String)
    infrared = read(joinpath(directory,"audit_charged_phase_low_energy.jl"),String)
    @test Meta.parseall(infrared) isa Expr
    @test occursin("B0_spectral_cut",infrared)
    @test occursin("bu_phase_integral_parts",infrared)
    @test occursin("refusing to overwrite",infrared)
    @test occursin("source changed during audit",infrared)
    @test occursin("full_propagator_from_oracle=false",infrared)
    @test !occursin("Models.solve",infrared)
    @test Meta.parseall(gap) isa Expr
    @test Meta.parseall(comparison) isa Expr
    @test occursin("continue_gap_roots",gap)
    @test occursin("split_bu_shell",gap)
    @test occursin("inverse parity failed",gap)
    @test occursin("source changed during audit",gap)
    @test occursin("full_state_count_certified=false",gap)
    @test !occursin("Models.solve",gap)
    @test occursin("refusing to overwrite",gap)
    @test occursin("refusing to overwrite",comparison)
    @test occursin("cross_background_parity_claim",comparison)
    @test occursin("current_mu_inv_fm",comparison)
    sandbox = Module(:ChargedGapCsvContract)
    Core.eval(sandbox,:(using CSV))
    top_level = map(Meta.parseall(gap).args) do expression
        expression isa Expr && expression.head === :macrocall ? expression.args[end] : expression
    end
    definition = only(filter(top_level) do expression
        expression isa Expr && expression.head === :module
    end).args[end]
    writer = only(filter(definition.args) do expression
        expression isa Expr && expression.head === :function &&
            expression.args[1].args[1] === :_write_root_rows
    end)
    Core.eval(sandbox,writer)
    mktempdir() do directory
        path = joinpath(directory,"no_roots.csv")
        sandbox._write_root_rows(path,NamedTuple[])
        @test isfile(path)
        @test length(sandbox.CSV.File(path)) == 0
        @test :omega_inv_fm in propertynames(sandbox.CSV.File(path))
    end
end
