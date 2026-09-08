"""Read-only source verification and component audit of a completed dense scan."""
module CausalGBUDenseSummary
include("audit_causal_gbu_freezeout_comparison.jl")
const A=CausalGBUFreezeoutComparison
const R=A.R
using CSV,JSON3

function main()
    base=joinpath(R.ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    input=get(ENV,"GBU_DENSE_INPUT",joinpath(base,"fig4_like_freezeout_ratio_dense_20260905_v3"))
    output=get(ENV,"GBU_DENSE_SUMMARY",joinpath(base,"fig4_like_freezeout_dense_audit_20260905"))
    ispath(output) && error("refusing to overwrite $(output)")
    m=JSON3.read(read(joinpath(input,"manifest.json"),String))
    m.status=="diagnostic_completed" || error("incomplete run")
    for (f,h) in pairs(m.source_hashes)
        R.hashfile(joinpath(input,"source_snapshot",String(f)))==String(h) || error("source snapshot mismatch: $(f)")
    end
    for (f,h) in pairs(m.output_hashes)
        R.hashfile(joinpath(input,String(f)))==String(h) || error("output mismatch: $(f)")
    end
    for (f,h) in pairs(m.input_hashes)
        R.hashfile(String(f))==String(h) || error("retained input mismatch: $(f)")
    end
    densities=vcat(collect(CSV.File(joinpath(input,"direct_densities.csv"))),collect(CSV.File(joinpath(input,"reference_densities.csv"))))
    shells=vcat(collect(CSV.File(joinpath(input,"direct_shells.csv"))),collect(CSV.File(joinpath(input,"reference_shells.csv"))))
    components=A.aggregate_components(shells;qmax=Float64(m.settings.qmax),nq=Int(m.settings.nq))
    length(densities)==length(components)==8length(m.energies_GeV) || error("incomplete density grid")
    maxerr=0.0
    for r in components
        d=only(filter(d->(d.sqrt_s_NN_GeV,d.channel,d.route)==(r.sqrt_s_NN_GeV,r.channel,r.route),densities))
        maxerr=max(maxerr,abs(d.density_inv_fm3-r.total_density_inv_fm3))
    end
    maxerr<1e-12 || error("component closure failed: $(maxerr)")
    factors=NamedTuple[]
    for e in Float64.(m.energies_GeV),channel in String.(R.CHANNELS)
        d=only(filter(d->d.sqrt_s_NN_GeV==e && d.channel==channel && d.route=="direct_finite_q",densities))
        r=only(filter(d->d.sqrt_s_NN_GeV==e && d.channel==channel && d.route=="q0_lambda_reference",densities))
        push!(factors,(sqrt_s_NN_GeV=e,channel=channel,reference_over_direct=r.density_inv_fm3/d.density_inv_fm3,
            passed=Bool(d.passed) && Bool(r.passed),production_authorized=false))
    end
    mkpath(output)
    CSV.write(joinpath(output,"component_totals.csv"),components)
    CSV.write(joinpath(output,"density_factors.csv"),factors)
    report=Dict("status"=>"diagnostic_data_audit","input_directory"=>input,"input_manifest_sha256"=>R.hashfile(joinpath(input,"manifest.json")),
        "source_snapshot_count"=>length(m.source_hashes),"input_hash_count"=>length(m.input_hashes),
        "density_rows"=>length(densities),"shell_rows"=>length(shells),"component_closure_max_abs"=>maxerr,
        "failed_density_rows"=>count(r->!r.passed,densities),"failed_shell_rows"=>count(r->!r.passed,shells),
        "negative_shell_rows_retained"=>count(r->r.shell_inv_fm2<0,shells),
        "max_abs_levinson_residual"=>maximum(abs(r.levinson_residual) for r in shells),
        "max_threshold_limit_change"=>maximum(r.threshold_limit_change for r in shells),
        "source_script_sha256"=>R.hashfile(@__FILE__),"production_authorized"=>false)
    open(joinpath(output,"manifest.json"),"w") do io
        JSON3.write(io,report)
    end
    println(JSON3.write(report))
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
