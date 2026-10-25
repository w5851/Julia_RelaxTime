if !isdefined(Main,:run_phase_guided_scan)
    include(joinpath(@__DIR__,"..","..","..","scripts","relaxtime","run_phase_guided_transport_scan.jl"))
end

module NativeTransportRegression
using Test, CSV, JSON3, SHA

const ROOT=normpath(joinpath(@__DIR__,"..","..",".."))
# The v3 production baseline uses full-RS effective couplings; v1/v2 and their
# qualified candidates and independent references remain historical evidence.
const BASELINE=joinpath(ROOT,"tests","baselines","relaxtime","baseline_thermal_weighted_transport_v3.csv")

function _model_config_relative(path,root)
    normalized=replace(String(path),'\\'=>'/')
    prefix=rstrip(replace(String(root),'\\'=>'/'),'/')*"/"
    startswith(normalized,prefix) || return nothing
    relative=chopprefix(normalized,prefix)
    any(part->part in ("",".",".."),split(relative,'/')) && return nothing
    return relative
end

"""Verify read-only replay across checkouts, retaining exact numerical identity."""
function replay_identity_matches(saved,identity;recorded_root,current_root=ROOT)
    hasproperty(saved,:identity) && hasproperty(saved,:identity_sha256) || return false
    hasproperty(saved.identity,:model_parameters) && hasproperty(identity,:model_parameters) || return false
    hasproperty(saved.identity.model_parameters,:path) && hasproperty(identity.model_parameters,:path) || return false
    old_path=String(saved.identity.model_parameters.path)
    old_relative=_model_config_relative(old_path,recorded_root)
    current_relative=_model_config_relative(identity.model_parameters.path,current_root)
    old_relative!==nothing && old_relative==current_relative || return false
    startswith(current_relative,"config/models/pnjl/") && endswith(current_relative,".toml") || return false
    count(entry->entry.path==current_relative,identity.sources)==1 || return false
    # The absolute location is provenance. Restore only that field for hashing;
    # every physical parameter, profile option and source hash must still match.
    restored=merge(identity,(model_parameters=merge(identity.model_parameters,(path=old_path,)),))
    raw=JSON3.write(restored)
    return saved.identity_sha256==bytes2hex(sha256(raw)) && JSON3.read(raw)==saved.identity
end

function verify(;full=false)
    targets=collect(CSV.File(BASELINE))
    metadata=JSON3.read(read(splitext(BASELINE)[1]*".json",String))
    @test bytes2hex(sha256(read(BASELINE)))==metadata.baseline_sha256
    profile=Main.Models.TransportSamplingProfiles.load_profile("thermal_weighted_v1")
    @test profile.source_sha256==metadata.profile_sha256
    ids=full ? unique(String.(getproperty.(targets,:point_id))) : ["muB0.0_alpha1.0_xi+0.3700"]
    replay=get(ENV,"TRANSPORT_REGRESSION_REPLAY_DIR","")
    replay_root=isempty(replay) ? nothing : String(JSON3.read(read(joinpath(replay,"manifest.json"),String)).root)
    for id in ids
        expected=filter(r->r.point_id==id,targets)
        point=first(expected)
        folder=isempty(replay) ? mktempdir() : joinpath(replay,"points",id)
        args=["--outdir",folder,"--mode","fixed-T-sparse-muB","--channel-diagnostics",
            "--T-list",string(point.T_MeV),"--muB-list",string(point.muB_MeV),"--xi-list",string(point.xi)]
        opts=Main.PhaseGuidedCLI.parse_args(args)
        path=joinpath(folder,"phase_guided_transport_scan.csv")
        scan_opts=Main._build_runtime_scan_opts(path,opts)
        if isempty(replay)
            context=Main.ProvenanceMetadata.new_run_context("native_transport_regression",args)
            Main.run_phase_guided_scan(opts,context;phase_reference=Main._load_runtime_phase_reference(opts))
        else
            # Reuse only actual production output tied to the current numerical
            # source and profile; no numerical baseline is treated as a run.
            saved=JSON3.read(read(Main.TransportSamplingScanConfig.sampling_contract_path(scan_opts),String))
            identity=Main.TransportSamplingScanConfig.calculation_identity(scan_opts)
            @test replay_identity_matches(saved,identity;recorded_root=replay_root)
        end
        actual=only(CSV.File(path;comment="#"))
        @test actual.converged
        @test (actual.T_MeV,actual.muB_MeV,actual.xi)==(point.T_MeV,point.muB_MeV,point.xi)
        @test actual.mode=="mode_b_fixed_T_sparse_muB"
        @test length(expected)==9
        for target in expected
            value=Float64(getproperty(actual,Symbol(target.observable)))
            @test isapprox(value,target.candidate;rtol=1e-6,atol=1e-8)
            @test abs(value/target.reference-1)<=0.01
        end
        sampling=filter(r->r.run_id==actual.run_id,
            collect(CSV.File(Main.TransportSamplingScanConfig.sampling_diagnostics_path(scan_opts))))
        required=Set(String.(Main.RelaxTime.RelaxationTime.REQUIRED_RATE_KEYS_FOR_TAU))
        @test length(sampling)==length(required)
        @test Set(String.(getproperty.(sampling,:channel)))==required
        @test all(r->r.converged && r.estimated_error<=profile.sampling.atol+profile.sampling.rtol*abs(r.rate),sampling)
    end
end
end # module
