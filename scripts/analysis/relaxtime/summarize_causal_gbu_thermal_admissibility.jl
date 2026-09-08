"""Read-only evidence verification; writes a new report, never rewrites a run."""
module CausalGBUThermalAdmissibilitySummary
using CSV,JSON3,SHA
const ROOT=normpath(joinpath(@__DIR__,"..","..",".."))
sha(p)=bytes2hex(sha256(read(p)))

function witness_verdict(r,tolerance)
    negative=r.k0_inv_fm>0 && r.vacuum_pair_inv_fm2==0 && r.landau_inv_fm2==0 &&
        r.thermal_pair_inv_fm2<0 && r.original_cut_inv_fm2<0 && r.hard_cut_inv_fm2==0
    numeric=0<=r.cut_error_inv_fm2<tolerance && 0<=r.node_change_inv_fm2<tolerance &&
        (r.q_inv_fm!=0 || 0<=r.analytic_q0_error_inv_fm2<tolerance)
    return (;negative,numeric)
end

function main()
    base=joinpath(ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    input=get(ENV,"GBU_THERMAL_INPUT",joinpath(base,"thermal_admissibility_v2_20260906"))
    output=get(ENV,"GBU_THERMAL_REPORT",joinpath(input,"audit"))
    ispath(output) && error("refusing to overwrite $(output)")
    m=JSON3.read(read(joinpath(input,"manifest.json"),String))
    m.all_evaluations_completed || error("incomplete run")
    for (f,h) in pairs(m.source_hashes)
        sha(joinpath(input,"source_snapshot",String(f)))==String(h) || error("source snapshot mismatch")
    end
    for (f,h) in pairs(m.output_hashes)
        sha(joinpath(input,String(f)))==String(h) || error("output mismatch")
    end
    for (f,h) in pairs(m.background.input_hashes)
        sha(joinpath(String(m.background.input_directory),String(f)))==String(h) || error("input mismatch")
    end
    for (f,h) in pairs(m.paper_hashes)
        sha(String(f))==String(h) || error("paper mismatch")
    end
    w=collect(CSV.File(joinpath(input,"spectral_witnesses.csv")))
    p=collect(CSV.File(joinpath(input,"continuous_probes.csv")))
    s=m.settings
    wk(r)=(String(r.channel),String(r.occupation),Float64(r.q_inv_fm),Float64(r.offset_inv_fm))
    expected_w=Set((String(c),String(o),Float64(q),d) for c in s.channels for o in s.occupations
        for q in s.witness_q for d in (0.05,0.25,1.))
    length(w)==length(expected_w) && Set(wk.(w))==expected_w || error("witness coverage mismatch")
    pk(r)=(String(r.channel),Float64(r.q_inv_fm),String(r.variant),String(r.region),Float64(r.eta_inv_fm))
    expected_p=Set((String(c),Float64(q),v,r,Float64(e)) for c in s.probe_channels for q in s.probe_q
        for v in ("all_hard","thermal24") for r in ("gap","unitary","tail") for e in s.etas)
    length(p)==length(expected_p) && Set(pk.(p))==expected_p || error("probe coverage mismatch")
    for r in w
        v=witness_verdict(r,s.cut_tolerance_inv_fm2)
        v.negative==r.negative_tail_witness && v.numeric==r.numerics_passed || error("witness gate mismatch")
        r.physical_full_correlator_sign_passed==!v.negative || error("sign flag mismatch")
        !r.production_authorized || error("invalid witness promotion")
    end
    for r in p
        ok=r.quadrature_converged && 0<=r.node_change_inv_fm2<s.direct_node_tolerance_inv_fm2
        ok==r.direct_numerics_passed || error("probe gate mismatch")
        r.quadrature_converged && !(0<=r.error_estimate_inv_fm2<=last(s.atols)) &&
            error("fine quadrature estimate exceeds its target")
        !r.rigorous_error_bound && !r.continuous_UHP_count_certified && !r.production_authorized ||
            error("invalid continuous certification")
        # PV tail rows verify the RPA sign directly, separately from GBU.
        if r.variant=="thermal24" && r.region=="tail" && r.eta_inv_fm==0
            r.k0_inv_fm>0 && r.direct_imag_inv_fm2<0 && r.imag_D_fm2<0 || error("tail sign mismatch")
        end
    end
    all(r.numerics_passed for r in w)==m.witness_numerics_passed || error("witness manifest mismatch")
    count(r->r.negative_tail_witness,w)==m.negative_tail_witnesses || error("sign manifest mismatch")
    all(r.direct_numerics_passed for r in p)==m.direct_probe_numerics_passed || error("probe manifest mismatch")
    !m.full_physical_correlator_interpretation_accepted && !m.continuous_UHP_stability_certified &&
        !m.Mott_Levinson_accepted && !m.solver_called && !m.density_computed && !m.production_authorized ||
        error("invalid run promotion")
    report=Dict("status"=>"thermal_integrity_and_reduction_checked","script_sha256"=>sha(@__FILE__),
        "input_manifest_sha256"=>sha(joinpath(input,"manifest.json")),
        "source_snapshot_count"=>length(m.source_hashes),"input_hash_count"=>length(m.background.input_hashes),
        "paper_hash_count"=>length(m.paper_hashes),"output_hash_count"=>length(m.output_hashes),
        "witness_rows"=>length(w),"negative_witnesses"=>count(r->r.negative_tail_witness,w),
        "probe_rows"=>length(p),"direct_numerics_passed"=>count(r->r.direct_numerics_passed,p),
        "max_cut_error_inv_fm2"=>maximum(r.cut_error_inv_fm2 for r in w),
        "max_q0_analytic_error_inv_fm2"=>maximum(r.analytic_q0_error_inv_fm2 for r in w if r.q_inv_fm==0),
        "max_direct_node_change_inv_fm2"=>maximum(r.node_change_inv_fm2 for r in p),
        "max_error_estimate_inv_fm2"=>maximum(r.error_estimate_inv_fm2 for r in p),
        "max_mesh_errors_inv_fm2"=>[maximum(getproperty(r,k) for r in p) for k in
            (:mesh128_error_inv_fm2,:mesh256_error_inv_fm2,:mesh512_error_inv_fm2)],
        "max_RPA_identity_error_fm2"=>maximum(r.rpa_identity_error_fm2 for r in p),
        "max_contact_identity_error_inv_fm2"=>maximum(r.contact_algebra_error_inv_fm2 for r in p),
        "full_continuum_stability_certified"=>false,"production_authorized"=>false,
        "limitation"=>"Hash/reduction checks are not a new continuous-kernel or physical proof")
    mkpath(output)
    open(joinpath(output,"manifest.json"),"w") do io
        JSON3.write(io,report)
    end
    println(JSON3.write(report))
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
