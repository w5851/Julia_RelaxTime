"""Independent integrity and gate audit of a completed regulator comparison."""
module CausalGBURegulatorSummary
using CSV,JSON3,SHA
const ROOT=normpath(joinpath(@__DIR__,"..","..",".."))
sha(path)=bytes2hex(sha256(read(path)))

function main()
    base=joinpath(ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    input=get(ENV,"GBU_CLOSURE_INPUT",joinpath(base,"regulator_routing_closure_refined_20260905"))
    output=get(ENV,"GBU_CLOSURE_SUMMARY",joinpath(input,"audit"))
    ispath(output) && error("refusing to overwrite $(output)")
    m=JSON3.read(read(joinpath(input,"manifest.json"),String))
    m.status=="diagnostic_closure_checks_completed" || error("incomplete run")
    for (f,h) in pairs(m.source_hashes)
        sha(joinpath(input,"source_snapshot",String(f)))==String(h) || error("source snapshot mismatch: $(f)")
    end
    for (f,h) in pairs(m.output_hashes)
        sha(joinpath(input,String(f)))==String(h) || error("output mismatch: $(f)")
    end
    for (f,h) in pairs(m.background.input_hashes)
        sha(joinpath(String(m.background.input_directory),String(f)))==String(h) || error("background mismatch")
    end
    rows=collect(CSV.File(joinpath(input,"closure_summary.csv")))
    profiles=collect(CSV.File(joinpath(input,"complex_profiles.csv")))
    gaps=collect(CSV.File(joinpath(input,"support_gap_counts.csv")))
    weak=collect(CSV.File(joinpath(input,"weak_eta_convergence.csv")))
    key(r)=(r.channel,r.regulator,r.q_inv_fm,r.mesh)
    expected=8length(m.configurations)
    length(unique(key.(rows)))==length(rows)==expected || error("missing or duplicate comparison rows")
    length(profiles)==5expected || error("missing profile probes")
    for r in rows
        ps=filter(p->key(p)==key(r),profiles)
        gs=filter(g->key(g)==key(r),gaps)
        length(ps)==5 && !isempty(gs) || error("missing probes or gaps")
        maximum(p.complex_difference for p in ps)==r.max_causal_difference || error("profile reduction mismatch")
        maximum(p.atom_node_change for p in ps)==r.max_atom_node_change || error("node reduction mismatch")
        sum(g.root_count for g in gs)==r.all_support_gap_count || error("root count mismatch")
        all(g.passed for g in gs)==r.all_support_gaps_passed || error("gap gate mismatch")
        causal=r.max_causal_difference<m.representation_target && r.max_atom_node_change<m.direct_node_target && r.max_contact_residual<1e-10
        causal==r.causal_passed || error("causal gate mismatch")
        passed=causal && r.normal_gap_passed && r.all_support_gaps_passed && r.normal_count==r.all_support_gap_count && r.phase_passed
        passed==r.conditional_window_passed || error("joint gate mismatch")
    end
    comparisons=NamedTuple[]
    for a in rows
        a.regulator=="two_line" || continue
        b=only(filter(b->b.channel==a.channel && b.q_inv_fm==a.q_inv_fm && b.mesh==a.mesh && b.regulator=="centered",rows))
        push!(comparisons,(channel=String(a.channel),q_inv_fm=a.q_inv_fm,mesh=a.mesh,
            two_line_root_inv_fm=a.normal_gap_root_inv_fm,centered_root_inv_fm=b.normal_gap_root_inv_fm,
            root_difference_inv_fm=a.normal_gap_root_inv_fm-b.normal_gap_root_inv_fm,
            two_line_count=a.all_support_gap_count,centered_count=b.all_support_gap_count,
            contact_difference_inv_fm2=a.contact_inv_fm2-b.contact_inv_fm2,
            both_conditionally_passed=a.conditional_window_passed && b.conditional_window_passed,production_authorized=false))
    end
    weakrows=NamedTuple[]
    for channel in unique(getproperty.(weak,:channel)),regulator in unique(getproperty.(weak,:regulator))
        rs=sort(filter(r->r.channel==channel && r.regulator==regulator,weak);by=r->-r.eta_inv_fm)
        length(rs)==6 || error("missing weak-eta probes")
        errors=[r.absolute_eta_error for r in rs]
        push!(weakrows,(channel=String(channel),regulator=String(regulator),eta_first=first(rs).eta_inv_fm,
            eta_last=last(rs).eta_inv_fm,error_first=first(errors),error_last=last(errors),
            monotone_decrease=all(<(0),diff(errors)),last_halving_ratio=errors[end-1]/errors[end],
            last_representation_difference=last(rs).representation_difference,gbu_density_convergence=false))
    end
    mkpath(output)
    CSV.write(joinpath(output,"regulator_root_comparison.csv"),comparisons)
    CSV.write(joinpath(output,"weak_eta_summary.csv"),weakrows)
    report=Dict("status"=>"diagnostic_integrity_and_gate_audit","input_directory"=>input,
        "input_manifest_sha256"=>sha(joinpath(input,"manifest.json")),"source_snapshot_count"=>length(m.source_hashes),
        "output_hash_count"=>length(m.output_hashes),"background_input_count"=>length(m.background.input_hashes),
        "comparison_rows"=>length(rows),"complex_profile_rows"=>length(profiles),"gap_rows"=>length(gaps),"weak_rows"=>length(weak),
        "conditional_passed"=>count(r->r.conditional_window_passed,rows),"failed_evaluations"=>m.failed_evaluations,
        "max_causal_difference"=>maximum(r.max_causal_difference for r in rows),
        "max_atom_node_change"=>maximum(r.max_atom_node_change for r in rows),
        "max_contact_residual"=>maximum(r.max_contact_residual for r in rows),
        "q0_root_max_difference"=>maximum(abs(r.root_difference_inv_fm) for r in comparisons if r.q_inv_fm==0),
        "weak_monotone_groups"=>count(r->r.monotone_decrease,weakrows),"density_computed"=>false,
        "production_authorized"=>false,"full_spectrum_certified"=>false,
        "script_sha256"=>sha(@__FILE__),"output_hashes"=>Dict(f=>sha(joinpath(output,f)) for f in readdir(output)))
    open(joinpath(output,"manifest.json"),"w") do io
        JSON3.write(io,report)
    end
    println(JSON3.write(report))
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
