"""Read-only input verification and automatic reduction of signed-gap audits.

Writes a new report directory only. A numerical pass never selects a regulator
or authorizes production; failed evaluations and failed gates are retained.
"""
module CausalGBUFullClosureSummary
using CSV,JSON3,SHA
const ROOT=normpath(joinpath(@__DIR__,"..","..",".."))
sha(path)=bytes2hex(sha256(read(path)))
key(r)=(String(r.channel),Float64(r.q_inv_fm),String(r.regulator))

function check_case(row,profiles,meshes,gaps,cuts)
    length(profiles)==14 && length(meshes)==2 && !isempty(gaps) && !isempty(cuts) || error("missing case evidence")
    ms=sort(meshes;by=r->r.mesh)
    for m in ms
        ps=filter(p->p.mesh==m.mesh,profiles)
        gs=filter(g->g.mesh==m.mesh,gaps)
        cs=filter(c->c.mesh==m.mesh,cuts)
        length(ps)==7 && !isempty(gs) && !isempty(cs) || error("missing mesh evidence")
        maximum(p.spectral_difference for p in ps)==m.max_spectral_difference || error("spectral reduction mismatch")
        maximum(c.interpolation_error for c in cs)==m.max_direct_cut_interpolation_error || error("cut reduction mismatch")
        maximum(c.reciprocity_error for c in cs)==m.max_cut_reciprocity_error || error("reciprocity reduction mismatch")
        sum(g.root_count for g in gs)==m.signed_real_root_count || error("root count mismatch")
        all(g.passed for g in gs)==m.all_interpolant_gaps_passed || error("root gate mismatch")
        geometry=max(m.support_excess_width_inv_fm,m.support_missing_width_inv_fm)<=m.support_roundoff_allowance_inv_fm
        geometry==m.geometry_passed || error("support gate mismatch")
        m.tail_passed==(m.tail_inverse_deviation_bound<1) || error("tail gate mismatch")
    end
    maximum(p.coordinate_difference for p in profiles)==row.max_coordinate_difference || error("coordinate reduction mismatch")
    maximum(p.radial_node_difference for p in profiles)==row.max_radial_node_difference || error("node reduction mismatch")
    coordinate=row.max_coordinate_difference<1e-6 && row.max_radial_node_difference<1e-6 && row.contact_coordinate_difference<1e-8
    algebra=row.max_contact_identity_residual<1e-10 && row.max_contact_moment_residual<1e-10
    representation=last(ms).max_spectral_difference<1e-4
    stable=first(ms).signed_real_root_count==last(ms).signed_real_root_count
    stable==row.mesh_root_counts_stable || error("mesh count mismatch")
    gap=stable && all(m.all_interpolant_gaps_passed && m.geometry_passed && m.tail_passed for m in ms)
    reciprocity=last(ms).max_cut_reciprocity_error<1e-10
    cut=last(ms).max_direct_cut_interpolation_error<1e-4
    flags=(coordinate,algebra,representation,gap,reciprocity,cut,last(ms).static_passed)
    recorded=(row.coordinate_passed,row.algebra_passed,row.representation_passed,
        row.all_real_gap_checks_passed,row.reciprocity_passed,row.cut_interpolation_passed,row.static_passed)
    flags==recorded && all(flags)==row.passed || error("joint gate mismatch")
    return row.passed
end

function main()
    base=joinpath(ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    input=get(ENV,"GBU_FULL_INPUT",joinpath(base,"full_signed_regulator_closure_20260906"))
    output=get(ENV,"GBU_FULL_REPORT",joinpath(input,"audit"))
    ispath(output) && error("refusing to overwrite $(output)")
    m=JSON3.read(read(joinpath(input,"manifest.json"),String))
    m.status=="preproduction_closure_audit_completed" || error("incomplete run")
    for (f,h) in pairs(m.source_hashes)
        sha(joinpath(input,"source_snapshot",String(f)))==String(h) || error("source snapshot mismatch: $(f)")
    end
    for (f,h) in pairs(m.output_hashes)
        sha(joinpath(input,String(f)))==String(h) || error("output mismatch: $(f)")
    end
    for (f,h) in pairs(m.background.input_hashes)
        sha(joinpath(String(m.background.input_directory),String(f)))==String(h) || error("background mismatch")
    end
    rows=collect(CSV.File(joinpath(input,"summary.csv")))
    profiles=collect(CSV.File(joinpath(input,"complex_probes.csv")))
    meshes=collect(CSV.File(joinpath(input,"spectral_mesh_checks.csv")))
    gaps=collect(CSV.File(joinpath(input,"signed_gap_counts.csv")))
    cuts=collect(CSV.File(joinpath(input,"real_axis_cut_checks.csv")))
    q0=collect(CSV.File(joinpath(input,"q0_checks.csv")))
    failures=isfile(joinpath(input,"failures.csv")) ? collect(CSV.File(joinpath(input,"failures.csv"))) : []
    length(unique(key.(rows)))==length(rows)==m.completed_cases || error("duplicate or missing rows")
    length(failures)==m.failed_evaluations || error("failure count mismatch")
    length(rows)+length(failures)==m.expected_cases==8length(m.settings.qs) || error("case coverage mismatch")
    expected=Set((ch,Float64(q),reg) for ch in ("pi_plus","pi_minus","K_plus","K_minus") for q in m.settings.qs for reg in ("two_line","centered"))
    Set(vcat(key.(rows),key.(failures)))==expected || error("case identity mismatch")
    length(profiles)==14length(rows) && length(meshes)==2length(rows) || error("extra or missing probes")
    for r in rows
        check_case(r,filter(p->key(p)==key(r),profiles),filter(p->key(p)==key(r),meshes),
            filter(p->key(p)==key(r),gaps),filter(p->key(p)==key(r),cuts))
    end
    length(q0)==4 && length(unique(r.channel for r in q0))==4 || error("missing q0 checks")
    for r in q0
        (r.q0_regulator_value_difference<1e-12 && r.q0_contact_A_difference<1e-10 &&
            r.q0_all_hard_old_B0_difference<1e-8)==r.passed || error("q0 gate mismatch")
    end
    passed=isempty(failures) && all(r.passed for r in rows) && all(r.passed for r in q0)
    passed==m.all_checks_passed && count(r->r.passed,rows)==m.passed_cases || error("manifest gate mismatch")
    comparisons=NamedTuple[]
    for a in rows
        a.regulator=="two_line" || continue
        bs=filter(b->b.channel==a.channel && b.q_inv_fm==a.q_inv_fm && b.regulator=="centered",rows)
        isempty(bs) && continue
        b=only(bs)
        push!(comparisons,(channel=String(a.channel),q_inv_fm=a.q_inv_fm,
            two_line_signed_roots=a.signed_real_root_count,centered_signed_roots=b.signed_real_root_count,
            vacuum_contact_difference_inv_fm2=a.vacuum_contact_inv_fm2-b.vacuum_contact_inv_fm2,
            thermal_contact_difference_inv_fm2=a.thermal_contact_inv_fm2-b.thermal_contact_inv_fm2,
            two_line_vacuum_boost_real=a.vacuum_boost_violation_real,
            centered_vacuum_boost_real=b.vacuum_boost_violation_real,
            both_passed=a.passed && b.passed,production_authorized=false))
    end
    mkpath(output)
    CSV.write(joinpath(output,"regulator_decomposition.csv"),comparisons)
    report=Dict("status"=>"signed_gap_integrity_and_gate_report","input_directory"=>input,
        "input_manifest_sha256"=>sha(joinpath(input,"manifest.json")),"script_sha256"=>sha(@__FILE__),
        "source_snapshot_count"=>length(m.source_hashes),"output_hash_count"=>length(m.output_hashes),
        "expected_cases"=>m.expected_cases,"completed_cases"=>length(rows),"failed_evaluations"=>length(failures),
        "passed_cases"=>count(r->r.passed,rows),"q0_passed"=>all(r.passed for r in q0),"all_checks_passed"=>passed,
        "complex_probe_rows"=>length(profiles),"mesh_rows"=>length(meshes),"gap_rows"=>length(gaps),"cut_rows"=>length(cuts),
        "max_coordinate_difference"=>maximum(r.max_coordinate_difference for r in rows),
        "max_node_difference"=>maximum(r.max_radial_node_difference for r in rows),
        "max_fine_representation_difference"=>maximum(r.fine_spectral_difference for r in rows),
        "max_contact_identity_residual"=>maximum(r.max_contact_identity_residual for r in rows),
        "max_root_mesh_shift_inv_fm"=>maximum(r.max_root_mesh_shift_inv_fm for r in rows),
        "density_computed"=>false,"production_authorized"=>false,"regulator_selected"=>false,
        "full_physical_spectrum_certified"=>false,"root_scope"=>m.root_scope,
        "output_hashes"=>Dict("regulator_decomposition.csv"=>sha(joinpath(output,"regulator_decomposition.csv"))))
    open(joinpath(output,"manifest.json"),"w") do io
        JSON3.write(io,report)
    end
    println(JSON3.write(report))
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
