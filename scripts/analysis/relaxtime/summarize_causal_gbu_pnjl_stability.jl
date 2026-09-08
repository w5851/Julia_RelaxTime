"""Verify retained health evidence without rerunning a solver or modifying inputs."""
module CausalGBUPNJLStabilitySummary
using CSV,JSON3,SHA
const ROOT=normpath(joinpath(@__DIR__,"..","..",".."))
sha(p)=bytes2hex(sha256(read(p)))
key(r)=(String(r.channel),Float64(r.q_inv_fm),String(r.variant))

function check_contour(c)
    finite=all(isfinite,(c.coarse_winding,c.fine_winding,c.max_step,c.minimum_inverse_abs))
    passed=finite && c.count>=0 && c.minimum_inverse_abs>1e-8 && c.max_step<pi/4 &&
        max(abs(c.coarse_winding-c.count),abs(c.fine_winding-c.count))<1e-7
    c.resolved==passed || error("contour gate mismatch")
    !c.full_UHP_certified && c.unresolved_near_axis_strip || error("invalid all-UHP claim")
    return c.resolved && c.count==0
end

function main()
    base=joinpath(ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    input=get(ENV,"GBU_HEALTH_INPUT",joinpath(base,"pnjl_routing_stability_v2_20260906"))
    output=get(ENV,"GBU_HEALTH_REPORT",joinpath(input,"audit"))
    ispath(output) && error("refusing to overwrite $(output)")
    m=JSON3.read(read(joinpath(input,"manifest.json"),String))
    m.all_evaluations_completed || error("incomplete run; inspect failures first")
    for (f,h) in pairs(m.source_hashes)
        sha(joinpath(input,"source_snapshot",String(f)))==String(h) || error("snapshot mismatch: $(f)")
    end
    for (f,h) in pairs(m.output_hashes)
        sha(joinpath(input,String(f)))==String(h) || error("output mismatch: $(f)")
    end
    for (f,h) in pairs(m.background.input_hashes)
        sha(joinpath(String(m.background.input_directory),String(f)))==String(h) || error("input mismatch")
    end
    readrows(f)=collect(CSV.File(joinpath(input,f*".csv")))
    rows,mesh,poles,contours,q0,routing=readrows.(("summary","mesh_checks","poles","contours","q0_contacts","routing"))
    channels=hasproperty(m.settings,:channels) ? String.(m.settings.channels) :
        ["pi_plus","pi_minus","K_plus","K_minus"]
    variants=hasproperty(m.settings,:variants) ? String.(m.settings.variants) : ["all_hard","thermal24"]
    expected=Set((ch,Float64(q),v) for ch in channels for q in m.settings.qs for v in variants)
    length(rows)==length(expected) && Set(key.(rows))==expected || error("case coverage mismatch")
    length(mesh)==2length(rows) && length(contours)==6length(rows) || error("mesh/contour coverage mismatch")
    for r in rows
        ms=sort(filter(x->key(x)==key(r),mesh);by=x->x.mesh)
        length(ms)==2 || error("missing mesh")
        for x in ms
            ps=filter(p->key(p)==key(r) && p.mesh==x.mesh,poles)
            cs=filter(c->key(c)==key(r) && c.mesh==x.mesh,contours)
            length(cs)==3 && Set(c.eta_floor_inv_fm for c in cs)==Set(Float64.(m.settings.etas)) ||
                error("missing contour")
            all(check_contour.(cs))==x.UHP_tested_bands_clear || error("UHP reduction mismatch")
            length(ps)==x.real_root_count || error("pole count mismatch")
            for p in ps
                (p.simple && p.k0_inv_fm*p.weight_fm>0)==p.sign_passed || error("pole sign mismatch")
                abs(p.derivative-p.finite_difference_derivative)/max(1.,abs(p.derivative))≈
                    p.derivative_relative_error || error("derivative reduction mismatch")
                if hasproperty(p,:derivative_check_passed)
                    (max(p.derivative_relative_error,p.derivative_step_change)<1e-5)==p.derivative_check_passed ||
                        error("finite difference step gate mismatch")
                end
            end
            all(p.sign_passed && p.derivative_relative_error<1e-5 &&
                (!hasproperty(p,:derivative_check_passed) || p.derivative_check_passed) for p in ps)==x.pole_weight_checks_passed ||
                error("pole check mismatch")
            (abs(x.static_cut_imag)<1e-10)==x.static_cut_passed || error("static cut mismatch")
            (x.minimum_k0_rho>=-1e-10 && abs(x.static_inverse_imag)<=1e-10 &&
                x.static_inverse_real>0)==x.numerical_passivity || error("passivity mismatch")
        end
        stable=ms[1].real_root_count==ms[2].real_root_count
        stable==r.real_root_counts_stable || error("root count stability mismatch")
        passed=stable && all(x.real_gap_checks_passed && x.pole_weight_checks_passed &&
            x.UHP_tested_bands_clear && x.static_cut_passed && x.static_inverse_real>0 for x in ms)
        passed==r.diagnostic_checks_passed || error("diagnostic reduction mismatch")
        last(ms).numerical_passivity==r.numerical_passivity || error("fine passivity mismatch")
    end
    length(q0)==2length(channels) && length(routing)==length(channels)*length(m.settings.qs) ||
        error("control coverage mismatch")
    all((r.error<1e-9)==r.passed for r in q0) || error("q0 gate mismatch")
    all((r.coordinate_error<1e-6 && r.node_change<1e-6)==r.numerics_passed &&
        (r.flavor_reflection_error<1e-9)==r.flavor_reflection_passed for r in routing) ||
        error("routing gate mismatch")
    all(r.diagnostic_checks_passed for r in rows)==m.main_diagnostic_checks_passed ||
        error("manifest diagnostic mismatch")
    all(r.numerical_passivity for r in rows)==m.all_numerical_passivity_passed ||
        error("manifest passivity mismatch")
    !m.production_authorized && !m.full_continuum_stability_certified || error("invalid promotion")
    report=Dict("status"=>"health_integrity_and_reduction_checked",
        "input_manifest_sha256"=>sha(joinpath(input,"manifest.json")),
        "script_sha256"=>sha(@__FILE__),"source_snapshot_count"=>length(m.source_hashes),
        "input_hash_count"=>length(m.background.input_hashes),"output_hash_count"=>length(m.output_hashes),
        "cases"=>length(rows),"diagnostic_passed_cases"=>count(r->r.diagnostic_checks_passed,rows),
        "pole_rows"=>length(poles),"contour_rows"=>length(contours),"mesh_rows"=>length(mesh),
        "all_tested_UHP_bands_clear"=>all(check_contour(c) for c in contours),
        "all_pole_signs_passed"=>all(p.sign_passed for p in poles),
        "root_mesh_drift_inv_fm"=>maximum(r.max_root_mesh_drift_inv_fm for r in rows),
        "max_q0_contact_error"=>maximum(r.error for r in q0),
        "max_single_sphere_coordinate_error"=>maximum(r.coordinate_error for r in routing),
        "max_single_sphere_node_change"=>maximum(r.node_change for r in routing),
        "max_derivative_error"=>maximum(p.derivative_relative_error for p in poles),
        "full_continuum_stability_certified"=>false,"production_authorized"=>false,
        "limitations"=>["Gap enclosures use the tested upstream gap auditor; raw gap rows are not repeated here",
            "Hash and reduction checks are not an independent continuous-kernel physics proof"])
    mkpath(output)
    open(joinpath(output,"manifest.json"),"w") do io
        JSON3.write(io,report)
    end
    println(JSON3.write(report))
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
