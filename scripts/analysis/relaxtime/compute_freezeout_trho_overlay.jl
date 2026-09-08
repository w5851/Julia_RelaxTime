"""Map saved quark-only BQS backgrounds to a historical T-rho plotting plane."""
module FreezeoutTRhoOverlay
using CSV,JSON3,SHA
const ROOT=normpath(joinpath(@__DIR__,"..","..",".."))
sha(p)=bytes2hex(sha256(read(p)))

function coordinates(rho,rho0)
    length(rho)==3 && all(isfinite,rho) && isfinite(rho0) && rho0>0 ||
        throw(ArgumentError("three finite net densities and positive rho0 required"))
    u,d,s=rho
    b=(u+d+s)/3
    return (rho_B=b,rho_Q=(2u-d-s)/3,rho_S=-s,rho_norm=b/rho0)
end

function coverage(T,rho)
    isfinite(T) && isfinite(rho) || return "nonfinite"
    reasons=String[]
    120<=T<=220 || push!(reasons,"temperature_outside")
    0.05<=rho<=1.0 || push!(reasons,"density_outside")
    return isempty(reasons) ? "inside_historical_grid" : join(reasons,';')
end

function main()
    isdefined(Main,:Models) || Base.include(Main,joinpath(ROOT,"src","models","Models.jl"))
    return Base.invokelatest(run_projection)
end

function run_projection()
    M=Main.Models
    base=joinpath(ROOT,"data","outputs","results","relaxtime","analysis","charged_rpa_phase_backend")
    input=joinpath(base,"fig4_like_freezeout_ratio_dense_20260905_v3")
    output=get(ENV,"FREEZEOUT_TRHO_OUTPUT",joinpath(base,"freezeout_on_historical_trho_20260905"))
    ispath(output) && error("refusing to overwrite $(output)")
    manifest=JSON3.read(read(joinpath(input,"manifest.json"),String))
    files=("backgrounds.csv","direct_densities.csv","manifest.json")
    input_hashes=Dict(f=>sha(joinpath(input,f)) for f in files)
    for f in files[1:2]
        input_hashes[f]==String(manifest.output_hashes[Symbol(f)]) || error("input hash mismatch: $(f)")
    end
    # The background-to-density kernel must match the saved run. Meson-only
    # changes are irrelevant to this projection and never enter the calculation.
    for (f,h) in pairs(manifest.source_hashes)
        p=replace(String(f),'\\'=>'/')
        if startswith(p,"src/models/") || p=="src/Constants_PNJL.jl" || startswith(p,"config/models/pnjl/")
            sha(joinpath(ROOT,p))==String(h) || error("upstream source changed: $(p)")
        end
    end
    model=M.create_model(:PNJL)
    rho0=Float64(Main.Constants_PNJL.ρ0_inv_fm3)
    hbarc=Float64(Main.Constants_PNJL.ħc_MeV_fm)
    backgrounds=collect(CSV.File(joinpath(input,"backgrounds.csv")))
    densities=collect(CSV.File(joinpath(input,"direct_densities.csv")))
    rows=NamedTuple[]
    for e in sort(unique(r.sqrt_s_NN_GeV for r in densities))
        T=first(filter(r->r.sqrt_s_NN_GeV==e,densities)).T_MeV
        bg=only(filter(r->r.T_MeV==T && r.p_nodes==48,backgrounds))
        state=[bg.phi_u,bg.phi_d,bg.phi_s,bg.Phi,bg.PhiBar]
        mu=[bg.mu_u,bg.mu_d,bg.mu_s]
        r48=M.model_rho(model,state,mu,T/hbarc;p_num=48,t_num=8,xi=0.)
        r96=M.model_rho(model,state,mu,T/hbarc;p_num=96,t_num=8,xi=0.)
        r192=M.model_rho(model,state,mu,T/hbarc;p_num=192,t_num=8,xi=0.)
        nd=M.number_densities(model,state,T/hbarc,mu;p_num=192,t_num=8,xi=0.)
        c48,c96,c192=coordinates(r48,rho0),coordinates(r96,rho0),coordinates(r192,rho0)
        residual=max(abs(c48.rho_Q-0.4c48.rho_B),abs(c48.rho_S))
        oracle_error=maximum(abs,r192-(nd.quark-nd.antiquark))
        drift=abs(c192.rho_norm-c48.rho_norm)/max(abs(c192.rho_norm),eps())
        passed=residual<1e-7 && oracle_error<1e-10 && drift<1e-3
        push!(rows,(sqrt_s_NN_GeV=e,T_MeV=T,muB_MeV=bg.muB_MeV,
            rho_u_fm3=r48[1],rho_d_fm3=r48[2],rho_s_fm3=r48[3],
            rho_B_fm3=c48.rho_B,rho_Q_fm3=c48.rho_Q,rho_S_fm3=c48.rho_S,
            rho0_fm3=rho0,rho_norm=c48.rho_norm,rho_norm_p96=c96.rho_norm,rho_norm_p192=c192.rho_norm,
            relative_p48_p192_change=drift,absolute_p96_p192_change=abs(c192.rho_norm-c96.rho_norm),
            distribution_oracle_error_fm3=oracle_error,bqs_residual_fm3=residual,
            coverage=coverage(T,c48.rho_norm),passed=passed,solver_called=false,production_authorized=false))
    end
    mkpath(output)
    sources=Dict{String,String}()
    for dir in ("src","config"),(folder,_,names) in walkdir(joinpath(ROOT,dir)),file in names
        endswith(file,".jl") || endswith(file,".toml") || continue
        p=joinpath(folder,file);rel=relpath(p,ROOT)
        sources[rel]=sha(p)
        target=joinpath(output,"source_snapshot",rel);mkpath(dirname(target));cp(p,target)
    end
    cp(@__FILE__,joinpath(output,"source_snapshot",basename(@__FILE__)))
    CSV.write(joinpath(output,"freezeout_trho_coordinates.csv"),rows)
    all(sha(joinpath(input,f))==h for (f,h) in input_hashes) || error("input changed")
    report=Dict("status"=>"diagnostic_background_projection","input_directory"=>input,"input_hashes"=>input_hashes,
        "source_hashes"=>sources,"script_sha256"=>sha(@__FILE__),"source_script"=>basename(@__FILE__),
        "coordinate_nodes"=>48,"check_nodes"=>[96,192],"angle_nodes"=>8,"rho0_fm3"=>rho0,
        "row_count"=>length(rows),"all_checks_passed"=>all(r.passed for r in rows),
        "solver_called"=>false,"meson_integrals_called"=>false,"production_authorized"=>false,
        "background"=>"quark-only BQS: rhoQ/rhoB=0.4,rhoS=0; saved p48 states",
        "coordinate"=>"net baryon density (rho_u+rho_d+rho_s)/(3rho0), not meson density",
        "output_hashes"=>Dict("freezeout_trho_coordinates.csv"=>sha(joinpath(output,"freezeout_trho_coordinates.csv"))))
    open(joinpath(output,"manifest.json"),"w") do io;JSON3.write(io,report);end
    println(rows)
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
end
