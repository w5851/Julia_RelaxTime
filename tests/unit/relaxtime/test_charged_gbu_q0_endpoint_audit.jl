using Test, JSON3, CSV

const Q0_ENDPOINT_PATH=joinpath(@__DIR__,"..","..","..","scripts","analysis","relaxtime",
    "audit_charged_gbu_q0_endpoint.jl")
isdefined(Main,:ChargedGBUQ0EndpointAudit) || include(Q0_ENDPOINT_PATH)
const Q0_ENDPOINT_AUDIT=Main.ChargedGBUQ0EndpointAudit

@testset "Endpoint diagnostic keeps original real/imaginary gate conditions separate" begin
    a=Q0_ENDPOINT_AUDIT
    @test a.inverse_record(1+0im).static_passed
    r=a.inverse_record(1+2e-8im)
    @test r.real_positive && !r.imaginary_passed && !r.static_passed
    r=a.inverse_record(-1+0im)
    @test !r.real_positive && r.imaginary_passed && !r.static_passed
    @test !a.inverse_record(Inf+0im).finite
    @test a.IMAG_TOL==1e-8
    @test a.inverse_record(1+1e-8im).imaginary_passed==false
    @test a.inverse_record(1+.2im).phase≈-atan(.2)
end

@testset "No new background: algebraic saved-seed restoration" begin
    a=Q0_ENDPOINT_AUDIT;model=Main.Models.create_model(:PNJL)
    seed=[-1.,-1.1,-1.4,.1,.2,.5,.6,.2]
    bg=a.restore_seed(model,seed,140.,400.,1e-14)
    state=Main.Models.meanfield_state(seed[1:5])
    @test collect(bg.m)==Main.Models.calculate_mass_vec(model,state.phi)
    @test bg.mu==(u=.5,d=.6,s=.2)
    @test bg.Phi==.1 && bg.PhiBar==.2
    @test bg.T_MeV==140. && bg.muB_MeV==400.
    @test bg.coupling[:K_plus]==bg.coupling[:K_minus]
    @test bg.coupling[:pi_plus]==bg.coupling[:pi_minus]
    @test_throws ArgumentError a.restore_seed(model,seed[1:7],140.,400.,1e-14)
    @test_throws ArgumentError a.restore_seed(model,fill(NaN,8),140.,400.,1e-14)
    @test_throws ArgumentError a.restore_seed(model,seed,-140.,400.,1e-14)
    source=read(Q0_ENDPOINT_PATH,String)
    @test !occursin("Models.solve(",source)
    @test !occursin("W.background(",source)
    @test !occursin("channel_density(",source)
    s=a.restore_settings(JSON3.read("{\"mesh\":512,\"np\":128,\"nx\":64,\"ne\":128,\"nw\":4800,\"nr\":128,\"nq\":24,\"qmax\":8,\"thermal\":24,\"lower\":0.00001,\"upper\":56,\"phase_tol\":0.005}"))
    @test s.thermal==24. && s.lower==1e-5
    @test s.mesh==512 && s.nq==24
    @test s.upper==56. && s.nw==4800
end

@testset "Mapping, valid-domain flags and both historical q meshes are retained" begin
    a=Q0_ENDPOINT_AUDIT
    qs=a.q_probes(.4)
    @test qs==sort(unique(qs))
    @test first(qs)==0. && .2 in qs
    @test first(a.R.gauleg(0.,8.,8)[1]) in qs
    @test first(a.R.gauleg(0.,8.,24)[1]) in qs
    bg=(m=(u=1.,d=1.,s=2.),mu=(u=.5,d=.6,s=.1))
    cs=a.coordinate_rows(bg,:K_plus)
    r=only(filter(r->r.q_inv_fm==.2,cs))
    @test r.lambda_inv_fm≈sqrt(.4^2-.2^2)
    @test r.in_timelike_map && r.in_q0_landau && !r.in_q0_analytic_gap
    r=last(cs)
    @test !r.in_timelike_map && r.lambda_inv_fm===nothing
    @test all(r->!r.in_timelike_map,a.coordinate_rows(bg,:K_minus))
end

@testset "Finite IR window is not a zero-endpoint yield: constant-phase counterexample" begin
    a=Q0_ENDPOINT_AUDIT;T=.7;q=.2;phase=w->.3
    lo=a.window_record(phase,T,q,1e-4,.02;nodes=128)
    hi=a.window_record(phase,T,q,1e-5,.02;nodes=128)
    @test iszero(lo.derivative_inv_fm2) && iszero(hi.derivative_inv_fm2)
    @test lo.bulk_inv_fm2>0 && hi.bulk_inv_fm2>9lo.bulk_inv_fm2
    @test hi.lower_boundary_inv_fm2>9lo.lower_boundary_inv_fm2
    @test abs(lo.identity_error_inv_fm2)<1e-12
    @test abs(hi.identity_error_inv_fm2)<1e-12
    varying=a.window_record(w->.3+.2w,T,q,1e-4,.02;nodes=256)
    @test abs(varying.identity_error_inv_fm2)<1e-12
    @test_throws ArgumentError a.window_record(phase,T,q,0.,.02)
    # Exercise actual mixed null/Float64 CSV rows without any physical kernel.
    mktempdir() do dir
        rows=[merge((lambda_inv_fm=nothing,),a.inverse_record(1+0im)),
              merge((lambda_inv_fm=.2,),a.inverse_record(1+.01im))]
        path=joinpath(dir,"rows.csv");a.writecsv(path,rows)
        actual=collect(CSV.File(path))
        @test ismissing(actual[1].lambda_inv_fm)
        @test actual[2].lambda_inv_fm==.2
    end
end

@testset "Audit CLI and Actions never dispatch a full scan in audit mode" begin
    a=Q0_ENDPOINT_AUDIT
    @test a.parse_args(["--help"])===nothing
    opts=a.parse_args(["--input-root","x","--output","y"])
    @test isabspath(opts.input_root) && isabspath(opts.output)
    @test_throws ArgumentError a.parse_args(["--output","y"])
    @test_throws ArgumentError a.parse_args(["--input-root","x","--scan","y"])
    workflow=read(joinpath(a.ROOT,".github","workflows","relaxtime-charged-gbu-contour-scan.yml"),String)
    @test occursin("endpoint_audit:",workflow)
    @test occursin("inputs.replot_run_id == '' && !inputs.endpoint_audit",workflow)
    @test occursin("!inputs.endpoint_audit &&",workflow)
    @test occursin("run-id: 36971423991",workflow)
    @test occursin("audit_charged_gbu_q0_endpoint.jl",workflow)
end
