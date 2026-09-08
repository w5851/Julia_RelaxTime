using Test
const _FC_ROOT=normpath(joinpath(@__DIR__,"..","..",".."))
include(joinpath(_FC_ROOT,"scripts","analysis","relaxtime","causal_gbu_coordinate_oracles.jl"))
include(joinpath(_FC_ROOT,"scripts","analysis","relaxtime","causal_gbu_gap_completeness.jl"))
const FCO=CausalGBUCoordinateOracles
const FCG=CausalGBUGapCompleteness
include(joinpath(_FC_ROOT,"scripts","analysis","relaxtime","summarize_causal_gbu_full_closure.jl"))
const FCS=CausalGBUFullClosureSummary

@testset "Closure reducer rejects missing evidence and false passes" begin
    row=(max_coordinate_difference=0.0,max_radial_node_difference=0.0,contact_coordinate_difference=0.0,
        max_contact_identity_residual=0.0,max_contact_moment_residual=0.0,mesh_root_counts_stable=true,
        coordinate_passed=true,algebra_passed=true,representation_passed=true,all_real_gap_checks_passed=true,
        reciprocity_passed=true,cut_interpolation_passed=true,static_passed=true,passed=true)
    profiles=[(mesh=m,spectral_difference=0.0,coordinate_difference=0.0,radial_node_difference=0.0) for m in (16,32) for _ in 1:7]
    meshes=[(mesh=m,max_spectral_difference=0.0,max_direct_cut_interpolation_error=0.0,
        max_cut_reciprocity_error=0.0,signed_real_root_count=1,all_interpolant_gaps_passed=true,
        support_excess_width_inv_fm=0.0,support_missing_width_inv_fm=0.0,support_roundoff_allowance_inv_fm=1e-12,
        geometry_passed=true,tail_passed=true,tail_inverse_deviation_bound=0.5,static_passed=true) for m in (16,32)]
    gaps=[(mesh=m,root_count=1,passed=true) for m in (16,32)]
    cuts=[(mesh=m,interpolation_error=0.0,reciprocity_error=0.0) for m in (16,32)]
    @test FCS.check_case(row,profiles,meshes,gaps,cuts)
    @test_throws ErrorException FCS.check_case(row,profiles[2:end],meshes,gaps,cuts)
    @test_throws ErrorException FCS.check_case(row,profiles,meshes,gaps[2:end],cuts)
    @test_throws ErrorException FCS.check_case(row,profiles,meshes,gaps,cuts[2:end])
    @test_throws ErrorException FCS.check_case(merge(row,(passed=false,)),profiles,meshes,gaps,cuts)
    broken=[merge(m,(support_excess_width_inv_fm=1e-6,)) for m in meshes]
    @test_throws ErrorException FCS.check_case(row,profiles,broken,gaps,cuts)
    slow=[merge(p,(radial_node_difference=2e-6,)) for p in profiles]
    failed=merge(row,(max_radial_node_difference=2e-6,coordinate_passed=false,passed=false))
    @test !FCS.check_case(failed,slow,meshes,gaps,cuts)
    @test_throws ErrorException FCS.check_case(merge(failed,(passed=true,)),slow,meshes,gaps,cuts)
end

@testset "Same regulator in independent radial coordinates" begin
    a,b,u,v,T,phi,phibar,L,H=1.,1.4,0.2,-0.1,0.8,0.3,0.4,3.,6.
    zs=[0.4+0.6im,3.0+0.8im]
    for q in (0.,0.5,7.),regulator in (:two_line,:centered)
        radial=FCO.radial_loop(q,a,u,b,v,T,phi,phibar,L,H,zs,regulator;np=64,nr=64)
        if regulator===:two_line
            grid=FCO.R.build_spectral_bubble(q,a,u,b,v,T;Phi=phi,PhiBar=phibar,
                vacuum_cutoff_inv_fm=L,thermal_cutoff_inv_fm=H,momentum_nodes=96,angle_nodes=64)
            values=[FCO.R.spectral_bubble(grid,real(z)-u+v;eta_inv_fm=imag(z)).value for z in zs]
            contact=grid.contact_inv_fm2
        else
            grid=FCO.C.centered_atoms(q,a,u,b,v,T,phi,phibar,L,H;np=96,nx=64)
            values=[FCO.C.atom_value(grid,z).value for z in zs]
            contact=grid.contact
        end
        @test radial.values≈values atol=2e-8
        @test radial.contact≈contact atol=1e-9
        @test radial.values==radial.vacuum.values+radial.thermal.values
        for part in (radial.vacuum,radial.thermal)
            @test abs(part.moment0)<1e-12
            @test part.moment_identity_residual<1e-11
            @test part.contact_identity_residual<1e-11
        end
        if q>2L && regulator===:two_line
            @test all(iszero,radial.vacuum.values)
            @test radial.vacuum.contact==0
        end
    end
    @test isempty(FCO.radial_domain(6.,3.,:two_line))
    @test_throws ArgumentError FCO.radial_domain(0.,3.,:unknown)
    @test_throws ArgumentError FCO.radial_loop(1.,a,u,b,v,T,phi,phibar,L,H,[1.],:centered)
end

@testset "Geometric endpoints cannot leak into interpolant gaps" begin
    # Synthetic near-degenerate flavors exercise both signs and both regulators.
    # Tiny but physical thermal values INSIDE support must still be retained.
    bg=(m=(u=1.6519019722802555,d=1.649895213451757,s=2.5),
        mu=(u=0.3742643948705271,d=0.4209954886862639,s=0.1),T=0.8615142220054972,
        Phi=0.3,PhiBar=0.4,coupling=Dict(c=>0.2 for c in FCO.R.CHANNELS),vacuum=3.0)
    for ch in (:pi_plus,:pi_minus),q in (0.0,1.0,3.0,6.0),reg in (:two_line,:centered)
        a,b=FCO.R.charged_rpa_spec(ch).pair
        m1,m2,u,v=bg.m[a],bg.m[b],bg.mu[a],bg.mu[b]
        cuts=FCO.C.support_intervals(q,m1,m2,3.0,24.0,reg)
        if reg===:two_line
            g=FCO.R.build_spectral_bubble(q,m1,u,m2,v,bg.T;Phi=bg.Phi,PhiBar=bg.PhiBar,
                vacuum_cutoff_inv_fm=3.0,thermal_cutoff_inv_fm=24.0,momentum_nodes=8,angle_nodes=8)
            p=FCO.R.build_bubble_dispersion(g;segment_nodes=16)
        else
            p=FCO.C.centered_profile(bg,ch,q;mesh=16,ne=32).profile
        end
        @test FCG.support_geometry_check(FCG.signed_cells(p),cuts).passed
        outer=last(cuts)[2]
        for l in (outer,nextfloat(outer),-outer,prevfloat(-outer))
            @test imag(FCO.R.cauchy_transform(p,l))==0
        end
    end
    # The thermal tail is ~1e-24, still inside support and nonzero.
    cut=FCO.C.centered_cut(47.0,3.0,bg.m.u,bg.mu.u,bg.m.d,bg.mu.d,0.4,
        bg.Phi,bg.PhiBar,24.0;component=:thermal)
    @test 0<abs(cut)<1e-15
end

@testset "Signed-cell enclosures include zeros, gaps and infinite tails" begin
    p=FCG.R.PiecewiseSpectralFunction([1.,2.,3.],[0.,1.,0.])
    cells=FCG.signed_cells(p)
    @test length(cells)==2
    for (left,right) in ((-10.,0.9),(3.1,7.),(-1.,-0.2))
        bound=FCG.inverse_enclosure(cells,1.,left,right)
        for w in range(left,right;length=11)
            v=real(1-4FCG.R.cauchy_transform(p,w))
            @test bound.lower<=v<=bound.upper
        end
    end
    @test_throws ArgumentError FCG.inverse_enclosure(cells,1.,1.5,2.5)
    tail=FCG.tail_exclusion(cells,1.)
    @test tail.passed && tail.inverse_deviation_bound<=0.5
    for z in (tail.radius_inv_fm, -2tail.radius_inv_fm, tail.radius_inv_fm*im)
        @test abs(4FCG.R.cauchy_transform(p,z))<=tail.inverse_deviation_bound
    end
    audit=FCG.count_all_real_gaps(p,1.,0.)
    @test audit.passed
    @test audit.count==1
    @test all(r.endpoints_excluded for r in audit.rows)
    @test count(r->r.count_method=="bisection_and_argument_principle",audit.rows)==1
    # The root inside the OLD 1e-6 margin cannot disappear from a passed result.
    near=1.0-1e-7
    K=1/(4real(FCG.R.cauchy_transform(p,near)))
    narrow=FCG.count_all_real_gaps(p,K,0.)
    @test narrow.passed && narrow.count==1
    @test first(narrow.rows).margin_inv_fm<1e-7
    @test FCG.support_geometry_check(cells,[(1.,3.)]).passed
    @test !FCG.support_geometry_check(cells,[(1.01,3.)]).passed
    @test !FCG.support_geometry_check(cells,[(0.9,3.)]).passed
    # Tiny but nonzero spectra are not thresholded into fake analytic gaps.
    tiny=FCG.R.PiecewiseSpectralFunction([-1.,0.,1.],[0.,1e-200,0.])
    @test FCG.profile_support(FCG.signed_cells(tiny))==[(-1.,1.)]
    signed=FCG.R.PiecewiseSpectralFunction([-2.,-1.,1.,2.],[0.,-1.,1.,0.])
    @test length(FCG.signed_cells(signed))==4
    @test FCG.profile_support(FCG.signed_cells(signed))==[(-2.,2.)]
    rounded=FCG.R.PiecewiseSpectralFunction([-2.,-1.,1.,2.],[0.,-1e-30,1.,0.])
    triangles=FCG.signed_cells(rounded)
    @test length(triangles)==4
    for w in (-4.,4.)
        @test sum(FCG.cell_real(c,w) for c in triangles)/pi≈real(FCG.R.cauchy_transform(rounded,w)) atol=1e-14
    end
    # Even root: winding finds multiplicity two; a simple-root scan MUST fail.
    even=FCG.R.PiecewiseSpectralFunction([-3.,-2.,-1.,1.,2.,3.],[0.,-1.,0.,0.,1.,0.])
    coupling=1/(4real(FCG.R.cauchy_transform(even,0.)))
    audit_even=FCG.count_all_real_gaps(even,coupling,0.)
    @test !audit_even.passed
    @test any(r->r.contour_count==2,audit_even.rows)
end
