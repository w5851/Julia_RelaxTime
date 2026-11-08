using Test
using ForwardDiff

if !isdefined(Main, :RelaxTime)
    include(joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "RelaxTime.jl"))
end

@testset "RS response conventions: independent analytic limits" begin
    # Exact fixed-background identities, no gap solve or production scan.
    nodes, weights = Main.RelaxTime.GaussLegendre.gauleg(-1.0, 1.0, 96)
    angular(f) = sum(w * f(c) for (c, w) in zip(nodes, weights)) / 2
    # Central difference truncation is O(h^2); h=2e-5 keeps it below 1e-9
    # for the analytic angular functions used here, without relaxing a gate.
    h = 2e-5
    centered(f) = (f(h) - f(-h)) / (2h)
    density(xi) = angular(c -> (1 + xi*c^2)^(-3/2))
    pressure(xi) = angular(c -> (1 + xi*c^2)^(-2))
    @test isapprox(density(0.3), inv(sqrt(1.3)); rtol=2e-13)
    @test isapprox(density(-0.3), inv(sqrt(0.7)); rtol=2e-13)
    @test isapprox(centered(density), -1/2; atol=4e-9)
    @test isapprox(centered(pressure), -2/3; atol=6e-9)
    @test isapprox(centered(x -> pressure(x)/density(x)), -1/6; atol=3e-9)

    # Massless FD with constant tau: radial integral factors cancel.
    # sigma_occ/sigma0 = <a^-3>, eta_occ/eta0 = <a^-5>.
    eta_occ(xi) = angular(c -> (1 + xi*c^2)^(-5/2))
    sigma_perp(xi) = 3/2 * angular(c -> (1-c^2)/(1+xi*c^2)^2)
    sigma_parallel(xi) = 3 * (1+xi) * angular(c -> c^2/(1+xi*c^2)^2)
    sigma_trace(xi) = (2sigma_perp(xi) + sigma_parallel(xi))/3
    @test isapprox(centered(eta_occ), -5/6; atol=9e-9)
    @test isapprox(centered(sigma_perp), -2/5; atol=4e-9)
    @test isapprox(centered(sigma_parallel), -1/5; atol=4e-9)
    @test isapprox(centered(sigma_trace), -1/3; atol=4e-9)
    @test isapprox(centered(x -> sigma_trace(x)-density(x)), 1/6; atol=4e-9)

    # An unweighted scalar shear contraction averages the five traceless
    # symmetric components. Individual components differ at O(xi).
    eta_xy(xi) = 15/8 * angular(c -> (1-c^2)^2/(1+xi*c^2)^(5/2))
    eta_xz(xi) = 15/2 * angular(c -> c^2*(1-c^2)/(1+xi*c^2)^(5/2))
    eta_axial(xi) = 45/4 * angular(c -> (c^2-1/3)^2/(1+xi*c^2)^(5/2))
    for xi in (-0.3, 0.0, 0.3)
        @test isapprox((2eta_xy(xi)+2eta_xz(xi)+eta_axial(xi))/5, eta_occ(xi); rtol=1e-12)
    end
    @test isapprox(centered(eta_xy), -5/14; atol=5e-9)
    @test isapprox(centered(eta_xz), -15/14; atol=1e-8)
    @test isapprox(centered(eta_axial), -55/42; atol=2e-8)

    # Massive thermal condensates: same physical integration domain and
    # independently expanded 1/E_xi.  The p=30 tail is exponentially small.
    ps, ws = Main.RelaxTime.GaussLegendre.gauleg(0.0, 30.0, 160)
    mass=0.8; temp_c=1.0; chem_c=0.3; nc=3
    qdist_c = Main.PNJLQuarkDistributions.quark_distribution
    adist_c = Main.PNJLQuarkDistributions.antiquark_distribution
    for phi in (0.4, 1.0)
        thermal(E) = qdist_c(E,chem_c,temp_c,phi,phi)+adist_c(E,chem_c,temp_c,phi,phi)
        function condensate_thermal(xi, variational)
            total=0.0
            for (p,w) in zip(ps,ws), (c,wc) in zip(nodes,weights)
                E=sqrt(p^2+mass^2); Ed=sqrt(E^2+xi*p^2*c^2)
                total += w*wc*p^2*thermal(Ed)/(variational ? Ed : E)
            end
            return nc*mass/(2pi^2)*total
        end
        slope = -nc*mass/(6pi^2)*sum(w*p^4/(p^2+mass^2)^(3/2)*thermal(sqrt(p^2+mass^2)) for (p,w) in zip(ps,ws))
        delta(xi)=condensate_thermal(xi,true)-condensate_thermal(xi,false)
        @test isapprox(delta(0.0),0.0;atol=1e-14)
        @test isapprox(centered(delta),slope;rtol=1e-7)
        a=Main.OneLoopIntegralsCorrection.A_aniso(mass,chem_c,temp_c,phi,phi,0.3,ps,ws,nodes,weights)
        vacuum=-nc*mass/pi^2*Main.OneLoopIntegrals.const_integral_term_A(mass)
        @test isapprox(nc*mass*a/(4pi^2)-vacuum,condensate_thermal(0.3,false);rtol=1e-12)
    end

    # Source gradient in the elementary FD case.
    p = [0.7, -0.2, 0.8]; m=0.35; T=0.8; mu=0.12; xi=0.2
    energy(v) = sqrt(sum(abs2, v)+m^2+xi*v[3]^2)
    occupation(v) = inv(1+exp((energy(v)-mu)/T))
    f=occupation(p)
    expected = -f*(1-f)/T/energy(p) .* [p[1],p[2],(1+xi)*p[3]]
    @test ForwardDiff.gradient(occupation,p) ≈ expected

    # Actual PNJL occupation: the color-average does not obey the FD identity.
    qdist = Main.PNJLQuarkDistributions.quark_distribution
    E=1.2; temp=0.7; chem=0.1
    for phi in (0.0, 0.4, 1.0)
        f0=qdist(E,chem,temp,phi,phi)
        W=-temp*ForwardDiff.derivative(e -> qdist(e,chem,temp,phi,phi),E)
        z=exp(-(E-chem)/temp)
        D=1+3phi*z+3phi*z^2+z^3
        algebra=(phi*z+4phi*z^2+3z^3)/D-3*f0^2
        @test isapprox(W,algebra;rtol=2e-12)
        if phi==0.0
            @test isapprox(W,3*f0*(1-f0);rtol=2e-12)
        elseif phi==1.0
            @test isapprox(W,f0*(1-f0);rtol=2e-12)
        else
            @test !isapprox(W,f0*(1-f0);rtol=1e-3)
        end
    end
end
