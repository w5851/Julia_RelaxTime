using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "Magnetic node caches under concurrent first use" for M in
    (Models.MagneticIntegrals, Models.MagneticThermodynamics.MagneticIntegrals)
    # Unusual orders exercise initially empty entries, including shared orders.
    orders = repeat([37, 41, 43, 47], 16)
    for node_function in (M.pz_nodes, M.zeta_nodes)
        tasks = [Threads.@spawn node_function(n) for n in orders]
        results = fetch.(tasks)
        for (n, (nodes, weights)) in zip(orders, results)
            cached_nodes, cached_weights = node_function(n)
            @test nodes === cached_nodes
            @test weights === cached_weights
            @test issorted(nodes) && all(x -> 0.0 < x < 1.0, nodes)
            @test all(>(0.0), weights)
            @test sum(weights) ≈ 1.0 atol=2e-14
            @test sum(weights .* nodes.^4) ≈ 0.2 atol=2e-14
        end
    end

    # Consumers integrate after releasing the node-cache lock.
    kernel(mass) = (M.omega_magnetic_mfir(mass, 2/3, 0.15; zeta_num=43),
        M.omegat_flavor_landau(mass, 0.2, 0.8, 0.3, 0.3, 2/3, 0.15;
            n_max=3, p_num=37))
    masses = collect(range(0.8, 1.4; length=8))
    serial = kernel.(masses)
    parallel = fetch.([Threads.@spawn kernel(mass) for mass in masses])
    @test parallel == serial
    @test all(value -> all(isfinite, value), parallel)
    @test_throws ArgumentError M.pz_nodes(3)
    @test_throws ArgumentError M.zeta_nodes(7)
end
Models.pnjl_module()
const PNJL = Models.pnjl_module()
const alpha_n = getproperty(PNJL, :alpha_n)
const energy_landau = getproperty(PNJL, :energy_landau)
const resolve_nmax_from_cutoff = getproperty(PNJL, :resolve_nmax_from_cutoff)
const smooth_cutoff = getproperty(PNJL, :smooth_cutoff)

@testset "magnetic energy levels" begin
    @test alpha_n(0) == 1.0
    @test alpha_n(1) == 2.0
    @test alpha_n(5) == 2.0

    m = 1.2
    pz = 0.8
    q = 2 / 3
    eB = 0.15

    E0 = energy_landau(m, pz, 0, q, eB)
    E1 = energy_landau(m, pz, 1, q, eB)
    @test isfinite(E0)
    @test isfinite(E1)
    @test E1 > E0

    nmax = resolve_nmax_from_cutoff(m, 1.0, q, eB)
    @test nmax >= 0

    val = smooth_cutoff(1.0; Λ=3.0, N=10)
    @test 0.0 < val <= 1.0
end
