include(joinpath(@__DIR__, "..", "..", "src", "relaxtime", "TransportCoefficients.jl"))
using .TransportCoefficients
using Statistics

qp = (m=(u=0.3, d=0.3, s=0.5), μ=(u=0.2, d=0.2, s=0.2))
tp = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)
tau = (u=1.0, d=1.0, s=1.0, ubar=1.0, dbar=1.0, sbar=1.0)
cfg = TransportIntegrationConfig(p_nodes=32, p_max=10.0)

shear_viscosity(qp, tp; tau=tau, config=cfg)
electric_conductivity(qp, tp; tau=tau, config=cfg)

n = 20
t_eta = Float64[]
t_sigma = Float64[]
for _ in 1:n
    push!(t_eta, @elapsed shear_viscosity(qp, tp; tau=tau, config=cfg))
    push!(t_sigma, @elapsed electric_conductivity(qp, tp; tau=tau, config=cfg))
end

println("eta_median_s=", median(t_eta))
println("sigma_median_s=", median(t_sigma))
println("eta_min_s=", minimum(t_eta))
println("sigma_min_s=", minimum(t_sigma))
