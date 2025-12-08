using Test
using BenchmarkTools
using Statistics: mean

include(joinpath(@__DIR__, "..", "..", "..", "src", "pnjl", "PNJL.jl"))
using .PNJL: ThermoDerivatives

const BENCH_CONFIG = (
    T_mev = 150.0,
    mu_mev = 10.0,
    xi = 0.0,
    p_num = 16,
    t_num = 8,
)

@testset "bulk_derivative_coeffs performance (average runtime)" begin
    ThermoDerivatives.bulk_derivative_coeffs(
        BENCH_CONFIG.T_mev,
        BENCH_CONFIG.mu_mev;
        xi = BENCH_CONFIG.xi,
        p_num = BENCH_CONFIG.p_num,
        t_num = BENCH_CONFIG.t_num,
    )

    bench = @benchmark ThermoDerivatives.bulk_derivative_coeffs(
        $(BENCH_CONFIG.T_mev),
        $(BENCH_CONFIG.mu_mev);
        xi = $(BENCH_CONFIG.xi),
        p_num = $(BENCH_CONFIG.p_num),
        t_num = $(BENCH_CONFIG.t_num),
    ) samples = 3 evals = 1

    avg_ms = mean(bench.times) / 1.0e6
    @info "bulk_derivative_coeffs average runtime (p=$(BENCH_CONFIG.p_num), t=$(BENCH_CONFIG.t_num))" average_ms = avg_ms samples = length(bench.times)
    @test avg_ms > 0.0
end
"""
测试结果（实际结果可能因硬件和环境不同而异）：
┌ Info: bulk_derivative_coeffs average runtime (p=16, t=8)
│   average_ms = 19.8912
│   samples = 3
"""