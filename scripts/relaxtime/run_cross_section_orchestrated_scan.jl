#!/usr/bin/env julia

module CrossSectionOrchestratedScan

export run_cross_section_orchestrated

function _energy_grid(xs::Dict{String,Any})
    energy = get(xs, "energy", Dict{String,Any}())
    mode = String(get(energy, "mode", "linspace"))
    if mode == "list"
        vals = Float64.(get(energy, "sqrt_s_list_MeV", Any[]))
        return unique(sort(vals))
    end
    smin = Float64(get(energy, "sqrt_s_min_MeV", 300.0))
    smax = Float64(get(energy, "sqrt_s_max_MeV", 2000.0))
    n = Int(get(energy, "sqrt_s_num", 64))
    return collect(range(smin, smax; length=n))
end

function run_cross_section_orchestrated(effective::Dict{String,Any}, outdir::String; run_id::String)
    scan = get(effective, "scan", Dict{String,Any}())
    xs = get(scan, "cross_section", Dict{String,Any}())

    T_list = Float64.(get(xs, "T_list_MeV", Any[]))
    xi_list = Float64.(get(xs, "xi_list", Any[]))
    processes = String.(get(xs, "processes", Any[]))
    muB = Float64(get(xs, "muB_MeV", 0.0))
    sqrt_s_values = _energy_grid(xs)

    out_csv = joinpath(outdir, "cross_section_orchestrated.csv")
    open(out_csv, "w") do io
        println(io, "run_id,T_MeV,muB_MeV,xi,process,sqrt_s_MeV")
        for T in T_list
            for xi in xi_list
                for p in processes
                    for s in sqrt_s_values
                        println(io, string(run_id, ',', T, ',', muB, ',', xi, ',', p, ',', s))
                    end
                end
            end
        end
    end

    return out_csv
end

end # module
