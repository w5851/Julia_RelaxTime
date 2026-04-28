using Test

const _REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const _SCRIPT_PATH = joinpath(_REPO_ROOT, "scripts", "relaxtime", "generate_xi_smoothness_params.jl")

function _write_anchor_csv(path::String, rows::Vector{String}; header::String)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, header)
        for row in rows
            println(io, row)
        end
    end
end

function _read_csv_rows(path::String)
    header = String[]
    rows = Vector{Dict{String, String}}()
    open(path, "r") do io
        header_seen = false
        for line in eachline(io)
            s = strip(line)
            isempty(s) && continue
            startswith(s, "#") && continue
            if !header_seen
                header = [strip(x) for x in split(s, ',')]
                header_seen = true
                continue
            end
            vals = split(line, ',')
            d = Dict{String, String}()
            for (i, k) in enumerate(header)
                d[k] = i <= length(vals) ? strip(vals[i]) : ""
            end
            push!(rows, d)
        end
    end
    return header, rows
end

@testset "generate xi smoothness params cli" begin
    @test isfile(_SCRIPT_PATH)

    tmp = mktempdir()
    boundary_csv = joinpath(tmp, "boundary.csv")
    crossover_csv = joinpath(tmp, "crossover.csv")
    _write_anchor_csv(
        boundary_csv,
        ["0.0,100.0,300.0,0.0,0.0"];
        header="xi,T_MeV,mu_transition_MeV,rho_hadron,rho_quark",
    )
    _write_anchor_csv(
        crossover_csv,
        ["0.0,120.0,180.0,165.0,0.0,0.0"];
        header="xi,mu_MeV,T_crossover_chiral_MeV,T_crossover_deconf_MeV,rho_chiral,rho_deconf",
    )

    @testset "random anchor fields are blank; near anchors numeric" begin
        output_csv = joinpath(tmp, "params.csv")
        run(`julia --project=. $_SCRIPT_PATH --seed 11 --output $output_csv --random-count 2 --near-count 1 --tmin 50 --tmax 270 --muqmin 0 --muqmax 360 --boundary-csv $boundary_csv --crossover-csv $crossover_csv`)
        @test isfile(output_csv)

        _, rows = _read_csv_rows(output_csv)
        @test length(rows) == 3

        random_rows = [r for r in rows if r["source"] == "random_uniform"]
        near_rows = [r for r in rows if r["source"] == "near_phase_line"]
        @test length(random_rows) == 2
        @test length(near_rows) == 1

        for r in random_rows
            @test r["anchor_type"] == ""
            @test r["anchor_T_MeV"] == ""
            @test r["anchor_muq_MeV"] == ""
            @test r["delta_T"] == ""
            @test r["delta_muq"] == ""
        end

        near = near_rows[1]
        @test near["anchor_type"] != ""
        @test parse(Float64, near["anchor_T_MeV"]) isa Float64
        @test parse(Float64, near["anchor_muq_MeV"]) isa Float64
        @test parse(Float64, near["delta_T"]) isa Float64
        @test parse(Float64, near["delta_muq"]) isa Float64
    end

    @testset "default output filename uses dynamic total" begin
        seed = 31415
        random_count = 3
        near_count = 4
        total = random_count + near_count
        expected = joinpath(
            _REPO_ROOT,
            "data",
            "outputs",
            "results",
            "relaxtime",
            "plan_c",
            "sampling",
            "params_$(total)_seed$(seed).csv",
        )

        isfile(expected) && rm(expected; force=true)
        run(`julia --project=. $_SCRIPT_PATH --seed $seed --random-count $random_count --near-count $near_count --boundary-csv $boundary_csv --crossover-csv $crossover_csv`)
        @test isfile(expected)
        rm(expected; force=true)
    end
end
