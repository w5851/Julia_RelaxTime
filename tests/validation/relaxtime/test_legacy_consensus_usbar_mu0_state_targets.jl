using Test

if !isdefined(Main, :_compute_relaxtime_literature_transport_point)
    include(joinpath(@__DIR__, "literature_validation_helpers.jl"))
end

const RELAXTIME_USBAR_LEGACY_STATE_DATA_PATH = validation_targets_path(
    "relaxtime",
    "legacy",
    "state",
    "relaxtime_usbar_mu0_legacy_state_targets_v1.csv",
)

function _load_usbar_legacy_state_targets(path::String)
    isfile(path) || error("legacy state validation target CSV not found: $path")
    lines = readlines(path)
    isempty(lines) && error("legacy state validation target CSV is empty: $path")

    rows = NamedTuple[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue
        cols = split(s, ',')
        length(cols) == 9 || error("invalid legacy state validation target row: $line")
        push!(rows, (
            record_id=strip(cols[1]),
            T_MeV=parse(Float64, strip(cols[2])),
            muB_MeV=parse(Float64, strip(cols[3])),
            field=strip(cols[4]),
            fortran_value=parse(Float64, strip(cols[5])),
            cpp_value=parse(Float64, strip(cols[6])),
            consensus_value=parse(Float64, strip(cols[7])),
            rtol=parse(Float64, strip(cols[8])),
            source=strip(cols[9]),
        ))
    end
    return rows
end

const RELAXTIME_USBAR_LEGACY_STATE_CACHE = Dict{Tuple{Float64, Float64}, NamedTuple}()

function _compute_usbar_legacy_state_point(T_MeV::Float64, muB_MeV::Float64)
    key = (T_MeV, muB_MeV)
    haskey(RELAXTIME_USBAR_LEGACY_STATE_CACHE, key) && return RELAXTIME_USBAR_LEGACY_STATE_CACHE[key]

    T_fm = T_MeV / RELAXTIME_LITERATURE_VALIDATION_HBARC_MEV_FM
    muq_fm = (muB_MeV / 3.0) / RELAXTIME_LITERATURE_VALIDATION_HBARC_MEV_FM
    equilibrium = _solve_relaxtime_literature_validation_equilibrium(T_fm, muq_fm, 0.0)
    Bool(equilibrium.converged) || error("equilibrium solve failed at T=$(T_MeV) MeV, muB=$(muB_MeV) MeV")

    RELAXTIME_USBAR_LEGACY_STATE_CACHE[key] = (
        phi_u=Float64(equilibrium.x_state[1]),
        phi_s=Float64(equilibrium.x_state[3]),
        Phi=Float64(equilibrium.x_state[4]),
        Phibar=Float64(equilibrium.x_state[5]),
        m_u_MeV=Float64(equilibrium.masses[1]) * RELAXTIME_LITERATURE_VALIDATION_HBARC_MEV_FM,
        m_s_MeV=Float64(equilibrium.masses[3]) * RELAXTIME_LITERATURE_VALIDATION_HBARC_MEV_FM,
    )
    return RELAXTIME_USBAR_LEGACY_STATE_CACHE[key]
end

@testset "RelaxTime legacy consensus usbar states" begin
    targets = _load_usbar_legacy_state_targets(RELAXTIME_USBAR_LEGACY_STATE_DATA_PATH)
    @test !isempty(targets)

    for row in targets
        actual_state = _compute_usbar_legacy_state_point(row.T_MeV, row.muB_MeV)
        actual = getproperty(actual_state, Symbol(row.field))
        @test isfinite(actual)
        @test isapprox(actual, row.consensus_value; rtol=row.rtol, atol=0.0)
    end
end