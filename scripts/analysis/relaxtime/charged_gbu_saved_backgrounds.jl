"""Restore immutable contour backgrounds algebraically, without an equilibrium solve."""
module ChargedGBUSavedBackgrounds

using JSON3, SHA
using Main.Models

const W = Main.ChargedGBUResearchWorkflow
const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
hashfile(path) = bytes2hex(sha256(read(path)))

function restore_seed(model, seed, T_MeV, muB_MeV, residual)
    length(seed) == 8 && all(isfinite, seed) || throw(ArgumentError("finite eight-component saved seed required"))
    all(isfinite, (T_MeV, muB_MeV, residual)) && T_MeV > 0 && muB_MeV >= 0 && residual >= 0 ||
        throw(ArgumentError("invalid saved background coordinates or residual"))
    state = Models.meanfield_state(Float64.(seed[1:5]))
    masses = Models.calculate_mass_vec(model, state.phi)
    interaction = Main.RelaxTime.MesonInteractionKernel.build_full_kmt_interaction(
        state.phi; G=model.params.G_fm2, K=model.params.K_fm5)
    coupling = Dict(ch => Main.RelaxTime.ChargedRPAKernel.charged_rpa_coupling(
        W.R.charged_rpa_spec(ch), interaction) for ch in W.R.CHANNELS)
    tri(x) = (u=Float64(x[1]), d=Float64(x[2]), s=Float64(x[3]))
    return (m=tri(masses), mu=tri(seed[6:8]), T=T_MeV/Main.Constants_PNJL.ħc_MeV_fm,
        Phi=Float64(state.Phi), PhiBar=Float64(state.PhiBar), coupling=coupling,
        vacuum=Float64(model.params.Λ_inv_fm), T_MeV=Float64(T_MeV),
        muB_MeV=Float64(muB_MeV), residual=Float64(residual))
end

function read_snapshot(input_root, T_values, muB_values; config_path=W.DEFAULT_CONFIG)
    root = abspath(input_root)
    isdir(root) || throw(ArgumentError("background input root does not exist"))
    manifests = sort([joinpath(d,f) for (d,_,fs) in walkdir(root) for f in fs if f == "manifest.json"])
    isempty(manifests) && throw(ArgumentError("no saved contour manifests"))
    points = Dict{Tuple{Int,Int},Any}()
    hashes = Dict{String,String}()
    heads = String[]
    for path in manifests
        m = JSON3.read(read(path, String))
        m.schema == "charged_gbu_contour_scan_v2" || throw(ArgumentError("unsupported saved background schema"))
        Float64.(m.T_grid) == T_values && Float64.(m.muB_grid) == muB_values ||
            throw(ArgumentError("saved background grid mismatch"))
        # Only density/postprocessing code may differ from the saved scan.
        for rel in ("src/models/Models.jl",
                "src/models/solver/runtime/ConstraintSolverFixedMuBConservedCharges.jl",
                "src/models/workflow_apps/ChargedGBUResearchWorkflow.jl",
                replace(relpath(config_path, ROOT), '\\'=>'/'))
            haskey(m.source_hashes, rel) && hashfile(joinpath(ROOT, rel)) == m.source_hashes[rel] ||
                throw(ArgumentError("saved background source mismatch: $(rel)"))
        end
        hashes[replace(relpath(path, root), '\\'=>'/')] = hashfile(path)
        push!(heads, String(m.git_head))
        files = sort(readdir(joinpath(dirname(path), "points"); join=true))
        for file in filter(p -> endswith(p, ".json"), files)
            d = JSON3.read(read(file, String))
            d.schema == "charged_gbu_contour_point_v2" && d.scan_identity == m.scan_identity ||
                throw(ArgumentError("saved background point identity mismatch"))
            i, j = Int(d.row_index), Int(d.col_index)
            1 <= i <= length(T_values) && 1 <= j <= length(muB_values) &&
                d.T_MeV == T_values[i] && d.muB_MeV == muB_values[j] ||
                throw(ArgumentError("saved background coordinates mismatch"))
            haskey(points, (i,j)) && throw(ArgumentError("duplicate saved background point"))
            points[(i,j)] = d
            hashes[replace(relpath(file, root), '\\'=>'/')] = hashfile(file)
        end
    end
    length(points) == length(T_values)*length(muB_values) || throw(ArgumentError("incomplete saved background grid"))
    fingerprint = bytes2hex(sha256(JSON3.write(sort!(collect(hashes); by=first))))
    return (points=points, provenance=(mode="saved_seed_algebraic_restoration", solver_called=false,
        source_git_heads=sort!(unique(heads)), input_hashes=hashes, fingerprint=fingerprint))
end

end
