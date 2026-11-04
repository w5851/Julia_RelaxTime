using Test, JSON3, SHA
if !isdefined(Main,:NativeTransportRegression)
    include(joinpath(@__DIR__,"..","..","regression","relaxtime","native_transport_fixture.jl"))
end

@testset "Read-only transport replay binds numerical content across source roots" begin
    old_root="D:/frozen/project"
    current_root="/tmp/current-project"
    relative="config/models/pnjl/default.toml"
    original=(schema_version=1,profile_id="thermal_weighted_v1",profile_sha256=repeat("a",64),
        propagator_xi_policy="match_thermo",compute_bulk=true,
        model_parameters=(path=old_root*"/"*relative,G_fm2=1.835,Λ_inv_fm=3.052),
        sources=[(path=relative,sha256=repeat("b",64)),(path="src/solver.jl",sha256=repeat("c",64))])
    saved_contract(identity)=JSON3.read(JSON3.write((identity_sha256=bytes2hex(sha256(JSON3.write(identity))),identity=identity)))
    saved=saved_contract(original)
    current=merge(original,(model_parameters=merge(original.model_parameters,(path=current_root*"/"*relative,)),))
    matches(s,i)=Main.NativeTransportRegression.replay_identity_matches(s,i;recorded_root=old_root,current_root=current_root)
    @test Main.NativeTransportRegression.replay_identity_matches(saved,original;recorded_root=old_root,current_root=old_root)
    @test matches(saved,current)
    for changed in (merge(current,(model_parameters=merge(current.model_parameters,(G_fm2=2.0,)),)),
                    merge(current,(profile_sha256=repeat("d",64),)),
                    merge(current,(compute_bulk=false,)),
                    merge(current,(propagator_xi_policy="isotropic",)),
                    merge(current,(sources=[current.sources[1],merge(current.sources[2],(sha256=repeat("e",64),))],)))
        @test !matches(saved,changed)
    end
    different_file=merge(original,(model_parameters=merge(original.model_parameters,(path=old_root*"/config/models/pnjl/other.toml",)),))
    @test !matches(saved_contract(different_file),current)
    outside=merge(original,(model_parameters=merge(original.model_parameters,(path="D:/another-project/"*relative,)),))
    @test !matches(saved_contract(outside),current)
    traversal=merge(original,(model_parameters=merge(original.model_parameters,(path=old_root*"/config/models/pnjl/../pnjl/default.toml",)),))
    @test !matches(saved_contract(traversal),current)
    missing_source=merge(original,(sources=original.sources[2:end],))
    @test !matches(saved_contract(missing_source),merge(current,(sources=current.sources[2:end],)))
    tampered=JSON3.read(JSON3.write((identity_sha256=saved.identity_sha256,
        identity=merge(original,(model_parameters=merge(original.model_parameters,(G_fm2=9.0,)),)))))
    @test !matches(tampered,current)
    @test !matches((identity_sha256=repeat("0",64),identity=saved.identity),current)
end
