"""
    ConstraintModes

PNJL 侧约束模式桥接层：复用 Models 域定义，避免双份类型实现漂移。
"""
module ConstraintModes

export ConstraintMode, FixedMu, FixedRho, FixedAsymmetricRho, FixedEntropy, FixedSigma
export state_dim, param_dim, constraint_description

const _MODELS_PATH = normpath(joinpath(@__DIR__, "..", "..", "Models.jl"))
if !isdefined(Main, :Models)
    Base.include(Main, _MODELS_PATH)
end

import Main.Models: ConstraintMode, FixedMu, FixedRho, FixedAsymmetricRho, FixedEntropy, FixedSigma
import Main.Models: state_dim, param_dim, constraint_description

end # module ConstraintModes

