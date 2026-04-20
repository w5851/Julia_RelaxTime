"""
    AbstractSusceptibilityProvider

守恒荷涨落提供器抽象层。
"""
module AbstractSusceptibilityProvider

export AbstractChiProvider, chi_provider_name

abstract type AbstractChiProvider end

chi_provider_name(::AbstractChiProvider) = :unknown

end # module AbstractSusceptibilityProvider
