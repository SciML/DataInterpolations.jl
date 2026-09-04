using ChainRulesCore
using LogExpFunctions
using StatsFuns
using Test

@testset "Dependency extensions" begin
    @test Base.get_extension(
        LogExpFunctions, :LogExpFunctionsChainRulesCoreExt
    ) isa Module
    @test Base.get_extension(StatsFuns, :StatsFunsChainRulesCoreExt) isa Module
end
