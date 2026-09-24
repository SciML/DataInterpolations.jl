using DataInterpolations, Test

# These names are documented for downstream use but not exported, so the `public`
# declaration is the only thing making them reachable as public API. Docstring coverage for
# them is checked by Aqua in the QA group.
@testset "Public API" begin
    @testset "$name" for name in
        (
            :AbstractInterpolation, :AbstractIntegralInverseInterpolation,
            :derivative, :integral, :invert_integral, :_interpolate, :_derivative,
            :_integral,
        )
        @test isdefined(DataInterpolations, name)
        if VERSION >= v"1.11"
            @test Base.ispublic(DataInterpolations, name)
        end
    end
end
