using DataInterpolations
using DataInterpolations: makima_init!, makima_init_slow!, makima_interpolate!
using StableRNGs
using Test

function coefficients(u::AbstractVector{T}, t::AbstractVector{T}) where {T}
    n = length(u)
    b = Vector{T}(undef, n)
    c = Vector{T}(undef, n - 1)
    d = Vector{T}(undef, n - 1)
    makima_init!(u, t, b, c, d)
    return b, c, d
end

function coefficients_slow(u::AbstractVector{T}, t::AbstractVector{T}) where {T}
    n = length(u)
    b = Vector{T}(undef, n)
    c = Vector{T}(undef, n - 1)
    d = Vector{T}(undef, n - 1)
    makima_init_slow!(u, t, b, c, d)
    return b, c, d
end

function strict_random_grid(rng, ::Type{T}, n::Int) where {T}
    gaps = T(0.1) .+ rand(rng, T, n - 1)
    t = Vector{T}(undef, n)
    t[1] = zero(T)
    @inbounds for i in 2:n
        t[i] = t[i - 1] + gaps[i - 1]
    end
    return t
end

@testset "Makima module" begin
    @testset "sin approximation" begin
        n = 48
        ntt = 1000
        t = collect(range(0.0, 2pi; length=n))
        u = sin.(t)
        tt = collect(range(first(t), last(t); length=ntt))
        r = similar(tt)

        makima_interpolate!(r, u, t, tt)

        @test maximum(abs.(r .- sin.(tt))) < 2e-3
        @test isapprox(r[1], u[1]; atol=1e-14, rtol=1e-12)
        @test isapprox(r[end], u[end]; atol=1e-14, rtol=1e-12)
    end

    @testset "fast init matches slow init" begin
        rng = StableRNG(1234)

        for T in (Float32, Float64), n in (3, 4, 8, 17, 64)
            for _ in 1:20
                t = strict_random_grid(rng, T, n)
                u = randn(rng, T, n)

                b, c, d = coefficients(u, t)
                b_ref, c_ref, d_ref = coefficients_slow(u, t)

                @test b ≈ b_ref
                @test c ≈ c_ref
                @test d ≈ d_ref
            end
        end
    end
end
