using BenchmarkTools
using DataInterpolations
using DataInterpolations: makima_interpolate!
using StableRNGs

seed = 1234
rng = StableRNG(seed)
n, n_tt, n_batch = 32, 512, 256

function akima_interpolation_vector(r, u, t, tt)
    itp = AkimaInterpolation(u, t; extrapolation=ExtrapolationType.Constant)
    itp(r, tt)
end

function akima_interpolation_matrix(r, u, t, tt)
    @inbounds @views for j in axes(u, 2)
        itp = AkimaInterpolation(u[:, j], t[:, j]; extrapolation=ExtrapolationType.Constant)
        itp(r[:, j], tt[:, j])
    end
end

let 
    println("="^8, "akima vs makima - vector", "="^8)

    u = rand(rng, n);
    t = sort!(rand(rng, n) * 10.0);
    tt = sort!(rand(rng, n_tt) * 10.0);
    r = zeros(n_tt);

    println("DataInterpolation.jl")
    @time akima_interpolation_vector(r, u, t, tt)
    display(@benchmark akima_interpolation_vector($r, $u, $t, $tt))

    println("makima vector:")
    @time makima_interpolate!(r, u, t, tt)
    display(@benchmark makima_interpolate!($r, $u, $t, $tt))
end

let 
    println("="^8, "akima vs makima - matrix", "="^8)

    u = rand(rng, n, n_batch);
    t = sort!(rand(rng, n, n_batch) * 10.0, dims=1);
    tt = sort!(rand(rng, n_tt, n_batch) * 10.0, dims=1);
    r = zeros(n_tt, n_batch);

    println("DataInterpolation.jl matrix version via loop")
    @time akima_interpolation_matrix(r, u, t, tt)
    display(@benchmark akima_interpolation_matrix($r, $u, $t, $tt))

    println("makima matrix:")
    @time makima_interpolate!(r, u, t, tt)
    display(@benchmark makima_interpolate!($r, $u, $t, $tt))
end


