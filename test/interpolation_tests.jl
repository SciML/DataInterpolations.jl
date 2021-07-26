using DataInterpolations, Test
using StableRNGs

# Linear Interpolation
u = 2.0collect(1:10)
t = 1.0collect(1:10)
A = LinearInterpolation(u,t)

@test A(1) == u[1]
@test A(5) == u[5]
@test A(5.5) == 11.0

u = vcat(2.0collect(1:10)', 3.0collect(1:10)')
A = LinearInterpolation(u,t)

@test A(1) == u[:,1]
@test A(5) == u[:,5]
@test A(5.5) == [11.0,16.5]

# Quadratic Interpolation
u = [1.0, 4.0, 9.0, 16.0]
t = [1.0, 2.0, 3.0, 4.0]
A = QuadraticInterpolation(u,t)

@test A(2.0) == 4.0
@test A(1.5) == 2.25
@test A(3.5) == 12.25
@test A(2.5) == 6.25

u = [1.0 4.0 9.0 16.0; 1.0 4.0 9.0 16.0]
A = QuadraticInterpolation(u,t)

@test A(2.0) == [4.0,4.0]
@test A(1.5) == [2.25,2.25]
@test A(3.5) == [12.25,12.25]
@test A(2.5) == [6.25,6.25]

# Lagrange Interpolation
u = [1.0, 4.0, 9.0]
t = [1.0, 2.0, 3.0]
A = LagrangeInterpolation(u,t)

@test A(2.0) == 4.0
@test A(1.5) == 2.25

u = [1.0, 8.0, 27.0, 64.0]
t = [1.0, 2.0, 3.0, 4.0]
A = LagrangeInterpolation(u,t)

@test A(2.0) == 8.0
@test A(1.5) ≈ 3.375
@test A(3.5) ≈ 42.875

u = [1.0 4.0 9.0 16.0; 1.0 4.0 9.0 16.0]
A = LagrangeInterpolation(u,t)

@test A(2.0) == [4.0,4.0]
@test A(1.5) ≈ [2.25,2.25]
@test A(3.5) ≈ [12.25,12.25]

# Akima Interpolation
u = [0.0, 2.0, 1.0, 3.0, 2.0, 6.0, 5.5, 5.5, 2.7, 5.1, 3.0]
t = collect(0.0:10.0)
A = AkimaInterpolation(u, t)

@test A(0.0) ≈ 0.0
@test A(0.5) ≈ 1.375
@test A(1.0) ≈ 2.0
@test A(1.5) ≈ 1.5
@test A(2.5) ≈ 1.953125
@test A(3.5) ≈ 2.484375
@test A(4.5) ≈ 4.1363636363636366866103344
@test A(5.1) ≈ 5.9803623910336236590978842
@test A(6.5) ≈ 5.5067291516462386624652936
@test A(7.2) ≈ 5.2031367459745245795943447
@test A(8.6) ≈ 4.1796554159017080820603951
@test A(9.9) ≈ 3.4110386597938129327189927
@test A(10.0) ≈ 3.0

@testset "ConstantInterpolation" begin

    t = [1.0, 2.0, 3.0, 4.0]

    @testset "Vector case" for u in
        [[1.0, 2.0, 0.0, 1.0], ["B", "C", "A", "B"]]

        A = ConstantInterpolation(u, t, dir=:right)
        @test A(0.5) == u[1]
        @test A(1.0) == u[1]
        @test A(1.5) == u[2]
        @test A(2.0) == u[2]
        @test A(2.5) == u[3]
        @test A(3.0) == u[3]
        @test A(3.5) == u[1]
        @test A(4.0) == u[1]
        @test A(4.5) == u[1]

        A = ConstantInterpolation(u, t) # dir=:left is default
        @test A(0.5) == u[1]
        @test A(1.0) == u[1]
        @test A(1.5) == u[1]
        @test A(2.0) == u[2]
        @test A(2.5) == u[2]
        @test A(3.0) == u[3]
        @test A(3.5) == u[3]
        @test A(4.0) == u[1]
        @test A(4.5) == u[1]
    end

    @testset "Matrix case" for u in
        [[1.0 2.0 0.0 1.0; 1.0 2.0 0.0 1.0], ["B" "C" "A" "B"; "B" "C" "A" "B"]]

        A = ConstantInterpolation(u, t, dir=:right)
        @test A(0.5) == u[:,1]
        @test A(1.0) == u[:,1]
        @test A(1.5) == u[:,2]
        @test A(2.0) == u[:,2]
        @test A(2.5) == u[:,3]
        @test A(3.0) == u[:,3]
        @test A(3.5) == u[:,1]
        @test A(4.0) == u[:,1]
        @test A(4.5) == u[:,1]

        A = ConstantInterpolation(u, t) # dir=:left is default
        @test A(0.5) == u[:,1]
        @test A(1.0) == u[:,1]
        @test A(1.5) == u[:,1]
        @test A(2.0) == u[:,2]
        @test A(2.5) == u[:,2]
        @test A(3.0) == u[:,3]
        @test A(3.5) == u[:,3]
        @test A(4.0) == u[:,1]
        @test A(4.5) == u[:,1]
    end
end

# QuadraticSpline Interpolation
u = [0.0, 1.0, 3.0]
t = [-1.0, 0.0, 1.0]

A = QuadraticSpline(u,t)

#             Solution ->
#             f(x) = (x+1)^2 for x -> [-1.0, 0.0]
#             f(x) = 1+2x    for x -> [0.0, 1.0]

@test A(-0.5) == 0.25
@test A(0.7) == 2.4
@test A(-1.0) == 0.0
@test A(0.0) == 1.0
@test A(1.0) == 3.0


# CubicSpline Interpolation
u = [0.0, 1.0, 3.0]
t = [-1.0, 0.0, 1.0]

A = CubicSpline(u,t)

#             Solution ->
#             f(x) = 1 + 1.5x + x^2 + 0.5x^3 for x -> [-1.0, 0.0]
#             f(x) = 1 + 1.5x + x^2 - 0.5x^3   for x -> [0.0, 1.0]

@test A(-0.5) == 0.4375
@test A(0.5) == 1.9375
@test A(-0.7) == 0.2685
@test A(0.3) == 1.5265
@test A(-1.0) == 0.0
@test A(0.0) == 1.0
@test A(1.0) == 3.0



# BSpline Interpolation and Approximation
t = [0,62.25,109.66,162.66,205.8,252.3]
u = [14.7,11.51,10.41,14.95,12.24,11.22]

A = BSplineInterpolation(u,t,2,:Uniform,:Uniform)

@test [A(25.0), A(80.0)] == [13.454197730061425, 10.305633616059845]
@test [A(190.0), A(225.0)] == [14.07428439395079, 11.057784141519251]
@test [A(t[1]), A(t[end])] == [u[1], u[end]]

A = BSplineInterpolation(u,t,2,:ArcLen,:Average)

@test [A(25.0), A(80.0)] == [13.363814458968486, 10.685201117692609]
@test [A(190.0), A(225.0)] == [13.437481084762863, 11.367034741256463]
@test [A(t[1]), A(t[end])] == [u[1], u[end]]

A = BSplineApprox(u,t,2,4,:Uniform,:Uniform)

@test [A(25.0), A(80.0)] ≈ [12.979802931218234, 10.914310609953178]
@test [A(190.0), A(225.0)] ≈ [13.851245975109263, 12.963685868886575]
@test [A(t[1]), A(t[end])] ≈ [u[1], u[end]]


# Curvefit Interpolation
rng = StableRNG(12345)
model(x, p) = @. p[1]/(1 + exp(x - p[2]))
t = range(-10, stop=10, length=40)
u = model(t, [1.0, 2.0]) + 0.01*randn(rng, length(t))
p0 = [0.5, 0.5]

A = Curvefit(u, t, model, p0, LBFGS())

ts = [-7.0, -2.0, 0.0, 2.5, 5.0]
vs = [1.0013468217936277, 0.9836755196317837, 0.8833959853995836, 0.3810348276782708, 0.048062978598861855]
us = A.(ts)

@test vs ≈ us

# missing values handling tests
u = [1.0, 4.0, 9.0, 16.0, 25.0, missing, missing]
t = [1.0, 2.0, 3.0, 4.0, missing, 6.0, missing]
A = QuadraticInterpolation(u,t)

@test A(2.0) == 4.0
@test A(1.5) == 2.25
@test A(3.5) == 12.25
@test A(2.5) == 6.25

u = copy(hcat(u, u)')
A = QuadraticInterpolation(u,t)

@test A(2.0) == [4.0, 4.0]
@test A(1.5) == [2.25, 2.25]
@test A(3.5) == [12.25,12.25]
@test A(2.5) == [6.25, 6.25]
