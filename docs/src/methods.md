# Interpolation using different methods

We will use the following data to demonstrate interpolation methods.

```@example tutorial
using DataInterpolations, Plots
gr() # hide

# Dependent variable
u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]

# Independent variable
t = [0.0, 62.25, 109.66, 162.66, 205.8, 252.3]
```

For each method, we will show how to perform the fit and use the plot recipe
to show the fitting curve.

## Linear Interpolation

This is a linear interpolation between the ends points of the interval of input data points.

```@example tutorial
A = LinearInterpolation(u, t)
plot(A)
```

## Quadratic Interpolation

This function fits a parabola passing through the two nearest points from the input
data point as well as the next-closest point on the right or left, depending on
whether the forward- or backward-looking mode is selected (default mode is
forward-looking). It is continuous and piecewise differentiable.

```@example tutorial
A = QuadraticInterpolation(u, t) # same as QuadraticInterpolation(u,t,:Forward)
# alternatively: A = QuadraticInterpolation(u,t,:Backward)
plot(A)
```

## Lagrange Interpolation

It fits a polynomial of degree d (=length(t)-1), and is thus a continuously
differentiable function.

```@example tutorial
A = LagrangeInterpolation(u, t)
plot(A)
```

## Akima Interpolation

This function fits piecewise cubic polynomials which forms a continuously differentiable function.
This differs from Cubic Spline as coefficients are computed using only neighbouring points and hence the
fit looks more natural.

```@example tutorial
A = AkimaInterpolation(u, t)
plot(A)
```

## Constant Interpolation

This function is constant between data points. By default,
it takes the value at the left end of the interval. One can change that behavior by
passing the keyword argument `dir = :right`.

```@example tutorial
A = ConstantInterpolation(u, t)
plot(A)
```

Or using the right endpoints:

```@example tutorial
A = ConstantInterpolation(u, t, dir = :right)
plot(A)
```

## Smoothed Constant Interpolation

This function is much like the constant interpolation above, but the transition
between consecutive values is smoothed out so that the function is continuously
differentiable. The smoothing is done in such a way that the integral of this function
is never much off from the same integral of constant interpolation without smoothing (because of the symmetry of the smoothing sections).
The maximum smoothing distance in the `t` direction from the data points can be set
with the keyword argument `d_max`.

```@example tutorial
A = ConstantInterpolation(u, t)
plot(A)
A = SmoothedConstantInterpolation(u, t; d_max = 10)
plot!(A)
```

Note that `u[end]` is ignored.

## Quadratic Spline

This is the quadratic spline. It is a continuously differentiable interpolation
which hits each of the data points exactly. Splines are a local interpolation
method, meaning that the curve in a given spot is only affected by the points
nearest to it.

```@example tutorial
A = QuadraticSpline(u, t)
plot(A)
```

## Cubic Spline

This is the cubic spline. It is a continuously twice differentiable interpolation
which hits each of the data points exactly.

```@example tutorial
A = CubicSpline(u, t)
plot(A)
```

## B-Splines

This is an interpolating B-spline. B-splines are a global method, meaning
that every data point is taken into account for each point of the curve.
The interpolating B-spline is the version which hits each of the points. This
method is described in more detail [here](https://pages.mtu.edu/%7Eshene/COURSES/cs3621/NOTES/INT-APP/CURVE-INT-global.html).
Let's plot a cubic B-spline (3rd order). Since the data points are not close to
uniformly spaced, we will use the `:ArcLen` and `:Average` choices:

```@example tutorial
A = BSplineInterpolation(u, t, 3, :ArcLen, :Average)
plot(A)
```

The approximating B-spline is a smoothed version of the B-spline. It again is
a global method. In this case, we need to give a number of control points
`length(t)>h` and this method fits a B-spline through the control points which
is a least square approximation. This has a natural effect of smoothing the
data. For example, if we use 4 control points, we get the result:

```@example tutorial
A = BSplineApprox(u, t, 3, 4, :ArcLen, :Average)
plot(A)
```

## Cubic Hermite Spline

This is the cubic (third order) Hermite interpolation. It matches the values and first order derivatives in the data points exactly.

```@example tutorial
du = [-0.047, -0.058, 0.054, 0.012, -0.068, 0.0011]
A = CubicHermiteSpline(du, u, t)
plot(A)
```

## PCHIP Interpolation

This is a type of `CubicHermiteSpline` where the derivative values `du` are derived from the input data in such a way that the interpolation never overshoots the data.

```@example tutorial
A = PCHIPInterpolation(u, t)
plot(A)
```

## Quintic Hermite Spline

This is the quintic (fifth order) Hermite interpolation. It matches the values and first and second order derivatives in the data points exactly.

```@example tutorial
ddu = [0.0, -0.00033, 0.0051, -0.0067, 0.0029, 0.0]
du = [-0.047, -0.058, 0.054, 0.012, -0.068, 0.0011]
A = QuinticHermiteSpline(ddu, du, u, t)
plot(A)
```

## Regularization Smoothing

Smoothing by regularization (a.k.a. ridge regression) finds a function ``\hat{u}``
that minimizes the objective function:

``Q(\hat{u}) = \int_{t_1}^{t_N} |\hat{u}(t) - u(t)|^2 \mathrm{d}t + \lambda \int_{\hat{t}_1}^{\hat{t}_N} |\hat{u}^{(d)}(\hat{t})|^2 \mathrm{d} \hat{t}``

where ``(d)`` denotes derivative order and ``\lambda`` is the regularization
(smoothing) parameter. The integrals are evaluated numerically at the set of
``t`` values for the first term and ``\hat{t}`` values for the second term
(equal to ``t`` if not provided). Regularization smoothing is a global method
that creates a smooth curve directly. See [Stickel (2010)
Comput. Chem. Eng. 34:467](http://dx.doi.org/10.1016/j.compchemeng.2009.10.007)
for details. The implementation in this package uses cubic splines to
interpolate between the smoothed points after they are determined.

```@example tutorial
using RegularizationTools
d = 2
λ = 1e3
A = RegularizationSmooth(u, t, d; λ = λ, alg = :fixed)
û = A.û
# interpolate using the smoothed values
N = 200
titp = collect(range(minimum(t), maximum(t), length = N))
uitp = A.(titp)
lw = 1.5
scatter(t, u, label = "data")
scatter!(t, û, marker = :square, label = "smoothed data")
plot!(titp, uitp, lw = lw, label = "smoothed interpolation")
```

## Dense Data Demonstration

Some methods are better suited for dense data. Let's generate such data to
demonstrate these methods.

```@example tutorial
import StableRNGs: StableRNG
rng = StableRNG(318)
t = sort(10 .* rand(rng, 100))
u = sin.(t) .+ 0.5 * randn(rng, 100);
```

## Regularization Smoothing

Although smoothing by regularization can be used to interpolate sparse data as
shown above, it is especially useful for dense as well as scattered data (unequally
spaced, unordered, and/or repeat-valued). Generalized cross validation (GCV) or
so-called L-curve methods can be used to determine an "optimal" value for the
smoothing parameter. In this example, we perform smoothing in two ways. In the
first, we find smooth values at the original ``t`` values and then
interpolate. In the second, we perform the smoothing for the interpolant
``\hat{t}`` values directly. GCV is used to determine the regularization
parameter for both cases.

```@example tutorial
d = 4
A = RegularizationSmooth(u, t, d; alg = :gcv_svd)
û = A.û
N = 200
titp = collect(range(minimum(t), maximum(t), length = N))
uitp = A.(titp)
Am = RegularizationSmooth(u, t, titp, d; alg = :gcv_svd)
ûm = Am.û
scatter(t, u, label = "simulated data", legend = :top)
scatter!(t, û, marker = (:square, 4), label = "smoothed data")
plot!(titp, uitp, lw = lw, label = "smoothed interpolation")
plot!(titp, ûm, lw = lw, linestyle = :dash, label = "smoothed, more points")
```

## Curve Fits

A curve fit works with both dense and sparse data. We will demonstrate the curve
fit on the dense data since we generated it based on `sin(t)`, so this is the
curve we want to fit through it. To do so, let's define a similar function
with parameters. Let's choose the form:

```@example tutorial
m(t, p) = @. p[1] * sin(p[2] * t) + p[3] * cos(p[4] * t)
```

Notice that this is a function on the whole array of `t` and expects an array
for the predicted `u` out. This choice of `m` is based on the assumption that our
function is of the form `p1*sin(p2*t)+p3*cos(p4*t)`. We want to find the `p` to
match our data. Let's start with the guess of every `p` being zero, that is
`p=ones(4)`. Then we would fit this curve using:

```@example tutorial
using Optim
A = Curvefit(u, t, m, ones(4), LBFGS())
plot(A)
```

We can check what the fitted parameters are via:

```@example tutorial
A.pmin
```

Notice that it essentially made `p3=0` with `p1=p2=1`, meaning it approximately
found `sin(t)`! But note that the ability to fit is dependent on the initial
parameters. For example, with `p=zeros(4)` as the initial parameters, the fit
is not good:

```@example tutorial
A = Curvefit(u, t, m, zeros(4), LBFGS())
plot(A)
```

And the parameters show the issue:

```@example tutorial
A.pmin
```
