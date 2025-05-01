# Extrapolation methods

We will use the following interpolation to demonstrate the various extrapolation methods.

```@example tutorial
using DataInterpolations, Plots

u = [0.86, 0.65, 0.44, 0.76, 0.73]
t = [0.0022, 0.68, 1.41, 2.22, 2.46]
t_eval_left = range(-1, first(t), length = 25)
t_eval_right = range(last(t), 3.5, length = 25)
A = QuadraticSpline(u, t)
plot(A)
```

Extrapolation behavior can be set left and right of the data simultaneously with the `extrapolation` keyword, or left and right separately with the `extrapolation_left` and `extrapolation_right` keywords respectively.

## `ExtrapolationType.None`

This extrapolation type will throw an error when the input `t` is beyond the data in the specified direction.

## `ExtrapolationType.Constant`

This extrapolation type extends the interpolation with the boundary values of the data `u`.

```@example tutorial
A = QuadraticSpline(u, t; extrapolation = ExtrapolationType.Constant)
plot(A)
plot!(t_eval_left, A.(t_eval_left); label = "extrapolation left")
plot!(t_eval_right, A.(t_eval_right); label = "extrapolation right")
```

## `ExtrapolationType.Linear`

This extrapolation type extends the interpolation with a linear continuation of the interpolation, making it $C^1$ smooth at the data boundaries.

```@example tutorial
A = QuadraticSpline(u, t; extrapolation = ExtrapolationType.Linear)
plot(A)
plot!(t_eval_left, A.(t_eval_left); label = "extrapolation left")
plot!(t_eval_right, A.(t_eval_right); label = "extrapolation right")
```

## `ExtrapolationType.Extension`

This extrapolation type extends the interpolation with a continuation of the expression for the interpolation at the boundary intervals for maximum smoothness.

```@example tutorial
A = QuadraticSpline(u, t; extrapolation = ExtrapolationType.Extension)
plot(A)
plot!(t_eval_left, A.(t_eval_left); label = "extrapolation down")
plot!(t_eval_right, A.(t_eval_right); label = "extrapolation up")
```

## `ExtrapolationType.Periodic`

this extrapolation type extends the interpolation such that `A(t + T) == A(t)` for all `t`, where the period is given by `T = last(A.t) - first(A.t)`.

```@example tutorial
T = last(A.t) - first(A.t)
t_eval_left = range(first(t) - 2T, first(t), length = 100)
t_eval_right = range(last(t), last(t) + 2T, length = 100)
A = QuadraticSpline(u, t; extrapolation = ExtrapolationType.Periodic)
plot(A)
plot!(t_eval_left, A.(t_eval_left); label = "extrapolation down")
plot!(t_eval_right, A.(t_eval_right); label = "extrapolation up")
```

## `ExtrapolationType.Reflective`

this extrapolation type extends the interpolation such that `A(t_ + t) == A(t_ - t)` for all `t_, t` such that `(t_ - first(A.t)) % T == 0` and `0 < t < T`, where `T = last(A.t) - first(A.t)`.

```@example tutorial
A = QuadraticSpline(u, t; extrapolation = ExtrapolationType.Reflective)
plot(A)
plot!(t_eval_left, A.(t_eval_left); label = "extrapolation down")
plot!(t_eval_right, A.(t_eval_right); label = "extrapolation up")
```

## Mixed extrapolation

You can also have different extrapolation types left and right of the data.

```@example tutorial
A = QuadraticSpline(u, t; extrapolation_left = ExtrapolationType.Reflective,
    extrapolation_right = ExtrapolationType.Periodic)
plot(A)
plot!(t_eval_left, A.(t_eval_left); label = "extrapolation left")
plot!(t_eval_right, A.(t_eval_right); label = "extrapolation right")
```
