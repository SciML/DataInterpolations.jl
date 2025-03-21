# Smooth arc length interpolation

Arc length interpolation is interpolation between points using a curve that is parameterized by arc length. That is: the curve parameterization has unit speed everywhere, and so the parameter `t` at each point on the curve is equal to the total distance traveled from the beginning of the curve. In this context, by 'smooth' we mean that the curve is continuously differentiable.

## Usage

`DataInteprolations.jl` offers an arc length interpolation method that approximates an existing non arc length interpolation by circle and line segments. This can be done by providing an interpolation object:

```@example tutorial
using DataInterpolations
using Plots
using Random

Random.seed!(2)

# Example from interpolation object
u = cumsum([rand(3) for _ in 1:10])
t = 1:10
A_shape = QuadraticSpline(u, t)
A = SmoothArcLengthInterpolation(A_shape; m = 10)

function plot_itp(itp)
    t_eval = range(itp.t[1], itp.t[end]; length = 1000)
    u_eval = zeros(3, 1000)
    itp(u_eval, t_eval)

    plot(eachrow(u_eval)...; label = "SmoothArcLengthInterpolation")
    scatter!(eachrow(u_eval[:, 1:50:end])...; label = "Equidistant points on the curve")
    scatter!(eachrow(hcat(A.shape_itp.u...))...; label = "Original data")
end

plot_itp(A)
```

Here `m` determines how fine the approximation is. It is also possible to just provide the data points, optionally providing `t` and a preferred interpolation type which determines the shape of the curve.

```@example tutorial
# Example from only u
A = SmoothArcLengthInterpolation(hcat(u...))
plot_itp(A)
```

## Docstrings

To do: add doc strings for the different constructors and add them here

```@docs
SmoothArcLengthInterpolation
```

## Method derivation

Say we have an ordered set of points $u_1, \ldots, u_n \in \mathbb{R}^N$ and we want to make a lightweight $C^1$ smooth interpolation by arc-length $\tilde{\gamma}: [0,T] \rightarrow \mathbb{R}^N$ through these points. The first part is easy, just pick your favorite established interpolation method that achieves $C^1$ smoothness. The arc-length part however turns out to be quite [nasty](https://ijpam.eu/contents/2006-31-3/10/10.pdf). Here I propose a method that is quite general and cheap to compute.

### The 2-dimensional case

2 is the smallest number of dimensions in which the posed problem is non-trivial. Say we use an established (non arc-length) interpolation method for our set of points to obtain the $C^1$ curve

```math
    \gamma : [0, T] \rightarrow \mathbb{R}^2
```

for which

```math
    \gamma(t_i) = u_i \quad i = 1, \ldots, n,
```

given a suitable set of 'time' values

```math
    0 = t_1 < t_2 < \ldots < t_n = T,
```

for instance

```math
    t_i = \sum_{k=1}^{i-1} \|u_{k+1} - u_k\|_2.
```

We now want to approximate $\gamma$ piecewise with sections that are trivially parameterizable by arc-length, namely line segments and circle segments. To do this, we fix some $m \in \mathbb{N}$ and define a refined set of time points $\left(\tilde{t}_j\right)_{j=1}^{m(n-1) + 1}$ given by

```math
    \tilde{t}_{m(k-1) + l} = t_k + \frac{l}{m + 1}(t_{k+1} - t_k), \quad k = 1 \ldots n-1, \; l = 1, \ldots m.
```

In these refined time points we evaluate $\gamma$ and its normalized derivative:

```math
    \tilde{u}_j = \gamma\left(\tilde{t}_j\right), \; \tilde{d}_j = \frac{\dot{\gamma}\left(\tilde{t}_j\right)}{\|\dot{\gamma}\left(\tilde{t}_j\right)\|_2}, \qquad j = 1, \ldots, m(n-1) + 1.
```

As a first step to create the interpolation by arc length $\tilde{\gamma}$, we make a piecewise linear curve which is tangent to $\gamma$ in $\tilde{u}_j$ for each line segment, where we denote the intersection of consecutive tangent lines by $\tilde{u}_{j, \text{int}}$:

```math
    \begin{align*}
    \tilde{u}_{j, \text{int}} &=& \tilde{u}_j + \frac{\langle\tilde{u}_{j+1}-\tilde{u}_j, \tilde{d}_j\rangle - \langle\tilde{d}_j, \tilde{d}_{j+1}\rangle \langle\tilde{u}_{j+1}-\tilde{u}_j, \tilde{d}_{j+1}\rangle}{1 - \langle\tilde{d}_j, \tilde{d}_{j+1}\rangle^2}\tilde{d}_j \\
    &=& \tilde{u}_{j+1} + \frac{\langle\tilde{d}_j, \tilde{d}_{j+1}\rangle\langle\tilde{u}_{j+1}-\tilde{u}_j, \tilde{d}_j\rangle - \langle\tilde{u}_{j+1}-\tilde{u}_j, \tilde{d}_{j+1}\rangle}{1 - \langle\tilde{d}_j, \tilde{d}_{j+1}\rangle^2}\tilde{d}_{j+1}.
    \end{align*}
```

As expected this doesn't work for $\langle\tilde{d}_j, \tilde{d}_{j+1}\rangle^2 = 1$, which means that the consecutive tangent lines are parallel. In fact, in the above equation we want the coefficient of $\tilde{d}_j$ to be positive and the coefficient of $\tilde{d}_{j+1}$ to be negative, to ensure that $\tilde{u}_{j, \text{int}}$ lies properly in between $\tilde{u}_j$ and $\tilde{u}_{j+1}$.

```@setup tutorial
using DataInterpolations
using Plots
using Random
using LinearAlgebra

Random.seed!(4)

n = 4
m = 2
u = [rand(2) for _ in 1:n]
t = cumsum(norm.(diff(u)))
pushfirst!(t, 0)

γ = QuadraticSpline(u, t)

t_eval = range(first(t), last(t), length = 250)
u_eval = hcat(γ.(t_eval)...)

t_tilde = zeros(m * (n - 1) + 1)

for k in 1:(n-1)
    t_tilde[m * (k - 1) + 1 : m * k] = range(t[k], t[k+1], length = m + 1)[1:(end - 1)]
end

t_tilde[end] = t[end]
u_tilde = hcat(γ.(t_tilde)...)

d_tilde = DataInterpolations.derivative.(Ref(γ), t_tilde)
normalize!.(d_tilde)
d_tilde = hcat(d_tilde...)

n_intervals = length(t_tilde) - 1
intersection_points = zeros(2, n_intervals)

u_tilde_j_int_1 = zeros(2)
u_tilde_j_int_2 = zeros(2)

for j in 1:n_intervals
    d_tilde_j = view(d_tilde, :, j)
    d_tilde_j_plus_1 = view(d_tilde, :, j + 1)
    d_inner = dot(d_tilde_j, d_tilde_j_plus_1)
    u_tilde_j = view(u_tilde, :, j)
    u_tilde_j_plus_1 = view(u_tilde, :, j + 1)
    Δu = u_tilde_j_plus_1 - u_tilde_j

    coef_1 = (dot(Δu, d_tilde_j) - d_inner * dot(Δu, d_tilde_j_plus_1))/(1 - d_inner^2)
    coef_2 = (d_inner * dot(Δu, d_tilde_j) - dot(Δu, d_tilde_j_plus_1))/(1 - d_inner^2)

    @. u_tilde_j_int_1 = u_tilde_j + coef_1 * d_tilde_j
    @. u_tilde_j_int_2 = u_tilde_j_plus_1 + coef_2 * d_tilde_j_plus_1
    @assert u_tilde_j_int_1 ≈ u_tilde_j_int_2

    intersection_points[:, j] .= u_tilde_j_int_1
end

function plot_tangent_curve()
    p = plot(; aspect_ratio = :equal, legend = :topleft, title = "m = $m")
  
    # Plot curve γ
    plot!(u_eval[1,:], u_eval[2,:]; label = raw"\gamma")

    # Plot original points
    u_ = hcat(u...)
    scatter!(u_[1,:], u_[2,:], label = raw"$u_i$"; markersize = 6, markerstrokewidth = 0)

    # Plot refined evaluation points
    scatter!(u_tilde[1,:], u_tilde[2,:]; label = raw"$\tilde{u}_j$", markerstrokewidth = 0)

    # Plot tangent curve
    scatter!(intersection_points[1, :], intersection_points[2, :]; 
    markerstrokewidth = 0, markersize = 3, label = raw"$\tilde{u}_{j, \mathrm{int}}$")
    plot!([u_tilde[1, 1], intersection_points[1, :]..., u_tilde[1, end]], 
          [u_tilde[2, 1], intersection_points[2, :]..., u_tilde[2, end]]; 
          label = "Tangent curve")
    p
end
```

```@example tutorial
plot_tangent_curve() # hide
```

As a last step to obtain our curve by arc length $\tilde{\gamma}$ we want to get rid of the kinks in the tangent curve. We do this by replacing sections of the tangent curve by circle arcs. For each $\tilde{u}_{j, \text{int}}$ we compute the shortest distance to the neighboring evaluation points on $\gamma$:

```math
    \delta_j = \min\left\{
        \|\tilde{u}_j - \tilde{u}_{j, \text{int}}\|_2, 
        \|\tilde{u}_{j + 1} - \tilde{u}_{j, \text{int}}\|_2
    \right\}.
```

From this we compute 2 points that are on the tangent curve and equidistant from $\tilde{u}_{j  + \frac{1}{2}}$:

```math
    \tilde{u}_{j, \text{start}} = \tilde{u}_{j, \text{int}} - \delta_j \tilde{d}_j, 
```

```math
    \tilde{u}_{j, \text{end}} = \tilde{u}_{j, \text{int}} + \delta_j \tilde{d}_{j+1}.
```

Note that by this definition

```math
    \tilde{u}_{j, \text{start}} = \tilde{u}_j \quad \vee \quad \tilde{u}_{j, \text{end}} = \tilde{u}_{j  + 1}.
```

Now we can define a circle arc from $\tilde{u}_{j, \text{start}}$ to $\tilde{u}_{j, \text{end}}$ given the center

```math
    c_j = \tilde{u}_{j, \text{int}} + \delta_j\frac{\tilde{d}_{j+1} - \tilde{d}_j}{1 - \langle\tilde{d}_j,\tilde{d}_{j+1}\rangle}
```

and radius

```math
    R_j = \delta_j\sqrt{\frac{1 + \langle\tilde{d}_j,\tilde{d}_{j+1}\rangle}{1 - \langle\tilde{d}_j,\tilde{d}_{j+1}\rangle}}.
```

We obtain the circle arc

```math
    c_j + \cos\left(\frac{t}{R_j}\right)v_{j, 1} + \sin\left(\frac{t}{R_j}\right)v_{j, 2}, \quad t \in [0, \Delta t_{j, \text{arc}}],
```

where

```math
    v_{j, 1} = -\delta_j \frac{\tilde{d}_{j+1} - \langle\tilde{d}_j,\tilde{d}_{j+1}\rangle\tilde{d}_j}{1 - \langle\tilde{d}_j,\tilde{d}_{j+1}\rangle},
    \quad
    v_{j, 2} = R_j \tilde{d}_j.
```

By this definition $\|v_{j, 1}\|_2 = \|v_{j, 2}\|_2 = R_j$ and $\langle v_{j, 1}, v_{j, 2}\rangle = 0$. Furthermore:

```math
    \Delta t_{j, \text{arc}} = R_j\theta_{j, \;\max}= 2R_j \arctan\left(\frac{\delta_j}{R_j}\right).
```

```@setup tutorial
function mark_right_angle!(corner, dir1, dir2; l = 0.05)
    points = zeros(3, 2)
    points[1, :] = corner + l * dir1
    points[2, :] = points[1, :] + l * dir2
    points[3, :] = points[2, :] - l * dir1
    
    plot!(points[:, 1], points[:, 2]; color = :black, label = "")
end

function plot_arc_construction(; indicate_delta = true)
    p = plot(; aspect_ratio = :equal, axis = false)

    tu_j = [0.0, 0.2]
    tu_j_int = [0.5, 0.6]
    tu_j_plus_1 = [1.2, 0.3]
    δⱼ = norm(tu_j_int - tu_j)
    tu_j_start = tu_j
    tu_j_end = tu_j_int + δⱼ * (tu_j_plus_1 - tu_j_int) / norm(tu_j_plus_1 - tu_j_int)

    td_j = tu_j_int - tu_j
    normalize!(td_j)
    td_j_plus_1 = tu_j_plus_1 - tu_j_int
    normalize!(td_j_plus_1)
    inner = dot(td_j, td_j_plus_1)
    Δtd_j = td_j_plus_1 - td_j

    tc_j = tu_j_int + δⱼ / (1 - inner) * Δtd_j

    Rⱼ = δⱼ * sqrt((1 + inner)/(1 - inner))
    vⱼ₁ = -δⱼ * (td_j_plus_1 - inner * td_j)/(1 - inner)
    vⱼ₂ = Rⱼ * td_j
    Δt_j_arc = 2 * Rⱼ * atan(δⱼ, Rⱼ)

    T_ = range(0, π/2; length = 100)
    x_arc = @. tc_j[1] + cos(T_) * vⱼ₁[1] + sin(T_) * vⱼ₂[1]
    y_arc = @. tc_j[2] + cos(T_) * vⱼ₁[2] + sin(T_) * vⱼ₂[2]
    plot!(x_arc, y_arc;  label = "", color = :gray, ls = :dash)
    
    T = range(0, Δt_j_arc, length = 100)
    X_arc = @. tc_j[1] + cos(T/Rⱼ) * vⱼ₁[1] + sin(T/Rⱼ) * vⱼ₂[1]
    Y_arc = @. tc_j[2] + cos(T/Rⱼ) * vⱼ₁[2] + sin(T/Rⱼ) * vⱼ₂[2]
    plot!(X_arc, Y_arc;  label = "", color = :black, linewidth = 2)

    x_arc = @. tc_j[1] + cos(T/Rⱼ) * vⱼ₁[1] / 7 + sin(T/Rⱼ) * vⱼ₂[1] / 7
    y_arc = @. tc_j[2] + cos(T/Rⱼ) * vⱼ₁[2] / 7 + sin(T/Rⱼ) * vⱼ₂[2] / 7
    plot!(x_arc, y_arc;  label = "", color = :black)

    z = tc_j + vⱼ₂

    annotate!(tu_j...,  "\n\n\n" * raw"$\tilde{u}_j = \tilde{u}_{j, \mathrm{start}}$")
    annotate!(tu_j_int..., raw"$\tilde{u}_{j, \mathrm{int}}$" * "\n\n\n")
    annotate!(tu_j_end...,  raw"$\tilde{u}_{j, \mathrm{end}}$" * "\n\n\n")
    annotate!(tu_j_plus_1..., "\n\n\n" * raw"$\tilde{u}_{j + 1}$")
    annotate!((tu_j + td_j/5)..., raw"$\tilde{d}_j$" * "         ")
    annotate!((tu_j_plus_1 + td_j_plus_1/5)..., raw"$\tilde{d}_{j+1}$" * "\n\n")
    annotate!(tc_j..., "       " * raw"$\tilde{c}_j$")
    annotate!(tc_j + vⱼ₁/2..., raw"$v_{j,1}$" * "     ")
    annotate!(tc_j + vⱼ₂/2..., "\n     " * raw"$v_{j,2}$")
    annotate!(tc_j..., raw"$\theta_{j, \max}$" * "  \n\n\n\n\n")
    indicate_delta && annotate!((tu_j + tu_j_int)/2..., raw"$\delta_j$" * "\n\n")

    mark_right_angle!(tu_j_start, td_j, normalize(tc_j - tu_j_start))
    mark_right_angle!(tu_j_end, td_j_plus_1, normalize(tc_j - tu_j_end))
    mark_right_angle!(tc_j, normalize(vⱼ₁), normalize(vⱼ₂))

    # u connections
    points = hcat(tu_j, tu_j_int, tu_j_end, tu_j_plus_1)
    plot!(points[1, :], points[2, :]; marker = :circle, c = :black, ls = :dash, label = "")

    # line segment
    plot!([tu_j_plus_1[1], tu_j_end[1]], [tu_j_plus_1[2], tu_j_end[2]], c = :black, label = "", linewidth = 2)

    # td_j and td_j_plus_1
    plot!([tu_j[1], tu_j[1] + td_j[1]/5], [tu_j[2], tu_j[2] + td_j[2]/5]; arrow=(:closed, 2.0), color = :black, label = "")
    plot!([tu_j_plus_1[1], tu_j_plus_1[1] + td_j_plus_1[1]/5], [tu_j_plus_1[2], tu_j_plus_1[2] + td_j_plus_1[2]/5]; arrow=(:closed, 2.0), color = :black, label = "")

    # Circle segment radii
    points = hcat(tu_j_start, tc_j, tu_j_end)
    plot!(points[1, :], points[2, :]; color = :gray, label = "", ls = :dash)

    # v₂ⱼ
    points = hcat(tc_j, tc_j + vⱼ₂)
    plot!(points[1, :], points[2, :]; color = :gray, label = "", ls = :dash)

    # tc_j
    scatter!([tc_j[1]], [tc_j[2]]; color = :gray, label = "")
    ylims!(-0.65, 0.8)
    return p, (; tu_j, tu_j_plus_1, tu_j_int, tu_j_start, tu_j_end, td_j, td_j_plus_1, inner, tc_j, δⱼ, Rⱼ)
end
```

```@example tutorial
plot_arc_construction()[1] # hide
```

```@setup tutorial
p = plot_tangent_curve()

δ = [min(
    norm(intersection_points[:, j] - u_tilde[:, j]), 
    norm(intersection_points[:, j] - u_tilde[:, j+1])
    ) 
    for j in 1:n_intervals
]

u_tilde_start = zeros(2, n_intervals)
u_tilde_end = zeros(2, n_intervals)

for j in 1:n_intervals
    u_tilde_intⱼ = intersection_points[:, j]
    u_tildeⱼ = u_tilde[:, j]
    u_tildeⱼ₊₁ = u_tilde[:, j + 1]
    d_tildeⱼ = d_tilde[:, j]
    d_tildeⱼ₊₁ = d_tilde[:, j + 1]

    u_tilde_start[:, j] = u_tilde_intⱼ - δ[j] * d_tildeⱼ
    u_tilde_end[:, j]   = u_tilde_intⱼ + δ[j] * d_tildeⱼ₊₁
end

scatter!(u_tilde_start[1, :], u_tilde_start[2, :]; markersize = 3, markerstrokewidth = 0, 
    label = raw"$\tilde{u}_{j, start}$")
scatter!(u_tilde_end[1, :], u_tilde_end[2, :]; markersize = 3, markerstrokewidth = 0,
    label = raw"$\tilde{u}_{j, end}$")

origins = zeros(2, n_intervals)

for j in 1:n_intervals
    u_tilde_intⱼ = intersection_points[:, j]
    d_tildeⱼ = d_tilde[:, j]
    d_tildeⱼ₊₁ = d_tilde[:, j + 1]
    inner = dot(d_tildeⱼ, d_tildeⱼ₊₁)

    origins[:, j] = u_tilde_intⱼ + δ[j] / (1 - inner) * (d_tildeⱼ₊₁ - d_tildeⱼ)

    plot!([u_tilde_start[1, j], origins[1, j], u_tilde_end[1, j]], 
          [u_tilde_start[2, j], origins[2, j], u_tilde_end[2, j]]; 
          label = (j ==1) ? "radii of circle arc" : "", ls = :dash, c = :gray)

    Rⱼ = δ[j] * sqrt((1 + inner)/(1 - inner))
    v₁ = -δ[j] * (d_tildeⱼ₊₁ - inner * d_tildeⱼ) / (1 - inner)
    v₂ = Rⱼ * d_tildeⱼ
    Δt = 2 * Rⱼ * atan(δ[j], Rⱼ)
    T = range(0, Δt, length = 25)
    x = @. origins[1, j] + cos(T/Rⱼ) * v₁[1] + sin(T/Rⱼ) * v₂[1]
    y = @. origins[2, j] + cos(T/Rⱼ) * v₁[2] + sin(T/Rⱼ) * v₂[2]
    plot!(x,y; c = :green, label = (j == 1) ? "circle arc" : "")
end

scatter!(origins[1, :], origins[2, :]; label = raw"$\tilde{c}_j$", c= :gray, markerstrokewidth = 0)
```

```@example tutorial
p # hide
```

That's pretty neat, but this method does not directly generalize to higher dimensional spaces. That is because in general the intersection points $\tilde{u}_{j, \text{int}}$ of the tangent lines do not exist.

### The higher dimensional case

Let's try to generalize the method above. The goal is to find a point $\tilde{u}_{j + \frac{1}{2}}$ and unit direction $\tilde{d}_{j + \frac{1}{2}}$ to add to the tangent curve between $\tilde{u}_j$ and $\tilde{u}_{j+1}$ such that:

  - the tangent line intersections $\tilde{u}_{j, \text{int left}}, \tilde{u}_{j, \text{int right}}$ exist. This means that the new line is fixed by these 2 points;
  - constructing $\tilde{\gamma}$ including this point gives gives an identical result to constructing $\tilde{\gamma}$ excluding this point if the tangent line intersection already existed. The latter implies that $\tilde{u}_{j + \frac{1}{2}}$ and $\tilde{d}_{j + \frac{1}{2}}$ yield a tangent line to the constructed circle arc.

Let's assume the tangent line intersection exists, and we define

```math
    \tilde{u}_{j, \text{int left}} = \tilde{u}_{j, \text{int}} - \delta_j^* \tilde{d}_{j},
```

```math
    \tilde{u}_{j, \text{int right}} = \tilde{u}_{j, \text{int}} + \delta_j^* \tilde{d}_{j+1}.
```

It turns out that if we then let

```math
\delta^*_j = \delta_j \frac{2 - \sqrt{2 + 2 \langle\tilde{d}_j, \tilde{d}_{j+1}\rangle}}{1 - \langle\tilde{d}_j, \tilde{d}_{j+1}\rangle},
```

The line between $\tilde{u}_{, \text{int left}}$ and $ \tilde{u}_{, \text{int right}}$ touches the circle arc as constructed before. It follows that

```math
    \tilde{u}_{j + \frac{1}{2}} = \frac{1}{2}\left[\tilde{u}_{j, \text{int left}} + \tilde{u}_{j, \text{int right}}\right],
    \qquad
    \tilde{d}_{j + \frac{1}{2}} = \frac{\tilde{u}_{j, \text{int right}} - \tilde{u}_{j, \text{int left}}}{\|\tilde{u}_{j, \text{int right}} - \tilde{u}_{j, \text{int left}}\|_2}.
```

```@setup tutorial
p, vars = plot_arc_construction(; indicate_delta = false)
xlims!(0.0, 1.3)
ylims!(0.0, 0.8)
δⱼ_star = vars.δⱼ * (2 - sqrt(2 + 2 * vars.inner)) / (1 - vars.inner)
tu_j_int_left = vars.tu_j_int - δⱼ_star * vars.td_j
tu_j_int_right = vars.tu_j_int + δⱼ_star * vars.td_j_plus_1
tu_j_plus_half = vars.tu_j_int + δⱼ_star / 2 * (vars.td_j_plus_1 - vars.td_j)
td_j_plus_half = (vars.td_j_plus_1 + vars.td_j) / sqrt(2 + 2 * vars.inner)

scatter!([tu_j_int_left[1], tu_j_plus_half[1], tu_j_int_right[1]], 
         [tu_j_int_left[2], tu_j_plus_half[2], tu_j_int_right[2]]; color = :black, label = "")
annotate!(tu_j_int_left..., raw"$\tilde{u}_{j, \mathrm{int \; left}}$" * "\n\n")
annotate!(tu_j_int_right..., raw"$\tilde{u}_{j, \mathrm{int \; right}}$" * "\n\n")
annotate!(tu_j_plus_half..., "\n\n" * raw"$\tilde{u}_{j + \frac{1}{2}}$")
annotate!(tu_j_plus_half..., "         " * raw"$\tilde{d}_{j + \frac{1}{2}}$" * "\n\n\n")
annotate!((vars.tu_j_int + δⱼ_star / 2 * normalize(vars.tu_j - vars.tu_j_int))..., raw"$\delta^*_j$" * "\n\n")
plot!(
    [vars.tu_j_int[1] - δⱼ_star * vars.td_j[1], vars.tu_j_int[1] + δⱼ_star * vars.td_j_plus_1[1]],
    [vars.tu_j_int[2] - δⱼ_star * vars.td_j[2], vars.tu_j_int[2] + δⱼ_star * vars.td_j_plus_1[2]];
    color = :black, ls = :dash, label = "")
plot!([tu_j_plus_half[1], tu_j_plus_half[1] + td_j_plus_half[1]/5], 
      [tu_j_plus_half[2], tu_j_plus_half[2] + td_j_plus_half[2]/5]; 
       arrow=(:closed, 2.0), color = :black, label = "")
```

```@example tutorial
p # hide
```

If we generalize the definition of $\tilde{u}_{j, \text{int}}$ then we can compute  $\tilde{u}_{j + \frac{1}{2}}$ and $\tilde{d}_{j + \frac{1}{2}}$ as above. Something we can always compute are the points on the tangent lines which are closest together, given by:

```math
    \argmin_{s,\; t \;\in\; \mathbb{R}} \|\tilde{u}_{j+1} + s\tilde{d}_{j+1} - (\tilde{u}_j + t\tilde{d}_j)\|_2.
```

This yields

```math
    \tilde{u}_{j, \text{close left}} = \tilde{u}_j + \frac{\langle\tilde{u}_{j+1}-\tilde{u}_j, \tilde{d}_j\rangle - \langle\tilde{d}_j, \tilde{d}_{j+1}\rangle \langle\tilde{u}_{j+1}-\tilde{u}_j, \tilde{d}_{j+1}\rangle}{1 - \langle\tilde{d}_j, \tilde{d}_{j+1}\rangle^2}\tilde{d}_j,
```

```math
    \tilde{u}_{j, \text{close right}} = \tilde{u}_{j+1} + \frac{\langle\tilde{d}_j, \tilde{d}_{j+1}\rangle\langle\tilde{u}_{j+1}-\tilde{u}_j, \tilde{d}_j\rangle - \langle\tilde{u}_{j+1}-\tilde{u}_j, \tilde{d}_{j+1}\rangle}{1 - \langle\tilde{d}_j, \tilde{d}_{j+1}\rangle^2}\tilde{d}_{j+1}.
```

This is the same as the two expressions for $\tilde{u}_{j, \text{int}}$ from before, except now these expressions aren't necessarily equal. We define $\tilde{u}_{j, \text{int}}$ as the average of these expressions:

```math
    \tilde{u}_{j, \text{int}} = \frac{\tilde{u}_{j, \text{close left}} + \tilde{u}_{j, \text{close right}}}{2}.
```

From this $\delta_j$ and $\delta_j^*$ follow, and

```math
    \tilde{u}_{j, \text{int left}} = \tilde{u}_{j, \text{close left}} - \delta_j^* \tilde{d}_{j},
```

```math
    \tilde{u}_{j, \text{int right}} = \tilde{u}_{j, \text{close right}} + \delta_j^* \tilde{d}_{j+1}.
```
