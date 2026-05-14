using DataInterpolations

# Dependent variable
u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]

# Independent variable
t = [0.0, 62.25, 109.66, 162.66, 205.8, 252.3]

t_eval_3 = [50.0, 100.0, 150.0]

# Matrix example (stacked form)
u_matrix = [14.7 11.51 10.41 14.95 12.24 11.22;
            7.35 5.76 5.21 7.48 6.12 5.61]
#A_matrix = CubicSpline(u_matrix, t)
A_matrix = AkimaInterpolation(u_matrix, t)

# Pre-allocate for 3 evaluation points
out_matrix = zeros(2, 3)
A_matrix(out_matrix, t_eval_3)

A_matrix(t_eval_3)

println(out_matrix)


M = [
    1.0 2.0 3.0;
    4.0 5.0 6.0;
    7.0 8.0 9.0
]

M3 = rand(3, 3, 3)


ind = findall(x -> x > 1/2, M3)
