# Define the function f(R) based on the equation in the image
function f(R, v)
    numerator = R * exp(1 / (3 - 50 * R))
    denominator = (2 + (200 * R) / 3)^(4 / 9) * (2 - (100 * R) / 3)^(5 / 9)
    return v - (numerator / denominator)
end

# Define the derivative of the function f(R) with respect to R
function df(R, v)
    # Approximate the derivative using central difference method
    h = 1e-6
    return (f(R + h, v) - f(R - h, v)) / (2 * h)
end

# Newton's method implementation
function newton(f, df, R0, v; tol=1e-6, maxiter=1000)
    R = R0
    for i in 1:maxiter
        R_new = R - f(R, v) / df(R, v)
        if abs(R_new - R) < tol
            return R_new
        end
        R = R_new
    end
    error("Newton's method did not converge")
end

# Array to store final converged values of R
# This is the value along the intial hypersurface
R_0 = Float64[] 
v = LinRange(0.0000, 150, 110001)
# Initial conditions
let R_initial = 0.025  # Initial guess for R
    for v in v
    R_solution = newton(f, df, R_initial, v)
    R_initial = R_solution
    push!(R_0, R_solution) # Update R0 to the last converged value of R
    end
end
using PlotlyJS
plot(R_0, v)

function Initial_Function(R)
    H_0 = R^2 * (0.0525 - R)^2
end

H_initial = Float64[]
for R in R_0
    if R <= 0.0525
        H = Initial_Function(R)
        push!(H_initial, H)
    end
    if R > 0.0525
        H = 0
        push!(H_initial, H)
    end
end
using PlotlyJS
plot(R_0, H_initial)
# H, R, v are defined
# We run the RK4 scheme to calculate the value of R across the grid.

u = range(-100, 10, length = 110001)
v_range = range(0,150,length = 110001)
using Distributed
addprocs(4)
@everywhere using SharedArrays
R_grid = SharedArray{Float64}(110001,110001)
R_grid[1,:] = R_0
# writing the ODE function
function dR_du(u, R)
    g = (R^2/6)*(u*R + (u^2 * R^2)/9)
    return (- 0.5 * g)
end

#RK4 Steps
function rk4(u, R, du)
    k1 = dR_du(u, R)
    k2 = dR_du(u + 0.5*du, R+0.5*du*k1)
    k3 = dR_du(u + 0.5*du, R + 0.5*du*k2)
    k4 = dR_du(u + du, R+ du*k3)
    return R + (du/6) * (k1 + 2*k2 + 2*k3 + k4)
end
Threads.@threads for j in 1:110001 
    for i in 2:110001
        u_vals = u[i-1]
        du = u[i] - u_vals
        R_grid[j, i] = rk4(u_vals, R_grid[j, i-1], du)
    end
end