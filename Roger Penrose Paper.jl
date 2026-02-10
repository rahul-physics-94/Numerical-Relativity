#Setting up all the required functions
###########################################
# 1. To calculate the absolute value of v
function f(R, v)
    numerator = R * exp(1 / (3 - 50 * R))
    denominator = (2 + (200 * R) / 3)^(4 / 9) * (2 - (100 * R) / 3)^(5 / 9)
    return v - (numerator / denominator)
end

# 2. Finding the derivative using the central difference method. Absolute derivative might make the equations not converge.
function df(R, v)
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

#3. Finding the value of g(u,R)
function g(u,R)
    return (1/6)*u*R^3 + (1/54)* u^2*R^4
end

#4. Finding the derivatives of g wrt R
function dg_dR(u, R)
    return 0.5*u*R^2 + (1/9)*u^2^R^3
end

#5. Finding the second derivative of g wrt R

function dg2_dR2(u, R)
    return u*R + (1/3)*u^2*R^2
end
#6. Derivative of R wrt u for the Runge Kutta scheme

function dR_du(u,R)
    return 0.5*g(u,R)
end

#7. Finding the initial value of the scalar function along R0

function Initial_Function(R)
    H_0 = R^2 * (0.0525 - R)^2
end

####################################################
#Creating the grid using values of u and v
u = LinRange(-100, 10, 20000)
v = LinRange(0,150,20000)

#Finding the values of R along the grid using Newton Method and Runge Kutta scheme
R_0 = Float64[]
let R_initial = 0.025
    for v in v
        R_solution = newton(f, df, R_initial, v)
        R_initial = R_solution
        push!(R_0, R_solution) # Update R0 to the last converged value of R
    end
end
using PlotlyJS
plot(R_0, v)

#Finding the values of the scalar function along the initial grid point

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
GC.gc()
R_grid = Array{Float64}(undef,20000,20000)
R_grid[:,1] .= R_0

#RK4 Steps
function rk4(u, R, du)
    k1 = dR_du(u, R)
    k2 = dR_du(u + 0.5*du, R+0.5*du*k1)
    k3 = dR_du(u + 0.5*du, R + 0.5*du*k2)
    k4 = dR_du(u + du, R+ du*k3)
    return R + (du/6) * (k1 + 2*k2 + 2*k3 + k4)
end
for j in 1:20000
    for i in 2:20000
        u_vals = u[i-1]
        du = u[i] - u_vals
        R_grid[j, i] = rk4(u_vals, R_grid[j, i-1], du)
    end
end