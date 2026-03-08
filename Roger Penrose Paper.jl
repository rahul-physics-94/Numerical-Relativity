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
    return (-f(R + 2*h, v) + (8*f(R + h, v))-(8*f(R - h,v))+f(R - 2*h, v)) / (12 * h)
end

# Newton's method implementation
function newton(f, df, R0, v; tol=1e-6, maxiter=10000)
    R = R0
    for i in 1:maxiter
        R_new = R - (f(R, v) / df(R, v))
        if abs(R_new - R) < tol
            return R_new
        end
        R = R_new
    end
    error("Newton's method did not converge")
end

#3. Finding the value of g(u,R)
@inline function g(u,R)
    return (1/6)*u*R^3 + (1/54)* u^2*R^4
end

#4. Finding the derivatives of g wrt R
@inline function dg_dR(u, R)
    return 0.5*u*R^2 + (1/9)*u^2*R^3
end

#5. Finding the second derivative of g wrt R

@inline function dg2_dR2(u, R)
    return u*R + (1/3)*u^2*R^2
end
#6. Derivative of R wrt u for the Runge Kutta scheme

function dR_du(u,R)
    return -0.5*g(u,R)
end

#7. Finding the initial value of the scalar function along R0

function Initial_Function(u0,R)
    H_0 = (R^2 * (0.0525 - R)^2)*10^8
    q_0 = (2*R*(0.0525)^2 - 6*0.0525*R^2 + 4*R^3)*10^8
    f = (q_0*dg_dR(u0,R) + H_0*(0.5*dg2_dR2(u0,R) - 1))*10^8
    return H_0, q_0, f
end

#RK4 Steps
@inline function rk4(u, R, du)
    k1 = dR_du(u, R)
    k2 = dR_du(u + 0.5*du, R+0.5*du*k1)
    k3 = dR_du(u + 0.5*du, R + 0.5*du*k2)
    k4 = dR_du(u + du, R+ du*k3)
    return R + (du/6) * (k1 + 2*k2 + 2*k3 + k4)
end

####################################################
#Creating the grid using values of u and v
# u = Vector{Float32}(LinRange(-100, 10, 22000))
# v = Vector{Float32}(LinRange(0,550,22000))

# #Finding the values of R along the grid using Newton Method and Runge Kutta scheme
# R_0 = Float32[]
# let R_initial = 0
#     for v in v
#         R_solution = newton(f, df, R_initial, v)
#         R_initial = R_solution
#         push!(R_0, R_solution) # Update R0 to the last converged value of R
#     end
# end
# #using PlotlyJS
# #PlotlyJS.plot(R_0, v)

# #Finding the values of the scalar function along the initial grid point

# H_initial = Float32[]
# q_inital = Float32[]
# f_initial = Float32[]
# for R in R_0
#     u0 = -100
#     if R <= 0.0525
#         H,q,f = Initial_Function(u0,R)
#         push!(H_initial, H)
#         push!(q_inital, q)
#         push!(f_initial, f)
#     else
#         push!(H_initial, 0)
#         push!(q_inital, 0)
#         push!(f_initial, 0)
#     end
# end
        
# using PlotlyJS
# display(PlotlyJS.plot(R_0, H_initial))
# # H, R, v are defined
# # We run the RK4 scheme to calculate the value of R across the grid.
# GC.gc()
# R_grid = Array{Float32}(undef,22000,22000)
# R_grid[:,1] = R_0



# @inbounds for i in 2:22000
#     du = u[i] - u[i-1]
#     @inbounds for j in 1:22000
#         R_grid[j, i] = rk4(u[i-1], R_grid[j, i-1], du)
#     end
# end

# #R_trace = @view R_grid[10000,:]
# #using PlotlyJS
# #plt = PlotlyJS.plot(
#  #   u,
#   #  R_trace
# #)
# #display(plt)

# # Writing the loop for calculating the value of the scalar field.
# H_grid = Array{Float64}(undef, 22000, 22000)
# H_grid[:,1] = H_initial
# q_grid = Array{Float64}(undef, 22000,22000)
# q_grid[:,1] = q_inital
# f_grid = Array{Float64}(undef,22000,22000)
# f_grid[:,1] = f_initial
# du = diff(u)
# @inbounds for i in 2:22000
#     @inbounds for j in 1:22000
#         q_grid[j,i] = q_grid[j,i-1] + f_grid[j,i-1]*du[i-1]
#     end
#     H_grid[1,i] = 0
#     f_grid[1,i] = 0

#     @inbounds for j in 2:22000
#         H_grid[j,i] = H_grid[j-1,i] + q_grid[j-1,i]*(R_grid[j,i]- R_grid[j-1,i])
#         f_grid[j,i] = q_grid[j,i]*dg_dR(u[i],R_grid[j,i]) + H_grid[j,i]*((0.5*dg2_dR2(u[i], R_grid[j, i]))-1)
#     end
# end


# index = findall(x -> abs(x) != 0, H_grid)

# index_2 = findall(x -> abs(x) > 0.0525, R_grid)

# index_3 = findfirst(x -> x >= 2, u)
# H_slice = H_grid[:,2]
# R_slice = R_grid[:,2]


# using PlotlyJS
# slice_plot = PlotlyJS.plot(R_slice, H_slice,
#      xlabel="R",
#      ylabel="H",
# )  
# display(slice_plot)

#using Plots; pyplot();
#plot(u,R_grid,H_grid,st=:surface,camera=(-30,30)
#####################################################################################################################################################################
#Rewriting the loop
#This includes working inside the loop with all the values and having to not store the whole grid values
#but slice them at specific point

using PlotlyJS
function Penrose_Smith()
    du = 0.001
    u = round.([(-100.0 + (i-1)*du) for i in 1:110001], digits = 3)
    dv = 0.005
    v = [(0 + (i-1)*dv) for i in 1:110001]
    uslice = [-100.0, -8.0, -0.340, 0.0, 1.200, 10.0]
    R_0 = Float32[]
    u0 = -100
    v0 = 0
    let R_initial = 0
        for i in v
            R_solution = newton(f, df, R_initial, i)
            R_initial = R_solution
            push!(R_0, R_solution)
        end
    end
    R_prev = zeros(Float32, 110001)
    R_curr = zeros(Float32, 110001)
    H_prev = zeros(Float32, 110001)
    H_curr = zeros(Float32, 110001)
    q_prev = zeros(Float32, 110001)
    q_curr = zeros(Float32, 110001)
    f_prev = zeros(Float32, 110001)
    f_curr = zeros(Float32, 110001)
    @inbounds for i in 1:110001
        @inbounds for j in 1:110001
            if i == 1
                if R_0[j] <= 0.0525
                    H, q, f= Initial_Function(u0, R_0[j])
                    H_prev[j] = H
                    q_prev[j] = q
                    f_prev[j] = f
                else
                    H_prev[j] = 0
                    q_prev[j] = 0
                    f_prev[j] = 0
                end
                R_prev[j] = R_0[j]
            end
            if i > 1
                R_curr[j] = rk4(u[i-1], R_prev[j], du)
                q_curr[j] = q_prev[j] + (f_prev[j]*du)
                H_curr[j] = H_prev[j] + (q_prev[j]*(R_curr[j] - R_prev[j]))
                f_curr[j] = q_curr[j]*dg_dR(u[i], R_curr[j]) + (H_curr[j]*(dg2_dR2(u[i], R_curr[j]) - 1))
                
            end
        end
        if u[i] in uslice
            display(plot(R_prev, H_prev, xlabel="R", ylabel="H", title="u = $(u[i])"))
        end
        copy!(R_prev, R_curr)
        copy!(H_prev,H_curr)
        copy!(q_prev,q_curr)
        copy!(f_prev,f_curr)
        println(i)
    end
    return H_curr, q_curr, R_curr, f_curr
end  

solution = Penrose_Smith()

# using PlotlyJS

# R_prev = zeros(Float64, 1100001)
# R_curr = zeros(Float64, 1100001)
# H_prev = zeros(Float64, 1100001)
# H_curr = zeros(Float64, 1100001)
# q_prev = zeros(Float64, 1100001)
# q_curr = zeros(Float64, 1100001)
# f_prev = zeros(Float64, 1100001)
# f_curr = zeros(Float64, 1100001)

# dv = 0.0005
# v = [(0 + (i-1)*dv) for i in 1:1100001]
# R_0 = Float64[]
# let R_initial = 0
#     for i in v
#         R_solution = newton(f, df, R_initial, i)
#         R_initial = round(R_solution, digits=4)
#         R_solution = round(R_solution, digits=4)
#         push!(R_0, R_solution)
#     end
# end
# R_prev = R_0
# @inbounds for j in 1:1100001
#     if R_0[j] <= 0.0525
#         H, q, f= Initial_Function(-100, R_0[j])
#         H_prev[j] = H
#         q_prev[j] = q
#         f_prev[j] = f
#     else
#         H_prev[j] = 0
#         q_prev[j] = 0
#         f_prev[j] = 0
#     end
# end
# display(PlotlyJS.plot(R_prev, H_prev))