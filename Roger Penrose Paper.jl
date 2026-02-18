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
    H_0 = R^2 * (0.0525 - R)^2
    q_0 = 2*R*(0.0525)^2 - 6*0.0525*R^2 + 4*R^3
    f = q_0*dg_dR(u0,R) + H_0*(0.5*dg2_dR2(u0,R) - 1)
    return H_0, q_0, f
end

####################################################
#Creating the grid using values of u and v
u = Vector{Float32}(LinRange(-100, 10, 22000))
v = Vector{Float32}(LinRange(0,550,22000))

#Finding the values of R along the grid using Newton Method and Runge Kutta scheme
R_0 = Float32[]
let R_initial = 0
    for v in v
        R_solution = newton(f, df, R_initial, v)
        R_initial = R_solution
        push!(R_0, R_solution) # Update R0 to the last converged value of R
    end
end
#using PlotlyJS
#PlotlyJS.plot(R_0, v)

#Finding the values of the scalar function along the initial grid point

H_initial = Float32[]
q_inital = Float32[]
f_initial = Float32[]
for R in R_0
    u0 = -100
    if R <= 0.0525
        H,q,f = Initial_Function(u0,R)
        push!(H_initial, H)
        push!(q_inital, q)
        push!(f_initial, f)
    end
    if R > 0.0525
        H = 0
        push!(H_initial, H)
        q = 0
        push!(q_inital,q)
        f = 0
        push!(f_initial,f)
    end
end
using PlotlyJS
display(PlotlyJS.plot(R_0, H_initial))
# H, R, v are defined
# We run the RK4 scheme to calculate the value of R across the grid.
GC.gc()
R_grid = Array{Float32}(undef,22000,22000)
R_grid[:,1] = R_0

#RK4 Steps
@inline function rk4(u, R, du)
    k1 = dR_du(u, R)
    k2 = dR_du(u + 0.5*du, R+0.5*du*k1)
    k3 = dR_du(u + 0.5*du, R + 0.5*du*k2)
    k4 = dR_du(u + du, R+ du*k3)
    return R + (du/6) * (k1 + 2*k2 + 2*k3 + k4)
end

@inbounds for i in 2:22000
    du = u[i] - u[i-1]
    @inbounds for j in 1:22000
        R_grid[j, i] = rk4(u[i-1], R_grid[j, i-1], du)
    end
end

#R_trace = @view R_grid[10000,:]
#using PlotlyJS
#plt = PlotlyJS.plot(
 #   u,
  #  R_trace
#)
#display(plt)

# Writing the loop for calculating the value of the scalar field.
H_grid = Array{Float64}(undef, 22000, 22000)
H_grid[:,1] = H_initial
q_grid = Array{Float64}(undef, 22000,22000)
q_grid[:,1] = q_inital
f_grid = Array{Float64}(undef,22000,22000)
f_grid[:,1] = f_initial

@inbounds for i in 2:22000
    du = u[i] - u[i-1]
    @inbounds for j in 1:22000
        q_grid[j,i] = q_grid[j,i-1] + f_grid[j,i-1]*du
        H_grid[j,i] = H_grid[j,i-1] + q_grid[j,i-1]*(R_grid[j,i]- R_grid[j,i-1])
        f_grid[j,i] = q_grid[j,i]*dg_dR(u[i],R_grid[j,i]) + H_grid[j,i]*((0.5*dg2_dR2(u[i], R_grid[j, i]))-1)
    end
end


index = findall(x -> abs(x) != 0, H_grid)

index_2 = findall(x -> abs(x) > 0.0525, R_grid)

index_3 = findall(x -> x <= -8, u)
H_slice = H_grid[1:21000,21000]
R_slice = R_grid[1:21000,21000]


using PlotlyJS
slice_plot = PlotlyJS.plot(R_slice, H_slice,
     xlabel="R",
     ylabel="H",
)  
display(slice_plot)

#using Plots; pyplot();
#plot(u,R_grid,H_grid,st=:surface,camera=(-30,30)

#Rewriting the loop
#This includes working inside the loop with all the values and having to not store the whole grid values
#but slice them at specific points
using PlotlyJS
function Penrose_Smith()
    du = 0.0001
    dv = 0.0005
    u_slice = [-8, -0.34, 0, 1.2, 10]
    R_0 = Float32[]
    let R_initial = 0
        for v in v
            R_solution = newton(f, df, R_initial, v)
            R_initial = R_solution
            push!(R_0, R_solution) # Update R0 to the last converged value of R
        end
    end
    H_initial = Float32[]
    q_inital = Float32[]
    f_initial = Float32[]
    for R in R_0
        u0 = -100
        if R <= 0.0525
            H,q,f = Initial_Function(u0,R)
            push!(H_initial, H)
            push!(q_inital, q)
            push!(f_initial, f)
        end
        if R > 0.0525
            H = 0
            push!(H_initial, H)
            q = 0
            push!(q_inital,q)
            f = 0
            push!(f_initial,f)
        end
    end
    plot_1 = display(PlotlyJS.plot(R_0, H_initial))
    #Filling up the R_grid
    R_grid = Array{Float32}(undef,1100001,1100001)
    R_grid[:,1] = R_0

    #RK4 Steps
    @inline function rk4(u, R, du)
        k1 = dR_du(u, R)
        k2 = dR_du(u + 0.5*du, R+0.5*du*k1)
        k3 = dR_du(u + 0.5*du, R + 0.5*du*k2)
        k4 = dR_du(u + du, R+ du*k3)
        return R + (du/6) * (k1 + 2*k2 + 2*k3 + k4)
    end

    @inbounds for i in 2:1100001
        @inbounds for j in 1:1100001
            R_grid[j, i] = rk4(u[i-1], R_grid[j, i-1], du)
        end
    end
    H_grid = Array{Float64}(undef, 1100001, 1100001)
    H_grid[:,1] = H_initial
    q_grid = Array{Float64}(undef, 1100001,1100001)
    q_grid[:,1] = q_inital
    f_grid = Array{Float64}(undef,1100001,1100001)
    f_grid[:,1] = f_initial

    @inbounds for i in 2:1100001
        @inbounds for j in 1:1100001
            q_grid[j,i] = q_grid[j,i-1] + f_grid[j,i-1]*du
            H_grid[j,i] = H_grid[j,i-1] + q_grid[j,i-1]*(R_grid[j,i]- R_grid[j,i-1])
            f_grid[j,i] = q_grid[j,i]*dg_dR(u[i],R_grid[j,i]) + H_grid[j,i]*((0.5*dg2_dR2(u[i], R_grid[j, i]))-1)
        end
        if u[i] in u_slice
            index = findfirst(x -> x == u[i], u)
            R_slice = R_grid[index,index]
            H_slice = H_slice[:,index]
            
        end
    end


end     
