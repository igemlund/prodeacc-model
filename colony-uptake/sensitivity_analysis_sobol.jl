using DifferentialEquations
using ODEInterfaceDiffEq
using ParameterizedFunctions
using Plots
using DiffEqSensitivity
include("colony-uptake.jl")

K = [V_max, k_m, k_CA, k_CG]
t = 5.0
c = 0.25
u0 = zeros(6)
u0[1] = N0
u0[2] = mcg2ions(c, 75)
tspan = (0.0, t)
prob = ODEProblem(bacteria_uptake,u0,tspan,K)
timestep = collect(range(0.0, stop = t, length = 100))
s1 = sobol_sensitivity(prob,Tsit5(),timestep,[[1e10,1e20],[0,5],[0,5],[0,5]],100,1)

p1 = bar(["V_max","k_m","k_CA","k_CG"],[s1[1][5,50],s1[2][5,50],s1[3][5,50],s1[4][5,50]])
plot(p1)
