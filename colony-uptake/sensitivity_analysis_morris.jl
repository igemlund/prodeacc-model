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
timestep = collect(range(0.0, stop = t, length = 200))
morris = morris_sensitivity(prob,Tsit5(),timestep,[[0,20],[0,20],[0,10],[0,10]],[200,200,100,100])

stdv1 = sqrt.(morris.variances[1])
stdv2 = sqrt.(morris.variances[2])

plot(morris.means[1]', layout = (3,2), label = "V_max")
plot!(morris.means[2]', label = "k_m")
plot!(morris.means[3]', label = "k_CA")
plot!(morris.means[2]', label = "k_CG")
