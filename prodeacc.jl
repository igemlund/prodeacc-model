using Plots
using DifferentialEquations
using ParameterizedFunctions

k_cp = 0.916    # Rate central --> periphiral compartement
k_pc = 0.0307   # Rate periphiral --> central compartement
k_ei = 1.45     # Rate clearence from central compartement
a_ka = 0.118
b_ka = 0.698
G    = 0.1

g = @ode_def begin
    dP = k_cp * C - k_pc * P
    dC = k_pc * P - k_cp* C - k_ei * C + a_ka * G
    #dG = -a_ka * G + 1
end k_cp k_pc k_ei a_ka G

u0 = [0.0;0.0;0]
tspan = (0.0, 10)
p = [k_cp, k_pc, k_ei, a_ka, G]
prob = ODEProblem(g,u0,tspan,p)
sol = solve(prob)[2:end]
plot(sol)
