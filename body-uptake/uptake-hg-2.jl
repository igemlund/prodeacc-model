# Based on Carrier et al 2000
using Plots
using DifferentialEquations
using ParameterizedFunctions

K = 12.9870         #Constant ratio Q⬚(t)/B⬚(t)
# ---- Rates [days^-1] ----
# Organic mercury
k_abs = 5.5440     # Oral absorption rate constant
k_QI  = 0.01347      # Metabolism rate constant of organic mercury to inorganic mercury
k_QF  = 9.0668e-5   # Whole body to feces transfer coefficient of organic mercury
k_QU  = 0           # Whole body to urine transfer coefficient of organic mercury
k_QH  = 2.3825e-4   # Whole body to hair transfer coefficient of organic mercury
k_elim = k_QU + k_QF + k_QI + k_QH  # Whole body elimination rate constant of
                                    # organic mercury

# Incorganic mercury
d_BL = 0.1750       # Blood to liver transfer coefficient combined with liver
                    # metabolism rate constant of organic mercury
d_BBr = d_BL/1000   # Blood to brain transfer coefficient combined with brain
                    # metabolism rate constant of organic mercury
k_LB = 0.8940       # Liver to blood transfer coefficient of inorganic mercury
k_BK = 17.1234      # Blood to kidney transfer coefficient of inorganic mercury
k_KB = 0.0010       # Kidney to blood transfer coefficient of inorganic mercury
k_KU = 0.006949     # Kidney to urine transfer coefficient of inorganic mercury
k_BH = 0.1400       # Blood to hair transfer coefficient of inorganic mercury
k_BU = 0.06994      # Blood to urine transfer coefficient of inorganic mercury
k_BF = 3.9917       # Blood to feces transfer coefficient of inorganic mercury
k_LF = 1.5476       # Liver to feces transfer coefficient of inorganic mercury
k_BBr = 0.0028      # Blood to brain transfer coefficient of inorganic mercury
k_BrB = 0.0520      # Brain to blood transfer coefficient of inorganic mercury
k_BL = 10           # TODO: Find the value for this. Was not included in the paper.

organic_simplified = @ode_def begin
    dg = - k_abs * g
    dQ = k_abs * g - (k_QI + k_QH + k_QF + k_QU) * Q
    dE = (k_QI + k_QH + k_QF + k_QU) * Q
end k_abs k_QI k_QH k_QF k_QU

inorganic_simplied = @ode_def begin
    dBₒ = - (d_BBr + d_BL) * Bₒ
    dBr = - k_BrB * Br + k_BBr * B
    dB  = - (k_BBr + k_BK + k_BL + k_BH + k_BF + k_BU) * B + k_KB * K + k_LB * L + k_BrB * Br
    dL  = - (k_LB + k_LF) * L + k_BL * B + d_BL * Bₒ
    dK  = - (k_KB + k_KU) * K + k_BK * B
    dU  = k_BU * B + k_KU * K
    dF  = k_BF * B + k_LF * L
    dH  = k_BH * B
end d_BBr d_BL k_BrB k_BBr k_BK k_BH k_BU k_KB k_LB k_LF k_BL k_BH k_BF k_KU

#carrier_II = @ode_def begin
#    dI = (k_QI * K) * Bₒ
#end

o_0 = [1.0;0.0;0]
o_p = [k_abs, k_QI, k_QH, k_QF, k_QU]
io_0 = [1.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0]
io_p = [d_BBr, d_BL, k_BrB, k_BBr, k_BK, k_BH, k_BU, k_KB, k_LB, k_LF, k_BL, k_BH, k_BF, k_KU]

tspan = (0.0, 365)
prob = ODEProblem(inorganic_simplied,io_0,tspan,io_p)
sol = solve(prob)[2:end]

# Plotting
plotly()
plot()                  
plot!(sol, vars=(1))
plot!(sol, vars=(2))
plot!(sol, vars=(3), yscale=:log)
plot!(sol, vars=(4))
plot!(sol, vars=(5))
plot!(sol, vars=(6))
plot!(sol, vars=(7))
plot!(sol, vars=(8))
