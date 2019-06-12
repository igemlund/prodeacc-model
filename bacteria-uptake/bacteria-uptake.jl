using DifferentialEquations
using ParameterizedFunctions
using Plots
using PyPlot
using PyCall
# http://homepages.ulb.ac.be/~dgonze/BIONUMBERS/bionumbers.html
V_cell = 6e-16                      # Cell volume   [l]
V_gut = 1.2                         # Gut volume    [l]
Av = 6.022e23                       # Avogrado's number
I_tsc_ = 0.02                       # Transcription rate [nM/h]
I_tsc = I_tsc_ * V_cell * 1e-9 * Av # Transcription rate [/h]
k_tsl = 0.1 * V_cell * 1e-9 * Av    # Translation rate [/h]
k_rna_deg = 50                      # RNA degradation rate [/h]
k_pro_deg = 0.2                     # Protein degradation rate [/h]

G0 = 0.25e-3 * 75 * Av      # [amount] gut amount
V_max = 500e-6 * Av         # [mol/h]    unknown
k_transport = 15000 * Av    # [mol/m^3]  unknown
k_CG = 1.0                  # unkown
k_CA = 1.0                  # unkonw

bacteria_uptake = @ode_def begin
    dmRNAᵀ = I_tsc - k_rna_deg * mRNAᵀ      # [amount] Transport mRNA
    dmRNAᴬ = I_tsc - k_rna_deg * mRNAᴬ      # [amount] Accumulation mRNA
    dTP = mRNAᵀ * k_tsl - k_pro_deg * TP    # [amount] Transport protein
    dAP_free = mRNAᴬ * k_tsl - k_pro_deg * AP_free - k_CA * AP_free * C    # [amount] Accumulation protein
    dAP_occupied = k_CA * AP_free * C - k_pro_deg * AP_occupied # Accumulated HM
    dG = -V_max / (1 + k_transport/(G/V_gut - C/V_cell))*TP + k_CG * C     # [amount] Gut concentration
    dC =  V_max / (1 + k_transport/(G/V_gut - C/V_cell))*TP - k_CA * AP_free * C + k_pro_deg * AP_occupied
      # [amount] Cytoplasm concentration
end a

u0 = [0;0;0;0;G0;0;0]
tspan = (0.0, 12.0)
prob = ODEProblem(bacteria_uptake,u0,tspan,1)
sol = solve(prob)
pyplot()
plot(sol)

# plot!(sol.t, sol[2,:])
# plot!(sol.t, sol[3,:])
# plot!(sol.t, sol[4,:])
# plot!(sol.t, sol[5,:])
# plot!(sol.t, sol[6,:])
# plot!(sol.t, sol[7,:])
