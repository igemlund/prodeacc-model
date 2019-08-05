# This is for modelling the Pb(II) accumulation in an E. coli cell using
# different proteins

using DifferentialEquations
using ParameterizedFunctions
using Plots
plotly()
plot()
Av = 6.022e23                       # Avogrado's number
V_cell = 8e-16                      # Cell volume   [l]

# Transcription Rates
tsc = 5.85              # [mRNA/min]

# Translation Rates ()
tlr_pbrR = 3.32         # [/min]
tlr_pbrD = 15.88        # [/min]
tlr_pbrT = 5.98         # [/min]

# Gene Lengths
pbrD_length = 726       # [bp]
pbrR_length = 438       # bp
pbrT_length = 1929      # bp

# Degradation
k_mrna_deg = 0.0693     # [/min]
k_pro_deg = 0.0333      #
5/k_mrna_deg
# Affinity
# aff_pbrR = 40e-6*1000*V_cell*207.2
aff_pbrD = 9.86e-5      # Fraction of removed pb per pbrR

pbrD_uptake = @ode_def begin
    dmRNA_pbrD = tsc - k_mrna_deg * mRNA_pbrD
    dmRNA_pbrT = tsc  - k_mrna_deg * mRNA_pbrT
    dpbrD = tlr_pbrD * mRNA_pbrD - k_pro_deg
    dpbrT = tlr_pbrT * dmRNA_pbrT - k_pro_deg
end a

pbrR_uptake = @ode_def begin
    dmRNA_pbrR = tsc - k_mrna_deg * mRNA_pbrR
    dmRNA_pbrT = tsc  - k_mrna_deg * mRNA_pbrT
    dpbrR = tlr_pbrR * mRNA_pbrR - k_pro_deg
    dpbrT = tlr_pbrT * dmRNA_pbrT - k_pro_deg
end a

u0 = zeros(4)
tspan = (0.0, 4000)
prob = ODEProblem(pbrD_uptake,u0,tspan,1)
pbrD_sol = solve(prob)

prob = ODEProblem(pbrR_uptake,u0,tspan,1)
pbrR_sol = solve(prob)

plotlyjs()
plot()
plot!(pbrD_sol, vars=[3])
plot!(pbrR_sol, vars=[3])
