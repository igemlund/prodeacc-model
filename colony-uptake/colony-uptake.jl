# This models the uptake of ions for a colony of E. coli inside the GI-trackt

using DifferentialEquations
using ODEInterfaceDiffEq
using ParameterizedFunctions
using Plots
pyplot()
# http://book.bionumbers.org/how-many-proteins-are-made-per-mrna-molecule/
# Constants
Av = 6.022e23                       # Avogrado's number
V_cell = 8e-16                      # Cell volume   [l]
V_gut = 1.2                         # Gut volume    [l]
thickness = 10e-9                    # Membrane thickness [meters]
ions2mcg(ions, molarmass) = ions/Av*molarmass*1e6
mcg2ions(mcg, molarmass) = mcg*1e-6/molarmass*Av

# Rates
k_tsl = 1*60*60                     # Translation rate [/h]
k_rna_deg = 20.79                   # RNA degradation rate [/h]
k_growth = log(2)/1                 # Growth rate, double every 30 min. [/h]
k_pro_deg = 1/(1000/60/60)          # Effective Protein degradation rate [/h]
k_pro_deg = k_growth          # Effective Protein degradation rate [/h]
V_max = 500e-6 * Av                 # unknown [/h]
k_transport = 1.5 * Av / V_cell     # unknown [mol/m^3]
k_CA = 1e3                          # unknown Accumulation protein capture rate
k_CG = 1.0                          # unknonw Ions Cytoplasm --> Gut

# Start conditions
G0 = mcg2ions(0.25, 75) * V_gut     # Gut start amount
N0 = 1e9                            # Bacteria start amount

# TODO: motvate this
mRNAᵀ = 4
mRNAᴬ = 4

bacteria_uptake = @ode_def begin
    # Colony
    dN = k_growth * N                                      # Number of bacteria in gut
    dG =(-V_max/(1+k_transport/((G/V_gut-C/V_cell)/thickness))*N*TP    # Number of ions in gut
        + k_CG * N * C)

    # Single Bacteria
    dTP = (k_tsl * mRNAᵀ                                 # Number of transporters
        - k_pro_deg * TP)
    dAP_free = (k_tsl * mRNAᴬ                            # Number of free accumulators
        - k_pro_deg * AP_free
        - k_CA * AP_free * C)
    dAP_occupied = (k_CA * AP_free * C                   # Number of occupied accumulators
        - k_pro_deg * AP_occupied)
    dC = (V_max/(1+k_transport/((G/V_gut-C/V_cell)/thickness))*TP    # Number of ions in cytoplasm
        - k_CG * C
        - k_CA * C * AP_free
        - k_pro_deg * C)
end a

u0 = zeros(6)
u0[1] = N0
u0[2] = G0
tspan = (0.0, 20.0)
prob = ODEProblem(bacteria_uptake,u0,tspan,1)
sol = solve(prob)

# Relative change of ions
ion_sum(sol) = (sol[2,:] + (sol[5,:]+sol[6,:]).*sol[1,:])./G0
plot(sol.t, ion_sum(sol), label="total")#, ylim=[0.0, 2.0])
plot!(sol.t, sol[2,:]/G0, label="in GI")
plot!(title="Relative change")
#plot(sol, vars=[5 6])