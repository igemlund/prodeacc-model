# This models the uptake of ions for a colony of E. coli inside the GI-trackt

using DifferentialEquations
using ODEInterfaceDiffEq
using ParameterizedFunctions
using Plots
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
k_pro_deg = 0.1          # Effective Protein degradation rate [/h]
V_max = 500e-6 * Av                 # unknown [/h]
k_m = 0.1 * Av / V_cell     # unknown [mol/m^3]
k_CA = 1.0                          # unknown Accumulation protein capture rate
k_CG = 1.0                          # unknonw Ions Cytoplasm --> Gut

# Start conditions
G0 = mcg2ions(0.25, 75) * V_gut     # Gut start amount
N0 = 1e9                            # Bacteria start amount
N_max = 1e13                        # Assumed max amount of prodeacc in human gut

# TODO: motvate this
mRNAᵀ = 4
mRNAᴬ = 4
V_max = V_max * 1e-5
bacteria_uptake = @ode_def begin
    # Colony
    dN = k_growth * N - k_growth * (N^2) / N_max                                    # Number of bacteria in gut
    dG =-V_max/(1+k_m/((G/V_gut-C/V_cell)/thickness))*N    # Number of ions in gut

    # Single Bacteria
    dAP_free = (k_tsl * mRNAᴬ           # Number of free accumulators
        - dN / N * AP_free
        - k_CA * AP_free * C
        - k_pro_deg * AP_free)
    dAP_occupied = (k_CA * AP_free * C  # Number of occupied accumulators
        - dN / N * AP_occupied
        - k_pro_deg * AP_occupied)
    dC = (V_max/(1+k_m/((V_gut-C/V_cell)/thickness))    # Number of ions in cytoplasm
        - k_CA * C * AP_free
        + k_pro_deg * AP_occupied
        - dN / N * C)
end V_max k_m

function colony_uptake(t, K, c)
    u0 = zeros(5)
    u0[1] = N0
    u0[2] = mcg2ions(c, 75)
    tspan = (0.0, t)
    prob = ODEProblem(bacteria_uptake,u0,tspan,K)
    sol = solve(prob)
end

function plot_init()
    plot(layout=(3,2))
end

function plot_model(sol)
    p1 = plot(sol, vars=[1], label="Number of cells", color="Red")
    p2 = plot(sol, vars=[2], label="Ions in medium", color="Blue")
    p3 = plot(sol, vars=[3,4,5], label=["Free accumulators" "Occupied accumulators" "Ions in cytoplasm"])
    sum_of_ions = sol[2,:].+sol[1,:].*(sol[4,:].+sol[5,:])
    control = plot(sol.t, sum_of_ions, label="Total ions (shoud be constant)", color="Black")
    plot(p1, p2, p3, control, layout=(2,2))
end

function plot_transport(G)
    flux(C) = V_max/(1+k_m/C)
    plot(flux.(G))
end

K = [V_max, k_m]
plotlyjs()
plot()
plot_init()
sol = colony_uptake(50.0, K, 0.25)
plot_model(sol)
#plot_transport(colony_uptake(18.0, K, 0.25)[2,:])
