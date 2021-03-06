# This models the uptake of ions for a colony of E. coli inside the GI-trackt

using DifferentialEquations
using ODEInterfaceDiffEq
using ParameterizedFunctions
using Plots
include("../append_solution.jl")

pyplot()

# Constants
const Av = 6.022e23           # Avogrado's number
const V_cell = 8e-16          # Cell volume   [l]
const V_gut = 1.2             # Gut volume    [l]
ions2mcg(ions, molarmass) = ions/Av*molarmass*1e6
mcg2ions(mcg, molarmass) = mcg*1e-6/molarmass*Av

bacteria_uptake = @ode_def begin
    # Colony
    dN = k_growth * N - k_growth * (N^2) / N_max    # Number of bacteria in gut
    dG =-V_max/(1+k_m/(G/V_gut))*TP*N                  # Number of ions in gut

    # Single Bacteria
    dTP = (k_tsl * mRNA
        - dN/N * TP
        - k_pro_deg * TP)
    dAP_free = (k_tsl * mRNA           # Number of free accumulators
        - dN/N * AP_free
        - k_CA * AP_free * C
        - k_pro_deg * AP_free)
    dAP_occupied = (k_CA * AP_free * C  # Number of occupied accumulators
        - dN/N * AP_occupied
        - k_pro_deg * AP_occupied)
    dC = (V_max/(1+k_m/(G/V_gut))*TP       # Number of ions in cytoplasm
        - k_CA * C * AP_free
        + k_pro_deg * AP_occupied
        - dN/N * C)
    dmRNA = tsc - k_rna_deg * mRNA
end k_growth N_max k_tsl tsc k_rna_deg k_pro_deg V_max k_m k_CA k_CG N0 G0 c lag ind TP0 AP0

function colony_model(t, symbols::Array{}=[])
    value(key::Symbol) = p[findfirst(x -> x[1] == key, p)][2]
    index(key::Symbol) = findfirst(x -> x[1] == key, p)

    # Parameters
    p = (
    [
    (:k_growth, log(2)/0.5)         # Growth rate, double every 30 min. [/h]
    (:N_max, 1e13)                      # Max number of cells
    (:k_tsl, 0.075*1e-9*Av*V_cell*60)   # Translation rate [/h]
    (:tsc, 15*60)                       # Transcription rate [/h]
    (:k_rna_deg, 20.79)                 # RNA degradation rate [/h]
    (:k_pro_deg, 0.1)                   # Protein degradation rate [/h]
    (:V_max, 3.5676e-1)                 # unknown [ion/h/TP]
    (:k_m, 2.3458e14)              # unknown [ion/m^3]
    (:k_CA, 0.1)                        # unknown Accumulation protein capture rate
    (:k_CG, 0.01)                       # unknonw Ions Cytoplasm --> Gut
    (:N0, 1e9)                          # Number of cells at start
    (:G0, mcg2ions(0.25, 75)*V_gut)     # Gut start amount
    (:lag_phase, 0.0)                   # Lag phase duration [h]
    (:induced_time, 0.0)
    (:TP0, 0.0)
    (:AP0, 0.0)
    ])              # Induction time [h]

    # Replace values
    for s in symbols p[index(s[1])] = (s[1], s[2]) end
    k = map(x -> x[2], p)

    # Set initial condidionsge
    u0_1 = zeros(7)
    u0_1[1] = value(:N0)
    u0_1[2] = value(:G0)
    u0_1[3] = value(:TP0)
    u0_1[4] = value(:AP0)
    u0_1[5] = 1e-10
    u0_1[6] = 1e-10

    t_lag = value(:lag_phase)
    t_ind = value(:induced_time)

    # Solve before induced_time if induced_time =! 0
    if t_ind != 0
        k[4]=0 # Transcription rate
        tspan1 = (0.0, t_ind-t_lag)
    else
        tspan1 = (0.0, t)
    end
    prob1 = ODEProblem(bacteria_uptake,u0_1,tspan1,k)
    sol1 = solve(prob1, isoutofdomain=(u,p,t) -> any(x -> x < 0, u), dense=false)

    # Ajust for lag phase
    sol1.t .+= t_lag
    pushfirst!(sol1.u, u0_1)
    pushfirst!(sol1.u, u0_1)
    pushfirst!(sol1.t, t_lag-0.01)
    pushfirst!(sol1.t, 0.0)
    pushfirst!(sol1.alg_choice, 1)
    pushfirst!(sol1.alg_choice, 1)
    if t_ind == 0
        return sol1
    end

    # Solve after induced_time
    k[4] = value(:tsc)
    tspan2 = (0, t-t_ind)
    u0_2 = sol1.u[end]
    prob2 = ODEProblem(bacteria_uptake,u0_2,tspan2,k)
    sol2 = solve(prob2,Rosenbrock23(autodiff=false),maxiters=3e5, isoutofdomain=(u,p,t) -> any(x -> x < 0, u), dense=false)
    sol2.t .+= t_ind

    # Append solutions
    sol_append(sol1, sol2)
end

function plot_init()
    plot(layout=(3,2))
end

function plot_model(sol)
    p1 = plot(sol, vars=[1],
        label="Number of cells",
        color="Red")
    p2 = plot(sol, vars=[2],
        label="Ions in medium",
        color="Blue",
        yscale=:identity)
    p3 = plot(sol,
        vars=[3,4,5],
        label=["Transporters" "Free accumulators" "Occupied accumulators"],
        yscale=:identity)
    sum_of_ions = sol[2,:].+sol[1,:].*(sol[5,:].+sol[6,:])
    control = plot(sol.t,sum_of_ions,
        label="Total ions (shoud be constant)",
        color="Black",
        linewidth=2,
        ylims=(0.99*sol[2,1], 1.01sol[2,1]))
    plot(p1, p2, p3, control, layout=(2,2))
end

sol = colony_model(5.0, [])

plot_init()
plot(sol, vars=[1])
plot_model(sol)
