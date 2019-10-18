using Plots
using StatsPlots
using DifferentialEquations
using ParameterizedFunctions
using IterableTables
using JLD2

include("../append_solution.jl")
include("../export_data.jl")
pyplot()

#Based on dede et al 2018
const Q_liv = 1084.3 #Blood plasma flow to liver (L/day)
const Q_kid = 737.3 #Blood plasma flow to kidney (L/day)
const Q_bone = 130.1
const Q_wp = 86.74 #Blood plasma flow to well-perfused kidney (L/day)
const Q_pp = 346.98 #Blood plasma flow to poorly-perfused kidney

const V_liv = 1.75 #Volume of liver (L)
const V_kid = 0.29 #kidney (L)
const V_bone = 9.8
const V_blood = 5.11
const V_plas = 2.81 #plasma
const V_wp = 1.96
const V_pp = 58.24

const P_liv = 100 #Partion coefficient, liver/plasma (no unit)
const P_kid = 100
const P_wp = 100
const P_pp = 20
const P_bone = 1000

#transfer rate from organs to blood-plasma (per day)
const k_LiB = Q_liv/(V_liv*P_liv)
const k_KB = Q_kid/(V_kid*P_kid)
const k_BoB = Q_bone/(V_bone*P_bone)
const k_WB = Q_wp/(V_wp*P_wp)
const k_PB = Q_pp/(V_pp*P_pp)

#transfer rate from blood-plasma to organs (per day)
const k_BLi = Q_liv/V_plas
const k_BK = Q_kid/V_plas
const k_BW = Q_wp/V_plas
const k_BBo = Q_bone/V_plas
const k_BP = Q_pp/V_plas

const A_gi = 0.06 #PB absorption coefficient from GI tract, ranges from 0.06 to 0.12 (unitless)
const eU = 0.47 #(/day)
const eB = 0.2 #(/day)

# Kinetic rates from parameter optimization
const k_GB = 0.44429068804813565 #(/day)
const k_GE = 7.027166796997966 #(/day)

const HCT = 0.45 # haematocrit fraction of whole blood
const R = 1.2 #ratio of unbound erythrocyte PB concentration to plasma PB concentration
const BIND = 0.437 #Pb binding capacity of erythrocytes (mg Pb L⁻¹ cell)
const KBIND = 3.72e-4 #Binding constant of erythrocytes (mg Pb L⁻¹ cell)

Av = 6.022e23
V_cell = 8e-16
V_gut = 1.2

mcg2ions(mcg, molarmass) = mcg*1e-6/molarmass*Av
mcg2mcm(mcg,molarmass) = mcg/molarmass
ions2mcg(ions, molarmass) = ions/Av*molarmass*1e6
mcm2ions(mcm) = mcm*1e-6*Av

k_growth = log(2)/40*24        # Growth rate [/day], double every 40 hours in human gut.
N_max = 1e13                 # Max number of cells
k_tsl = 0.075*1e-9*Av*V_cell*60*24   # Translation rate [/day]
tsc =  15*60*24                       # Transcription rate [/day]
k_rna_deg = log(2)/5*60*24                # RNA degradation rate [/day]
k_pro_deg = log(2)/(2/3)*24                   # Protein degradation rate [/day]
#V_max = 3.5676e-1*24                  # unknown [ion/day/TP] (there's no difference for this V_max but there is when V_max is increased by say 50 times)
V_max = 40.1*24   # from fitting [ion/day/TP] (there's no difference for this V_max but there is when V_max is increased by say 50 times)
#k_m = 3.9*1e-4              # unknown [ion/m^3]
k_m = 1e4              # unknown [ion/m^3]
k_CA = 0.1                        # unknown Accumulation protein capture rate
k_CG = 0.01                       # unknonw Ions Cytoplasm --> Gut
N0 = 1e11                       # Number of cells at start

lead_coupled = @ode_def begin #ODE system for constant intake of lead per day over several days
    dG =  -k_GB * G - k_GE * G - V_max/(1+k_m/(G/V_gut))*TP*N
    dLi = k_GB * G - (k_LiB + eB) * Li + k_BLi * B #amount in liver
    dK = - (k_KB + eU) * K + k_BK * B #amount in kidney
    dBo = - k_BoB * Bo + k_BBo * B #amount in bone
    dWP = - k_WB * WP + k_BW * B #amount in well-perfused tissues
    dPP = - k_PB * PP + k_BP * B #amount in poorly-perfused tissues
    dB = (- (k_BLi + k_BK + k_BBo + k_BW + k_BP) * B #amount in blood plasma
     + k_LiB * Li + k_KB * K + k_BoB * Bo + k_WB * WP + k_PB * PP)
    dU = eU * K #urine
    dBi = eB * Li #Biliary excretion
    dF = k_GE * G #Faeces

    # Colony
    dN = k_growth * N - k_growth * (N^2) / N_max    # Number of bacteria in gut

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
    dC = (V_max/(1+k_m/(G/V_gut))*TP # Number of ions in cytoplasm
        - k_CA * C * AP_free
        + k_pro_deg * AP_occupied
        - dN/N * C)
    dmRNA = tsc - k_rna_deg * mRNA
end a

lead_once = @ode_def begin #ODE system for constant intake of lead per day over several days
    dG =  -k_GB * G - k_GE * G  #IR_gi = oral intake rate  of Pb (mg/day) (given in solving function) #amount in gut
    dLi = k_GB * G - (k_LiB + eB) * Li + k_BLi * B #amount in liver
    dK = - (k_KB + eU) * K + k_BK * B #amount in kidney
    dBo = - k_BoB * Bo + k_BBo * B #amount in bone
    dWP = - k_WB * WP + k_BW * B #amount in well-perfused tissues
    dPP = - k_PB * PP + k_BP * B #amount in poorly-perfused tissues
    dB = (- (k_BLi + k_BK + k_BBo + k_BW + k_BP) * B #amount in blood plasma
     + k_LiB * Li + k_KB * K + k_BoB * Bo + k_WB * WP + k_PB * PP)
    dU = eU * K #urine
    dBi = eB * Li #Biliary excretion
    dF = k_GE * G #Faeces
end a

function CB(CPLASMA) #Function to calculate whole blood concentration from plasma concentration
    CB = (((1 - HCT) * CPLASMA)
        + (HCT * CPLASMA * (R + (BIND/(KBIND + CPLASMA)))))
    return CB
end

function wholeblood_conc(solution) #take an array of plasma concentration and turn it into an array of wholeblood concentration (μg/dL)
    wholeblood = []
    for i in solution[7,:]
        push!(wholeblood,CB(ions2mcg(i,207.2)/(V_plas*10)) * 1000) #(μg dL⁻¹)
    end
    return wholeblood
end

function pb_body_uptake(days, intake, init_bac, coupled = false, u0=zeros(16))
    tspan = (0.0,1.0)
    if coupled
        model = lead_coupled
        #u0 = zeros(16)
        u0[1] = intake
        u0[11] = init_bac
    else
        model = lead_once
        u0 = zeros(10)
        u0[1] = intake
    end
    sol = solve(ODEProblem(model, u0, tspan,1), dense=false)
    for t in 2:days
        tspan = tspan .+ 1
        u0 = sol.u[end]
        u0[1] += intake #
        if coupled u0[11] = init_bac end
        sol1 = sol
        sol2 = solve(ODEProblem(model, u0, tspan,1), dense=false)
        sol = sol_append(sol1, sol2)
    end
    sol
end

function wrapping()
    u2 = zeros(16)
    for i in 1:8*4
        global sol2 = pb_body_uptake(91,G0,N0,true, u2)

        u2 = sol2.u[end]
        #export_to_csv(sol2, "pb_11mcg_per_day_8y_PRODEACC$i", 2000)
    end
    #@save "pb_11mcg_per_day_final_coupled.sol" sol2
    return sol2
end

plotlyjs()
G0 = mcg2ions(114,207.2)            # Gut start amount
sol1 = pb_body_uptake(365*8,G0,N0,false)
sol2 = wrapping()
sol1_end = sol1[end]
sol2_end = sol2[end]
print(sol2_end[1:10]./sol1_end)
mn = [sol1_end[2:7], sol2_end[2:7]]

data = []
append!(data, mn[1])
append!(data, mn[2])
data = float.(data)
g = repeat(["Without Probiotic", "With Probiotic"], inner = 6)
nam = [
    "Liver", "Kidney", "Bone", "WP Tissue", "PP-tissue", "Blood",
    "Liver", "Kidney", "Bone", "WP Tissue", "PP-tissue", "Blood",
    ]
groupedbar(nam, data,  group = g, ylabel="Number of Pb(II) ions")
savefig("pb_compare_8y_114mcg.svg")
savefig("pb_compare_8y_114mcg.png")
#export_to_after_effects(sol, "pb_114mcg_per_day_8y_PRODEACC1", 2000)
#average daily lead intake through diet is about 114 μg/day
