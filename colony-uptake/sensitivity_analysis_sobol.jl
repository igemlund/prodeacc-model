using DifferentialEquations
using ODEInterfaceDiffEq
using ParameterizedFunctions
using Plots
using DiffEqSensitivity
include("colony-uptake.jl")

value(key::Symbol) = p[findfirst(x -> x[1] == key, p)][2]
index(key::Symbol) = findfirst(x -> x[1] == key, p)


# Parameters
p = ([(:k_growth, log(2)/1)         # Growth rate, double every 30 min. [/h]
    (:N_max, 1e13)                      # Max number of cells
    (:k_tsl, 0.075*1e-9*Av*V_cell*60)   # Translation rate [/h]
    (:tsc, 15*60)                       # Transcription constant [amount/h]
    (:k_rna_deg, 20.79)                 # RNA degradation rate [/h]
    (:k_pro_deg, 0.1)                   # Protein degradation rate [/h]
    (:V_max, 3.5676e-1)                 # unknown [ion/h/TP]
    (:k_m, 3.9e-6*Av*1e-4)              # unknown [ion/m^3]
    (:k_CA, 0.1)                        # unknown Accumulation protein capture rate
    (:k_CG, 0.01)                       # unknonw Ions Cytoplasm --> Gut
    (:N0, 1e9)                          # Number of cells at start
    (:G0, mcg2ions(50, 75)*V_gut)     # Gut start amount
    (:lag_phase, 0.0)])                 # Lag phase duration [h]
# Replace values
k = map(x -> x[2], p)

# Set initial condidionsge
u0_1 = zeros(7)
u0_1[1] = value(:N0)
u0_1[2] = value(:G0)
u0_1[3] = 1e-10
u0_1[4] = 1e-10
u0_1[5] = 1e-10
u0_1[6] = 1e-10
c = 0.25

factor = 10
t = 3.0
parameter_range(x) = [x/factor, x*factor]
pr = parameter_range.(k)
pr[1] = [log(2)/40, log(2)/0.5] #range for k_growth
pr[2] = [1e10 , 1e14] #range for N_max
pr[4] = [120 , 1000] #range for tsc
pr[7] = [1 , 1e2] #range for V_max
pr[8] = [1 , 1e14] #range for k_m
tspan = (0.0 , t)
prob = ODEProblem(bacteria_uptake,u0_1,tspan,k)
timestep = collect(range(0.0, stop = t, length = 1000))
sobol = DiffEqSensitivity.Sobol(N=5000, order=[0,1, 2])
s1 = DiffEqSensitivity.gsa(prob,Tsit5(),timestep,pr,sobol)
s1
s2 = s1.ST

plotlyjs()
i = 2
names = ["k_growth", "N_max", "k_tsl", "tsc", "k_rna_deg", "k_pro_deg", "V_max", "k_m"]
p1 = bar(names,
        [
         s2[1][i,end-2],
         s2[2][i,end-2],
         s2[3][i,end-2],
         s2[4][i,end-2],
         s2[5][i,end-2],
         s2[6][i,end-2],
         s2[7][i,end-2],
         s2[8][i,end-2],
        ],
        legend=false)
savefig("sobol_2.svg")
plot(p1)
