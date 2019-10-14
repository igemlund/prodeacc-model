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
    (:tsc, 15*60)                       # Transcription rate [/h]
    (:k_rna_deg, 20.79)                 # RNA degradation rate [/h]
    (:k_pro_deg, 0.1)                   # Protein degradation rate [/h]
    (:V_max, 3.5676e-1)                 # unknown [ion/h/TP]
    (:k_m, 3.9e-6*Av*1e-4)              # unknown [ion/m^3]
    (:k_CA, 0.1)                        # unknown Accumulation protein capture rate
    (:k_CG, 0.01)                       # unknonw Ions Cytoplasm --> Gut
    (:N0, 1e9)                          # Number of cells at start
    (:G0, mcg2ions(0.25, 75)*V_gut)     # Gut start amount
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

factor = 20
t = 3.0
parameter_range(x) = [x/factor, x*factor]
pr = parameter_range.(k)
tspan = (0.0, t)
prob = ODEProblem(bacteria_uptake,u0_1,tspan,k)
timestep = collect(range(0.0, stop = t, length = 1000))
sobol = DiffEqSensitivity.Sobol(N=4000, order=[0,1])
s1 = DiffEqSensitivity.gsa(prob,Tsit5(),timestep,pr,sobol)
s2 = s1.S1

plotlyjs()
i = 7
p1 = bar(["k_tsl", "tsc", "k_rna\n_deg", "k_pro\n_deg", "V_max", "k_m", "k_CA", "k_CG", "N0", "G0", "lag_phase"],
        [
    #     s2[1][i,end-1],
         s2[2][i,end-1],
         s2[3][i,end-1],
         s2[4][i,end-1],
         s2[5][i,end-1],
         s2[6][i,end-1],
         s2[7][i,end-1],
         s2[8][i,end-1],
         s2[9][i,end-1],
        s2[10][i,end-1],
        s2[11][i,end-1],
        s2[12][i,end-1]])
savefig("sobo.png")
plot(p1)
