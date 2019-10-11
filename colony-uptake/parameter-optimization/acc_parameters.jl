using CSV
using Optim
include("../colony-uptake.jl")
include("../optimize_model.jl")

c = [6,	0.32].*1e-3/74.9*Av # Pb Tac 200M pH3
t = [0, 2.5]

f(x)=loss_function(target_data(t, 2, [(:G0, c[1]), (:N0, 1e11), (:k_growth, 0.0), (:V_max, x[1]), (:k_m, x[2]), (:k_CA, x[3]), (:TP0, 5000), (:AP0, 5000)]), c)
options = Optim.Options(show_trace=true,
                        show_every=100,
                        iterations=8000,
                        extended_trace=true,
                        store_trace=true,
                        time_limit=1000)

# k_growth, N_max
k_upper =   [1e6, 1e12, 1.0]
k_init =    [1e2, 1e8, 0.1]
k_lower =   [1e1, 1e6, 1e-2]


#alg = ParticleSwarm(n_particles=6, lower=k_lower, upper=k_upper)
alg = NelderMead()
time = @elapsed res = optimize(f, k_init,alg, options)
V_max, k_m, k_CA = res.minimizer
sol = colony_model(12, [(:G0, c[1]),(:N0, 1e11),(:k_growth, 0.0), (:lag_phase, 0.0),(:V_max, V_max), (:k_m, k_m), (:k_CA, k_CA)])
plotlyjs()
plot(sol, vars=[2],
    label="CFU according to model",
    linewidth=4,
    #linecolor= 'b',
    ylabel="CFU",
    xlabel="hours")
scatter!(t ,c, label="CFU measured in wetlab")
