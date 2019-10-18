using CSV
using Optim
include("../colony-uptake.jl")
include("../optimize_model.jl")

c0 = sum([0.52, 1.1, 1])/3*1e-3/70.9*Av
c1 = 0.8 * c0
c = [c0, c1]
t = [0.0, 2.5]

f(x)=loss_function(target_data(t, 2,
    [(:G0, c[1]),
    (:N0, 1e13),
    (:V_max, x[1]),
    (:k_m, x[2]),
    (:k_CA, x[3]),
    (:TP0, 1000),
    (:AP0, 1000)]), c)

options = Optim.Options(show_trace=true,
                        show_every=100,
                        iterations=8000,
                        extended_trace=true,
                        store_trace=true,
                        time_limit=1000)

# k_growth, N_max
V_max, k_m, k_CA = [40.1, 1e-4, 1e-1]
k_init =    [V_max, k_m, k_CA]
k_upper = [1e3, 1e1, 3]
k_lower =  [1e1, 1e-4, 0.1]

alg = ParticleSwarm(n_particles=5, lower=k_lower, upper=k_upper)
#alg = NelderMead()
time = @elapsed res = optimize(f, k_init,alg, options)
plotlyjs()
V_max, k_m, k_CA = res.minimizer
sol = colony_model(3.0, [(:G0, c[1]),(:N0, 1e13),(:V_max, V_max), (:k_m, k_m), (:k_CA, k_CA), (:TP0, 1000), (:AP0, 1000)])
plot(sol, vars=[2],
    label="Modeled concentration",
    xlims=[0, 3.0],
    linewidth=4,
    ylabel="Number of As ions",
    xlabel="hours")
scatter!(t , c, label="Measured concentration")
savefig("fit_as.svg")
savefig("fit_as.png")
