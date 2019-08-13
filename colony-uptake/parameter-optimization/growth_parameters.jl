using CSV
using Optim

EcN = CSV.read("data/EcN.csv")
EcNT = CSV.read("data/EcNT.csv")
EcN_od = [EcN[:,1], EcN[:,5]]
EcNT_od = [EcNT[:,1], EcNT[:,5]]
deleteat!(EcN_od[1], 10)
deleteat!(EcN_od[2], 10)
scatter(EcN_od[1], EcN_od[2].*1e9)

f(x)=loss_function(target_data(float.(EcN_od[1]), 1, [(:k_growth, x[1]), (:N_max, x[2])]), EcN_od[2].*1e9)
options = Optim.Options(show_trace=true,
                        show_every=100,
                        iterations=8000,
                        extended_trace=true,
                        store_trace=true,
                        time_limit=1000)

# k_growth, N_max
upper = [log(2)/0.1, 3e10]
lower = [log(2)/2, 1e10]
k_init = [log(2)/1, 1e10]

alg = ParticleSwarm(n_particles=6, lower=lower, upper=upper)
alg = NelderMead()

time = @elapsed res = optimize(f, k_init,alg, options)
k_growth, N_max = res.minimizer
plotlyjs()
plot!(colony_model(12, [(:k_growth, k_growth), (:N_max, N_max)]), vars=[1])
