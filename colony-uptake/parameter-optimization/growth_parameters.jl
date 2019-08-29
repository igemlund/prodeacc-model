using CSV
using Optim
include("../optimize_model.jl")


T1 = CSV.read("T1.csv")
ind = [1, 5, 8]
time = T1[1][ind]
cfu = T1[13][ind]

scatter(time, cfu)

f(x)=loss_function(target_data(float.(time), 1, [(:N0, cfu[1]), (:k_growth, x[1]), (:N_max, x[2]), (:lag_phase, x[3])]), cfu)
options = Optim.Options(show_trace=true,
                        show_every=100,
                        iterations=8000,
                        extended_trace=true,
                        store_trace=true,
                        time_limit=1000)

# k_growth, N_max
upper = [log(2)/0.3, 7e8, 6.0]
lower = [log(2)/2, 1e6, 3.0]
k_init = [log(2)/1, 2e6, 4.0]

alg = ParticleSwarm(n_particles=6, lower=lower, upper=upper)
#alg = NelderMead()
time = @elapsed res = optimize(f, k_init,alg, options)
k_growth, N_max, lag_phase = res.minimizer
plotlyjs()
sol = colony_model(12, [(:N0, cfu[1]),(:k_growth, k_growth), (:N_max, N_max), (:lag_phase, lag_phase)])
sol(0)[1]
plot!(sol, vars=[1])
