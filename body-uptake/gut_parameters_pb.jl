# This script is used for determining the kinetic rates for Gut-> Blood and
# for Gut -> Faeces
using Optim
include("uptake-pb.jl")

function loss_function(k)
    global k_GB = k[1]
    global k_GE = k[2]
    days = 40
    dose = 105e-3
    t = 1:2:days
    try
        sol = pb_body_uptake(days, dose, false)(t)
        sol_mod = pb_body_uptake(days, dose, true)(t)
        return sum(abs.(sol .- sol_mod))/length(sol)
    catch e
        @warn "Did not converge."
        @warn e
        return Inf
    end
end

l = [0.0, 0.0]
u = [10.0, 20.0]
timelimit = 1000
iterations = 1000
options = Optim.Options(show_trace=true,
                        show_every=100,
                        iterations=iterations,
                        extended_trace=true,
                        store_trace=true,
                        time_limit=timelimit)
k_init = [-log(A_gi), -log(A_gi)]
alg = ParticleSwarm(lower=l, upper=u,n_particles=6)
time = @elapsed res = optimize(loss_function, k_init, alg, options)
k_GB, k_GE = res.minimizer
display(k_GE)

sol1 = pb_body_uptake(400, 205e-2, false)
sol2 = pb_body_uptake(400, 205e-2, true)
plot(sol1, vars=[2,3])
plot!(sol2, vars=[2,3], linestyle=:dash)
plot(sol2, vars=[1], xlims=[0, 1])
