using Optim
using Plots
using PrettyTables
include("colony-uptake.jl")

# Returns the number of ions in the growth medium according to the model.
function ions_in_medium(t, k, c)
    sol = colony_uptake(10.0, k, c)
    return sol(t)[2,:]
end

loss_function(data1, data2) = 1.0/length(data1)*sum(abs.(data1 - data2))

# Initial guess for the unknown constants
k_init = [V_max, k_m]

# Generated data using other constants
t = 0.0:1:5.0
k_data = [V_max*7, k_m]
init = ions_in_medium(t, k_init, 0.25)
data1 = ions_in_medium(t, k_data, 0.25)
data2 = ions_in_medium(t, k_data, 0.50)
data3 = ions_in_medium(t, k_data, 0.125)

# Cost function for optimizing model, using different concentrations
f(x) = (loss_function(data1, ions_in_medium(t,x, 0.25))
    + loss_function(data2, ions_in_medium(t,x,0.5))
    + loss_function(data3, ions_in_medium(t,x,0.125)))

# Optimizers
function opt_particle_swarm(k_init)
    options = Optim.Options(show_trace=true,
                            f_tol=1e-2,
                            extended_trace=true,
                            iterations=160)
    res = optimize(f, k_init, ParticleSwarm(lower=[1e18,1e18],
                                            upper=k_init.*7,
                                            n_particles=3),
                                            options)
    return Optim.minimizer(res)
end

function opt_auto(k_init)
    res = optimize(f, k_init)
    return Optim.minimizer(res)
end

time1 = @elapsed k_opt1 = opt_auto(k_init)
time2 = @elapsed k_opt2 = opt_particle_swarm(k_init)
time3 = @elapsed k_opt3 = opt_auto(opt_particle_swarm(k_init))

# Plot
plot_init()
plot_model(colony_uptake(15.0, k_init, 0.50))
plot_model(colony_uptake(15.0, k_data, 0.50))
plot_model(colony_uptake(15.0, k_opt2, 0.50))

# Print results
results = ["Real" f(k_data) k_data[1] k_data[2];
            "Initial guess" f(k_init) k_init[1] k_init[2];
            "Optimized model (defualt)" f(k_opt1) k_opt1[1] k_opt1[2];
            "Optimized model (PSO)" f(k_opt2) k_opt2[1] k_opt2[2];]

pretty_table(results , ["State" "Cost" "V_max" "k_m"])
