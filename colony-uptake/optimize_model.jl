using Optim
using Plots
using PrettyTables
include("colony-uptake.jl")

# Returns the number of ions in the growth medium according to the model.
function ions_in_medium(t, k, c)
    sol = colony_uptake(18.0, k, c)
    return sol(t)[2,:]
end

# Using leat squares
loss_function(data1, data2) = sum((data1 - data2).^2)

# Initial gues for the unknown constants
k_init = [V_max, k_m]

# Generated data using other constants
t = 0.0:1:5.0
k_data = [V_max*5, k_m]
init = ions_in_medium(t, k_init, 0.25)
data1 = ions_in_medium(t, k_data, 0.25)
data2 = ions_in_medium(t, k_data, 0.50)
data3 = ions_in_medium(t, k_data, 0.125)

# Cost function for optimizing model, using different concentrations
f(x) = (loss_function(data1, ions_in_medium(t,x, 0.25))
    + loss_function(data2, ions_in_medium(t,x,0.5))
    + loss_function(data3, ions_in_medium(t,x,0.125)))

# Optimize
res = optimize(f, k_init, NelderMead(), Optim.options(  g_tol=1e-9,
                                                        show_trace=true))
k_opt = Optim.minimizer(res)

# Plot
plot_init()
plot_model(colony_uptake(15.0, k_init, 0.50))
plot_model(colony_uptake(15.0, k_data, 0.50))
plot_model(colony_uptake(15.0, k_opt, 0.50))

# Print results
results = ["Real" f(k_data) k_data[1] k_data[2];
            "Initial guess" f(k_init) k_init[1] k_init[2];
            "Optimized model" f(k_opt) k_opt[1] k_opt[2];]
print(res)
pretty_table(results , ["State" "Cost" "V_max" "k_m"])
