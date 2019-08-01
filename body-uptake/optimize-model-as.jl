using Optim
using Plots
using PrettyTables
include("uptake-as.jl")

As_molecular_weight = 74.92159 #g/mol

loss_function(data1,data2) = sum(abs.(data1 - data2)) / length(data1)

#experimental data, total dose of arsenic (μg) excretion in urine after a single dose of ingestion of 500μg in day 1,2,3,4 (J. P. Buchet et al 1981)
As3_urine = [110.7,66.0,45.1,32.1] #after ingestion of NaAsO₂
AsM_urine = [359.1,30.0,18.5,12.5] #after ingestion of CH₃AsO₃Na₂
AsN_urine = [286.4,72.5,32.6,12.6] #after ingestion of (CH₃)₂AsO₂Na

function accumulating(array) #(also converts from μg to μmol)
    for i = 2:length(array)
        array[i] = array[i] + array[i-1]
    end
    return array / As_molecular_weight
end

As3_real = accumulating(As3_urine)
AsM_real = accumulating(AsM_urine)
AsN_real = accumulating(AsN_urine)

function urine_single_dose(intake,type,days,timestep,K) #intake in μmol
    initial = zeros(48)
    if type == "As3"
        initial[1] = intake
    elseif type == "As5"
        initial[13] = intake
    elseif type == "AsM"
        initial[25] = intake
    elseif type == "AsN"
        initial[37] = intake
    end
    tspan = (0.0,days)
    prob = ODEProblem(arsenic_once,initial,tspan,K)
    sol = solve(prob)
    urine = []
    for i in timestep
        a = sol(i)[10] + sol(i)[22] + sol(i)[34] + sol(i)[46]
        push!(urine,a)
    end
    return urine  #in μmol
end

function total_loss(K)
    timestep = [1,2,3,4]
    As3_model = urine_single_dose(500 / As_molecular_weight, "As3", 4, timestep,K) #500μg
    AsM_model = urine_single_dose(500 / As_molecular_weight, "AsM", 4, timestep,K)
    AsN_model = urine_single_dose(500 / As_molecular_weight, "AsN", 4, timestep,K)
    total_loss = loss_function(As3_real,As3_model) + loss_function(AsM_real,AsM_model) + loss_function(AsN_real,AsN_model)
    return total_loss

end

#total_loss([k_GBₘ, k_BG])

function opt_particle_swarm(k_init)
    res = optimize(total_loss, k_init, ParticleSwarm(lower=[0.0],
                                            upper=[1e4],
                                            n_particles=7))

    return Optim.minimizer(res)
end

print( opt_particle_swarm([100.0]))
print(total_loss(opt_particle_swarm([100.0])))

function plotting()
    tspan = (0.0,4)
    initial = zeros(48)
    initial[25] = 500 / As_molecular_weight
    prob1 = ODEProblem(arsenic_once,initial,tspan,eUₘ)
    sol1 = solve(prob1)
    urine1 = sol1[10,:] + sol1[22,:] + sol1[34,:] + sol1[46,:]
    plot()
    plot!(sol1.t, urine1, label = "model")

    timestep = [1,2,3,4]
    plot!(timestep, AsM_real, label = "real")

    prob2 = ODEProblem(arsenic_once,initial,tspan,opt_particle_swarm([100.0]))
    sol2 = solve(prob2)
    urine2 = sol2[10,:] + sol2[22,:] + sol2[34,:] + sol2[46,:]
    plot!(sol2.t, urine2, label = "optimized model")
end
plotting()
