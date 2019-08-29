using Optim
using Plots
using PrettyTables
using JLD2
include("colony-uptake.jl")
V_max, k_m, k_CA = [3.5676e-1, 3.9*Av*V_cell*60, 0.1]

function target_data(time_points, output_int, p::Array{})
    # Returns inf if the solution does not converge.
    data = []
    try
        data = colony_model(time_points[end]+1, p)(time_points)[output_int,:]
    catch e
        @warn "Did not converge."
        @warn e
        data = repeat([Inf], length(time_points))
    end
    return data
end

function opt_particle_swarm(k_init, iterations=1000, timelimit=1000)
    lower_bound = [V_max, k_m, k_CA].*1e-2
    upper_bound = [V_max, k_m, k_CA].*1e2
    options = Optim.Options(show_trace=true,
                            show_every=100,
                            iterations=iterations,
                            extended_trace=true,
                            store_trace=true,
                            time_limit=timelimit)

    alg = ParticleSwarm(lower=lower_bound, upper=upper_bound,n_particles=6)
    time = @elapsed res = optimize(f, k_init, alg, options)
    return Optim.minimizer(res), res.trace, time
end

function opt_auto(k_init, iterations=1000, timelimit=1000)
    options = Optim.Options(show_trace=true,
                            extended_trace=true,
                            store_trace=true,
                            show_every=100,
                            iterations=iterations,
                            time_limit = timelimit)
    time = @elapsed res = optimize(f, k_init, options)
    return Optim.minimizer(res), res.trace, time
end

function plot_trace(trace::Array{Optim.ParticleSwarmState}, label="", t0=0)
    particles = trace[1].n_particles
    len = length(trace)
    x = Array{Float64, 2}(undef, len, particles)
    y = Array{Float64, 2}(undef, len, particles)
    z = Array{Float64, 2}(undef, len, particles)
    p = plot(xlims=[trace[1].lower[1], trace[1].upper[1]],
        ylims=[trace[1].lower[2], trace[1].upper[2]],
        zlims=[trace[1].lower[3], trace[1].upper[3]])
    for i in 1:len
        for p in 1:particles
            x[i,p] = trace[i].X[1,p]
            y[i,p] = trace[i].X[2,p]
            z[i,p] = trace[i].X[3,p]
        end
    end
    plot!(x,y,z)
    return p, x, y, z
end

function plot_trace(trace, label="", t0=0)
    len = length(trace)
    t = zeros(len)
    x = zeros(len)
    for i in 1:len
        t[i] = trace[i].metadata["time"] + t0
        x[i] = trace[i].value
    end
    plot!(t,x, label=label, yscale=:log10)
end

function get_cost(args...)
    time = []
    cost = []
    for trace in args
        for it in trace
            if (trace != args[1])
                t0 = args[1][end].metadata["time"]
                display(t0)
            else t0 = 0 end
            push!(time, it.metadata["time"] + t0)
            push!(cost, it.value)
        end
    end
    return [time,cost]
end

loss_function(data1, data2) = 1.0/length(data1)*sum(abs.(data1 - data2))

function opt_bech_test(iterations)
    nm_trace = []
    pso_trace = []
    pso_nm_trace = []
    k_initOLD = k_init
    p = plot()
    for i = 1:iterations
        k_init .*= rand(3)/10
        push!(nm_trace, opt_auto(k_init, 5000, 5000)[2])
        push!(pso_trace, opt_particle_swarm(k_init, 5000, 5000)[2])
        #k_init3, trace0, time31 = opt_particle_swarm(k_init, 600)
        #push!(pso_nm_trace, [trace0, opt_auto(k_init3, 400)[2]])
    end

    return nm_trace, pso_trace, pso_nm_trace
end
k_init = [V_max, k_m, k_CA]
k_data = [V_max, k_m, k_CA]*10
t = collect(0.0:1:3.0)
data1 =  target_data(t, 2, [(:V_max, k_data[1]), (:k_m, k_data[2]), (:k_CA, k_data[3])])

# Cost function for optimizing model, using different concentrations
f(x) = loss_function(data1, target_data(t, 2, [(:V_max, x[1]), (:k_m, x[2]), (:k_CA, x[3])]))
# Optimize using different strategies
nm, pso, pso_nm = opt_bech_test(10)
pso_cost = get_cost.(pso)
nm_cost = get_cost.(nm)
limit(costs, t) = (costs[1][1:t], costs[2][1:t])

@gif for i in 2:1000
    plot(limit.(pso_cost,i), xscale=:log10, yscale=:log10, color=:blue, legend=false, ylabel="Cost", xlabel="Seconds")
    plot!(limit.(nm_cost,i) ,yscale=:log10, color=:red)
end

k_opt1, trace1, time1 = opt_auto(k_init, 60)
k_opt2, trace2, time2 = opt_particle_swarm(k_init, 1000)
k_init3, trace31, time31 = opt_particle_swarm(k_init, 1000)
k_opt3, trace32, time32 = opt_auto(k_init3, 1000)
p, x, y, z = plot_trace(swarm_trace)
x
trace = swarm_trace
swarm_trace[200]
plot(xlims=[trace[1].lower[1], trace[1].upper[1]],
    ylims=[trace[1].lower[2], trace[1].upper[2]],
    zlims=[trace[1].lower[3], trace[1].upper[3]])
pyplot()
plot(z, yscale=:log10)
plot_trace(trace2)
plot()
@gif for i in 1:length(x[:,1])
    scatter!(x[i,:]',y[i,:]',z[i,:]', legend=false)
end every 1000

# Plot
plot()
plot_trace(trace1, "Neader-Mead")
plot_trace(trace2, "Particle Swarm")
plot_trace(trace31, "Particle Swarm")
plot_trace(trace32, "Neader-Mead", trace32[end].metadata["time"])
plot!(xlims=[0,2000])

# Print results
results = ["Real" f(k_data) k_data[1] k_data[2] k_data[3] "-";
            "Guess" f(k_init) k_init[1] k_init[2] k_init[3] "-";
            "Defualt" f(k_opt1) k_opt1[1] k_opt1[2] k_opt1[3] time1;
            "PSO" f(k_opt2) k_opt2[1] k_opt2[2] k_opt2[3] time2;
            "PSO+NM)" f(k_opt3) k_opt3[1] k_opt3[2] k_opt3[3] (time31+time32);]

function Base.show(io::IO, x::Union{Float64,Float32})
    Base.Grisu._show(io, round(x, sigdigits=3), Base.Grisu.SHORTEST, 0, get(io, :typeinfo, Any) !== typeof(x), false)
end
pretty_table(results , ["State" "Cost" "V_max" "k_m" "k_CA" "time"])
