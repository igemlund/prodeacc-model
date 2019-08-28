using CSV
using DataFrames

# Takes a diffeq solution, name, and n_values and export the data  as csv
# Used for exporting data for animation in After Effects.
function export_to_after_effects(sol, name, n_values)
    delta_t = sol.t[end]/n_values
    t = collect(0:delta_t:sol.t[end])[1:end-1]'
    u = sol.(t)
    max_values = zeros(length(u[1]))
    for step in u
        for index in 1:length(max_values)
            if max_values[index] < step[index]
                max_values[index] = step[index]
            end
        end
    end

    display(u)
    foreach(x -> x ./=max_values, u)
    display(u)
    df = DataFrame(sol')
    [df[col] ./= max_values[col]' for col in 1:length(max_values)]
    [df[i] = map(x -> round(x, digits=10), df[i]) for i in 1:length(max_values)]
    display(df)
    CSV.write("$name.csv", df)#, head=false)
end

export_to_after_effects(sol, "test_round", 400)
