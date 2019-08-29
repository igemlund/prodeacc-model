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

    df = DataFrame(sol')
    [df[col] ./= max_values[col]' for col in 1:length(max_values)]
    [df[i] = map(x -> float_to_string(round(x, digits=4)), df[i]) for i in 1:length(max_values)]
    display(df)
    CSV.write("$name.csv", df)#, head=false)
end

function float_to_string(float)
    s = string(float)
    i = findfirst("e", s)

    if s == "NaN"
        return "0.0"
    end

    if i == nothing
        return s
    else
        i = i[1]
    end

    if float >= 1
        exp = parse(Int8, s[i+1:end])
        s = s[1:i-1]
        s *= repeat("0", exp)
    else
        exp = parse(Int8, s[i+2:end])
        if s[1] != "0"
            exp -=1
            s = "0."*repeat("0", exp)*s[1]*s[3:i-1]
        end
    end
    return s
end
