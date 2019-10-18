using CSV
using DataFrames

# Takes a diffeq solution, name, and n_values and export the data  as csv
# Used for exporting data for animation in After Effects.
function export_to_csv(sol, name, n_values)
    delta_t = sol.t[end]/n_values
    t = collect(0:delta_t:sol.t[end])[1:end-1]'
    u = sol.(t)

    df = DataFrame(u[1]')
    for i in 2:n_values
        push!(df, u[i]')
    end
    CSV.write("csv/$name.csv", df)#, head=false)
end

function normalize_csv(name)
    df = CSV.read(name, copycols=true)
    max_values = zeros(length(df[1,:]))
    println(length(df[1,:]))
    println(length(df[:,1]))
    for i in 2:length(df[:,1])
        step = df[i,1:end]
        for col in 1:length(max_values)
            if max_values[col] < step[col]
                max_values[col] = step[col]
            end
        end
        #CSV.write("$name.csv", df)#, head=false)
    end

    [df[col] ./= max_values[col]' for col in 1:length(max_values)]
    [df[i] = map(x -> float_to_string(round(x, digits=4)), df[i]) for i in 1:length(max_values)]
    CSV.write("$name-normalized.csv", df)
end
export_to_csv(sol1, "11mcg_pb", 8000)
normalize_csv("csv/11mcg_pb.csv")

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

#export_to_after_effects(sol1, "pb_114mcg_per_day_8y_PRODEACC", 16000)
#export_to_after_effects(sol2, "pb_114mcg_per_day_8y", 16000)
