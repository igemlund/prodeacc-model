# Takes two solution objects from DifferentialEquations.jl and retuns
# a merged one.
function sol_append(sol1, sol2)
    #sol2.t .+= sol1.t[end]
    if sol1.dense || sol2.dense
        @warn "Option dense=true, set it to false"
    end
    append!(sol1.u, sol2.u[2,:])
    append!(sol1.t, sol2.t[2,:])
    append!(sol1.alg_choice, sol2.alg_choice[2,:])
    return sol1
end
