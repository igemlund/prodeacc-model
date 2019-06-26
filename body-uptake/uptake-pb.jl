using DifferentialEquations
using ParameterizedFunctions
using Plots

#Based on dede et al 2018
Q_liv = 1084.3 #Blood plasma flow to liver (L/day)
Q_kid = 737.3 #Blood plasma flow to kidney (L/day)
Q_bone = 130.1
Q_wp = 86.74 #Blood plasma flow to well-perfused kidney (L/day)
Q_pp = 346.98 #Blood plasma flow to poorly-perfused kidney

V_liv = 1.75 #Volume of liver (L)
V_kid = 0.29 #kidney (L)
V_bone = 9.8
V_blood = 5.11
V_plas = 2.81 #plasma
V_wp = 1.96
V_pp = 58.24

P_liv = 100 #Partion coefficient, liver/plasma (no unit)
P_kid = 100
P_wp = 100
P_pp = 20
P_bone = 1000

#transfer rate from organs to blood-plasma (per day)
k_LiB = Q_liv/(V_liv*P_liv)
k_KB = Q_kid/(V_kid*P_kid)
k_BoB = Q_bone/(V_bone*P_bone)
k_WB = Q_wp/(V_wp*P_wp)
k_PB = Q_pp/(V_pp*P_pp)

#transfer rate from blood-plasma to organs (per day)
k_BLi = Q_liv/V_plas
k_BK = Q_kid/V_plas
k_BW = Q_wp/V_plas
k_BBo = Q_bone/V_plas
k_BP = Q_pp/V_plas

A_gi = 0.1 #PB absorption coefficient from GI tract, ranges from 0.06 to 0.12 (unitless)
eU = 0.47 #(/day)
eB = 0.2 #(/day)

HCT = 0.45 # haematocrit fraction of whole blood
R = 1.2 #ratio of unbound erythrocyte PB concentration to plasma PB concentration
BIND = 0.437 #Pb binding capacity of erythrocytes (mg Pb L⁻¹ cell)
KBIND = 3.72e-4 #Binding constant of erythrocytes (mg Pb L⁻¹ cell)

lead = @ode_def begin #ODE system for constant intake of lead per day over several days
    dG = - A_gi * IR_gi - (1 - A_gi) * IR_gi  #IR_gi = oral intake rate  of Pb (mg/day) (given in solving function) #amount in gut
    dLi = A_gi * IR_gi  - (k_LiB + eB) * Li + k_BLi * B #amount in liver
    dK = - (k_KB + eU) * K + k_BK * B #amount in kidney
    dBo = - k_BoB * Bo + k_BBo * B #amount in bone
    dWP = - k_WB * WP + k_BW * B #amount in well-perfused tissues
    dPP = - k_PB * PP + k_BP * B #amount in poorly-perfused tissues
    dB = (- (k_BLi + k_BK + k_BBo + k_BW + k_BP) * B #amount in blood plasma
     + k_LiB * Li + k_KB * K + k_BoB * Bo + k_WB * WP + k_PB * PP)
    dU = eU * K #urine
    dBi = eB * Li #Biliary excretion
    dF = (1 - A_gi) * IR_gi #Faeces
end a

lead_once = @ode_def begin #ODE system for one time intake of lead
    dG = - A_gi * G - (1 - A_gi) * G
    dLi = A_gi * G  - (k_LiB + eB) * Li + k_BLi * B
    dK = - (k_KB + eU) * K + k_BK * B
    dBo = - k_BoB * Bo + k_BBo * B
    dWP = - k_WB * WP + k_BW * B
    dPP = - k_PB * PP + k_BP * B
    dB = (- (k_BLi + k_BK + k_BBo + k_BW + k_BP) * B
     + k_LiB * Li + k_KB * K + k_BoB * Bo + k_WB * WP + k_PB * PP) #blood plasma amount
    dU = eU * K
    dBi = eB * Li
    dF = (1 - A_gi) * G
end b



function CB(CPLASMA) #Function to calculate whole blood concentration from plasma concentration
    CB = (((1 - HCT) * CPLASMA)
        + (HCT * CPLASMA * (R + (BIND/(KBIND + CPLASMA)))))
    return CB
end

function wholeblood_conc(solution) #take an array of plasma concentration and turn it into an array of wholeblood concentration (μg/dL)
    wholeblood = []
    for i in solution[7,:]
        push!(wholeblood,CB(i/(V_plas*10)) * 1000) #(μg dL⁻¹)
    end
    return wholeblood
end

function solving(intake,initial,tspan) #solving for lead ODE
    global IR_gi =  intake
    prob = ODEProblem(lead,initial,tspan,1)
    sol = solve(prob)
    return sol
end

#plotting function for lead ODE
function plotting(intake, duration, range, compartment) #intake is in mg/day, duration is in days, compartment is the order of compartment in the differential equations (from 1 to 11)
    A = solving(intake,zeros(10),(0.0,duration))
    B = solving(0.0, A[end],(duration,range))

    if compartment == 11   #compartment 11 is whole blood concentration
        plot()
        plot!(A.t, wholeblood_conc(A))
        plot!(B.t, wholeblood_conc(B))
    else
        plot()
        plot!(A.t, A[compartment,:])
        plot!(B.t, B[compartment,:])
    end
end

function plotting_once(intake,tspan) #plotting function for lead_once
    x_0 = zeros(10)
    x_0[1] = intake
    prob = ODEProblem(lead_once,x_0,tspan,1)
    sol = solve(prob)
    plot()
    plot!(sol.t,wholeblood_conc(sol))
end

#plotting_once(68e-3,(0.0,10.0)) #figure 10 subject C in dede et al
plotting(105e-3,82.0,300.0,11) #figure 10 subject D in dede et al
