# Based on Dede et al 2018
using Plots
using DifferentialEquations
using ParameterizedFunctions
using IterableTables, DataFrames

Q_gi = 1440  #Blood flow to GI tract (l/day)
Q_liv = 446.4 #Blood flow to liver (l/day)
Q_kid = 1440 #Blood flow to kidney (L/day)
Q_skin = 374.4 #Blood flow to skin (l/day)
Q_lung = 7488 #Blood flow to lung (L/day)
Q_mus = 2592 #Blood flow to muscle (L/day)
Q_hrt = 288 #Blood flow to heart (L/day)
Q_br = 907.2 #Blood flow to brain (L/day)

V_gi = 1.2 #Volume of GI tract (L)
V_liv = 1.82 #Volume of liver (L)
V_kid = 0.28 #Volume of kidney (L)
V_skin = 2.6 #Volume of skin (L)
V_lung = 0.56 #Volume of lung (L)
V_blood = 5.53 #Volume of blood (L)
V_mus = 55.5 #Volume of muscle (L)
V_hrt = 0.35 #Volume of heart (L)
V_br = 1.4 #Volume of brain (L)

#Inorganic patition coefficient
P_gi3 = 8.3 #Gi tract/blood partition coefficient for As(III) (no unit)
P_gi5 = 2.7 #Gi tract/blood partition coefficient for As(V) (no unit)
P_liv3 = 16.5 #Same as above but for liver
P_liv5 = 15.8
P_kid3 = 20 #Same as above but for kidney
P_kid5 = 8.3
P_skin3 = 7.4 #Same as above but for skin
P_skin5 = 7.9
P_lung3 = 6.7 #Same as above but for lung
P_lung5 = 2.1
P_mus3 = 7.4 #Same as above but for muscle
P_mus5 = 7.9
P_hrt3 = 7.4 #Same as above but for heart
P_hrt5 = 7.9
P_br3 = 2.4 #Same as above but for brain
P_br5 = 2.4

#Organic partition coefficient (MMA and DMA are two forms of organic arsenic)
P_giM = 2.2 #GI tract/blood partition coefficient for MMA
P_giN = 2.1 #GI tract/blood partition coefficient for DMA
P_livM = 3.3 #Same as above but for liver
P_livN = 3.3
P_kidM = 100 #Same as above but for kidney
P_kidN = 1
P_skinM = 2.61 #Same as above but for skin
P_skinN = 2.4
P_lungM = 1.3 #Same as above but for skin
P_lungN = 1.3
P_musM = 2.61 #Same as above but for muscle
P_musN = 2.4
P_hrtM = 2.61 #Same as above but for heart
P_hrtN = 2.4
P_brM = 2.2 #Same as above but for brain
P_brN = 3.3

eF = 0.03 #faecal excretion rate (/day)
eB = 0.43 #Biliary excretion rate (/day)
eUᵢ = 57.6 #Urinary excretion rate for As3 and As5 (/day)
eUₘ = 432 #Urinary excretion rate for MMA (/day)
eUₙ = 100.8 #Urinary excretion rate for DMA (/day)


#Rate constant of As(III) (/day)
k_GB₃ = Q_gi/(V_gi*P_gi3) #from GI tract to blood
k_LiB₃ = Q_liv/(V_liv*P_liv3) #from liver to blood
k_KB₃ = Q_kid/(V_kid*P_kid3) #from kidney to blood
k_SB₃ = Q_skin/(V_skin*P_skin3) #from skin to blood
k_LuB₃ = Q_lung/(V_lung*P_lung3) #from lung to blood
#R stand for the rest, it's a combined compartment for muscle, heart and brain since arsenic don't have as big of an effect on them as the other organs so we are not interested on how they behave
k_RB₃ = Q_mus/(V_mus*P_mus3) + Q_hrt/(V_hrt*P_hrt3) + Q_br/(V_br*P_br3)

#Rate constant of As(V) (/day)
k_GB₅ = Q_gi/(V_gi*P_gi5) #from GI tract to blood
k_LiB₅ = Q_liv/(V_liv*P_liv5) #from liver to blood
k_KB₅ = Q_kid/(V_kid*P_kid5) #from kidney to blood
k_SB₅ = Q_skin/(V_skin*P_skin5) #from skin to blood
k_LuB₅ = Q_lung/(V_lung*P_lung5) #from lung to blood
k_RB₅ = Q_mus/(V_mus*P_mus5) + Q_hrt/(V_hrt*P_hrt5) + Q_br/(V_br*P_br5) #from the rest to blood

#Rate constant of MMA (/day)
k_GBₘ = Q_gi/(V_gi*P_giM) #from GI tract to blood
k_LiBₘ = Q_liv/(V_liv*P_livM) #from liver to blood
k_KBₘ = Q_kid/(V_kid*P_kidM) #from kidney to blood
k_SBₘ = Q_skin/(V_skin*P_skinM) #from skin to blood
k_LuBₘ = Q_lung/(V_lung*P_lungM) #from lung to blood
k_RBₘ = Q_mus/(V_mus*P_musM) + Q_hrt/(V_hrt*P_hrtM) + Q_br/(V_br*P_brM) #from the rest to blood

#Rate constant of DMA (/day)
k_GBₙ = Q_gi/(V_gi*P_giN) #from GI tract to blood
k_LiBₙ = Q_liv/(V_liv*P_livN) #from liver to blood
k_KBₙ = Q_kid/(V_kid*P_kidN) #from kidney to blood
k_SBₙ = Q_skin/(V_skin*P_skinN) #from skin to blood
k_LuBₙ = Q_lung/(V_lung*P_lungN) #from lung to blood
k_RBₙ = Q_mus/(V_mus*P_musN) + Q_hrt/(V_hrt*P_hrtN) + Q_br/(V_br*P_brN) #from the rest to blood

#Rate constant from blood to other organs for both organic and inorganic arsenic
k_BG = Q_gi/V_blood#from Blood to GI tract
k_BLi = Q_liv/V_blood #from blood to liver
k_BK = Q_kid/V_blood #from blood to kidney
k_BS = Q_skin/V_blood #from blood to skin
k_BLu = Q_lung/V_blood #from blood to lung
k_BR = (Q_mus+Q_hrt+Q_br)/V_blood #from blood to the rest

#Reduction and Oxidation of inorganic arsenic
R_Tis = 32.88 #As(V) reduction in tissues (/day)
R_K = 42 #As(V) reduction in kidney (/day)
O_Tis = 43.92 #As(III) oxidation is tissues (/day)

#Metabolism constant for methylation of arsenic
M_Li3M = 763 #from As3 to MMA in the liver (μmol/day)
M_Li3N = 2880 #from As3 to DMA in the liver (μmol/day)
M_LiMN = 950.4 #from MMA to DMA (μmol/day)
M_K3M = 305.2 #from As3 to MMA in the kidney (μmol/day)
M_K3N = 1152 #from As3 to DMA in the kidney (μmol/day)
M_KMN = 380.16 #from MMA to DMA in the kidney (μmol/day)
M = 3 #another metabolism constant (μmol/L), all of them have value 3 so it can be combined into one variable
#G₅ = 1.67
dede = @ode_def begin

#Differential equations for As(III)
    dG₃ = - (k_GB₃ + eF + O_Tis) * G₃ + k_BG * B₃ + R_Tis * G₅
    dLi₃ = - (k_LiB₃ + eB + O_Tis) * Li₃ + k_BLi * B₃ + R_Tis * Li₅ - (M_Li3M * (Li₃ / V_liv)) / (M + Li₃ / V_liv) - (M_Li3N * (Li₃ / V_liv)) / (M + Li₃ / V_liv)
    dK₃ = - (k_KB₃ + eUᵢ + O_Tis) * K₃ + k_BK * B₃ + R_K * K₅ - (M_K3M *( K₃ / V_kid))/(M + (K₃ / V_kid)) - (M_K3N * (K₃ / V_kid)) / (M + K₃ / V_kid)
    dS₃ = - (k_SB₃ + O_Tis) * S₃ + k_BS * B₃ + R_Tis * S₅
    dLu₃ = - (k_LuB₃ + O_Tis) * Lu₃ + k_BLu * B₃ + R_Tis * Lu₅
    dR₃ = - (k_RB₃ + O_Tis) * R₃ + k_BR * B₃ + R_Tis * R₅
    dB₃ = - (k_BG + k_BLi + k_BK + k_BS + k_BLu + k_BR + O_Tis) * B₃ + k_GB₃ * G₃ + k_LiB₃ * Li₃ + k_KB₃ * K₃ + k_SB₃ * S₃ + k_LuB₃ * Lu₃ + k_RB₃ * R₃ + R_Tis * B₅
    dU₃ = eUᵢ * K₃
    dF₃ = eF * G₃
    dBil₃ = eB * Li₃

#Differential equations for As(V)
    dG₅ = - (k_GB₅ + eF + R_Tis) * G₅ + k_BG * B₅ + O_Tis * G₃
    dLi₅ = - (k_LiB₅ + eB + R_Tis) * Li₅ + k_BLi * B₅ + O_Tis * Li₃
    dK₅ = - (k_KB₅ + eUᵢ + R_K) * K₅ + k_BK * B₅ + O_Tis * K₃
    dS₅ = - (k_SB₅ + R_Tis) * S₅ + k_BS * B₅ + O_Tis * S₃
    dLu₅ = - (k_LuB₅ + R_Tis) * Lu₅ + k_BLu * B₅ + O_Tis * Lu₃
    dR₅ = - (k_RB₅ + R_Tis) * R₅ + k_BR * B₅ + O_Tis * R₃
    dB₅ = - (k_BG + k_BLi + k_BK + k_BS + k_BLu + k_BR + R_Tis) * B₅ + k_GB₅ * G₅ + k_LiB₅ * Li₅ + k_KB₅ * K₅ + k_SB₅ * S₅ + k_LuB₅ * Lu₅ + k_RB₅ * R₅ + O_Tis * B₃
    dU₅ = eUᵢ * K₅
    dF₅ = eF * G₅
    dBil₅ = eB * Li₅

#Differential equations for MMA
    dGₘ = - (k_GBₘ + eF) * Gₘ + k_BG * Bₘ
    dLiₘ = - (k_LiBₘ + eB) * Liₘ + k_BLi * Bₘ + (M_Li3M * (Li₃ / V_liv)) / (M + Li₃ / V_liv) - (M_LiMN * (Liₘ / V_liv)) / (M + Liₘ / V_liv)
    dKₘ = - (k_KBₘ + eUₘ) * Kₘ + k_BK * Bₘ + (M_K3M * (K₃ / V_kid)) / (M + K₃ / V_kid) - (M_KMN * (Kₘ / V_kid)) / (M + Kₘ / V_kid)
    dSₘ = - k_SBₘ * Sₘ + k_BS * Bₘ
    dLuₘ = - k_LuBₘ * Luₘ + k_BLu * Bₘ
    dRₘ = - k_RBₘ * Rₘ + k_BR * Bₘ
    dBₘ = - (k_BG + k_BLi + k_BK + k_BS + k_BLu + k_BR) * Bₘ + k_GBₘ * Gₘ + k_LiBₘ * Liₘ + k_KBₘ * Kₘ + k_SBₘ * Sₘ + k_LuBₘ * Luₘ + k_RBₘ * Rₘ
    dUₘ = eUₘ * Kₘ
    dFₘ = eF * Gₘ
    dBilₘ = eB * Liₘ

#Differential equations for DMA
    dGₙ = - (k_GBₙ + eF) * Gₙ + k_BG * Bₙ
    dLiₙ = - (k_LiBₙ + eB) * Liₙ + k_BLi * Bₙ + (M_Li3N * (Li₃ / V_liv)) / (M + Li₃ / V_liv) + (M_LiMN * (Liₘ / V_liv)) / (M + Liₘ / V_liv)
    dKₙ = - (k_KBₙ + eUₙ) * Kₙ + k_BK * Bₙ + (M_K3N * (K₃ / V_kid)) / (M + K₃ / V_kid) + (M_KMN * (Kₘ / V_kid)) / (M + Kₘ / V_kid)
    dSₙ = - k_SBₙ * Sₙ + k_BS * Bₙ
    dLuₙ = - k_LuBₙ * Luₙ + k_BLu * Bₙ
    dRₙ = - k_RBₙ * Rₙ + k_BR * Bₙ
    dBₙ = - (k_BG + k_BLi + k_BK + k_BS + k_BLu + k_BR) * Bₙ + k_GBₙ * Gₙ + k_LiBₙ * Liₙ + k_KBₙ * Kₙ + k_SBₙ * Sₙ + k_LuBₙ * Luₙ + k_RBₙ * Rₙ
    dUₙ = eUₙ * Kₙ
    dFₙ = eF * Gₙ
    dBilₙ = eB * Liₙ
end a

#dose = 1.67
#time = 5 #days
#dede_0 = [dose;zeros(39)]
#tspan = (0.0,time)
#prob = ODEProblem(dede,dede_0,tspan,1)
#sol = solve(prob)
#display(sol[28,:])
y = []
my_x = Array{Float64,1}[]
dose = [1.67]
for i = 0:1
    dede_0 = [zeros(10);dose[trunc(Int,i+1)];zeros(29)]
    tspan = (i, i + 0.99)
    prob = ODEProblem(dede,dede_0,tspan,1)
    sol = solve(prob)
    push!(dose, sol[11,length(sol.t)]+1.67)
    append!(y,sol[28,:])
    #append!(x,sol.t)

end

display(y)
#plot!(sol.t,sol[28,:])
#plot!(sol.t, sol[38,:])
#plot!(sol.t, sol[8,:] + sol[18,:])
#plot!(sol.t, sol[8,:] + sol[18,:] + sol[28,:] + sol[38,:])
