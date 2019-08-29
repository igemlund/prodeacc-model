# Based on Dede et al 2018
using Plots
using DifferentialEquations
using ParameterizedFunctions
using IterableTables
include("../export_data.jl")

const Q_gi = 1440  #Blood flow to GI tract (l/day)
const Q_liv = 446.4 #Blood flow to liver (l/day)
const Q_kid = 1440 #Blood flow to kidney (L/day)
const Q_skin = 374.4 #Blood flow to skin (l/day)
const Q_lung = 7488 #Blood flow to lung (L/day)
const Q_mus = 2592 #Blood flow to muscle (L/day)
const Q_hrt = 288 #Blood flow to heart (L/day)
const Q_br = 907.2 #Blood flow to brain (L/day)

const V_gi = 1.2 #Volume of GI tract (L)
const V_liv = 1.82 #Volume of liver (L)
const V_kid = 0.28 #Volume of kidney (L)
const V_skin = 2.6 #Volume of skin (L)
const V_lung = 0.56 #Volume of lung (L)
const V_blood = 5.53 #Volume of blood (L)
const V_mus = 55.5 #Volume of muscle (L)
const V_hrt = 0.35 #Volume of heart (L)
const V_br = 1.4 #Volume of brain (L)

#Inorganic patition coefficient
const P_gi3 = 8.3 #Gi tract/blood partition coefficient for As(III) (no unit)
const P_gi5 = 2.7 #Gi tract/blood partition coefficient for As(V) (no unit)
const P_liv3 = 16.5 #Same as above but for liver
const P_liv5 = 15.8
const P_kid3 = 20 #Same as above but for kidney
const P_kid5 = 8.3
const P_skin3 = 7.4 #Same as above but for skin
const P_skin5 = 7.9
const P_lung3 = 6.7 #Same as above but for lung
const P_lung5 = 2.1
const P_mus3 = 7.4 #Same as above but for muscle
const P_mus5 = 7.9
const P_hrt3 = 7.4 #Same as above but for heart
const P_hrt5 = 7.9
const P_br3 = 2.4 #Same as above but for brain
const P_br5 = 2.4

#Organic partition coefficient (MMA and DMA are two forms of organic arsenic)
#partition coefficient is the ratio of concentrations of a compound in a mixture of two immiscible phases at equilibrium
const P_giM = 2.2 #GI tract/blood partition coefficient for MMA ((GI/V_gi)/(B/V_blood))
const P_giN = 2.1 #GI tract/blood partition coefficient for DMA
const P_livM = 3.3 #Same as above but for liver
const P_livN = 3.3
const P_kidM = 100 #Same as above but for kidney
const P_kidN = 1
const P_skinM = 2.61 #Same as above but for skin
const P_skinN = 2.4
const P_lungM = 1.3 #Same as above but for skin
const P_lungN = 1.3
const P_musM = 2.61 #Same as above but for muscle
const P_musN = 2.4
const P_hrtM = 2.61 #Same as above but for heart
const P_hrtN = 2.4
const P_brM = 2.2 #Same as above but for brain
const P_brN = 3.3

#calculate absorption coefficient (A_gi) from P_gi, A_gi is the percentage of amount of arsenic being obsorbed to blood
const A_gi3 = (V_blood/V_gi)/(P_gi3 + (V_blood/V_gi))
const A_gi5 = (V_blood/V_gi)/(P_gi5 + (V_blood/V_gi))
const A_giM = (V_blood/V_gi)/(P_giM + (V_blood/V_gi))
const A_giN = (V_blood/V_gi)/(P_giN + (V_blood/V_gi))

const eF = 0.03 #faecal excretion rate (/day)
const eB = 0.43 #Biliary excretion rate (/day)
const eUᵢ = 100.8 #Urinary excretion rate for As3 and As5 (/day)
const eUₘ = 432 #Urinary excretion rate for MMA (/day)
const eUₙ = 187.2 #Urinary excretion rate for DMA (/day)


#Rate constant of As(III) (/day)
const k_GB₃ = Q_gi/(V_gi*P_gi3) #from GI tract to blood
const k_LiB₃ = Q_liv/(V_liv*P_liv3) #from liver to blood
const k_KB₃ = Q_kid/(V_kid*P_kid3) #from kidney to blood
const k_SB₃ = Q_skin/(V_skin*P_skin3) #from skin to blood
const k_LuB₃ = Q_lung/(V_lung*P_lung3) #from lung to blood
const k_MB₃ = Q_mus/(V_mus*P_mus3)
const k_HB₃ = Q_hrt/(V_hrt*P_hrt3)
const k_BrB₃ = Q_br/(V_br*P_br3)

#Rate constant of As(V) (/day)
const k_GB₅ = Q_gi/(V_gi*P_gi5) #from GI tract to blood
const k_LiB₅ = Q_liv/(V_liv*P_liv5) #from liver to blood
const k_KB₅ = Q_kid/(V_kid*P_kid5) #from kidney to blood
const k_SB₅ = Q_skin/(V_skin*P_skin5) #from skin to blood
const k_LuB₅ = Q_lung/(V_lung*P_lung5) #from lung to blood
const k_MB₅ = Q_mus/(V_mus*P_mus5)
const k_HB₅ = Q_hrt/(V_hrt*P_hrt5)
const k_BrB₅ = Q_br/(V_br*P_br5)

#Rate constant of MMA (/day)
const k_GBₘ = Q_gi/(V_gi*P_giM) #from GI tract to blood
const k_LiBₘ = Q_liv/(V_liv*P_livM) #from liver to blood
const k_KBₘ = Q_kid/(V_kid*P_kidM) #from kidney to blood
const k_SBₘ = Q_skin/(V_skin*P_skinM) #from skin to blood
const k_LuBₘ = Q_lung/(V_lung*P_lungM) #from lung to blood
const k_MBₘ = Q_mus/(V_mus*P_musM)
const k_HBₘ = Q_hrt/(V_hrt*P_hrtM)
const k_BrBₘ = Q_br/(V_br*P_brM)

#Rate constant of DMA (/day)
const k_GBₙ = Q_gi/(V_gi*P_giN) #from GI tract to blood
const k_LiBₙ = Q_liv/(V_liv*P_livN) #from liver to blood
const k_KBₙ = Q_kid/(V_kid*P_kidN) #from kidney to blood
const k_SBₙ = Q_skin/(V_skin*P_skinN) #from skin to blood
const k_LuBₙ = Q_lung/(V_lung*P_lungN) #from lung to blood
const k_MBₙ = Q_mus/(V_mus*P_musN)
const k_HBₙ = Q_hrt/(V_hrt*P_hrtN)
const k_BrBₙ = Q_br/(V_br*P_brN)

#Rate constant from blood to other organs for both organic and inorganic arsenic
const k_BG = Q_gi/V_blood#from Blood to GI tract
const k_BLi = Q_liv/V_blood #from blood to liver
const k_BK = Q_kid/V_blood #from blood to kidney
const k_BS = Q_skin/V_blood #from blood to skin
const k_BLu = Q_lung/V_blood #from blood to lung
const k_BM = Q_mus/V_blood
const k_BH = Q_hrt/V_blood
const k_BBr = Q_br/V_blood

#Reduction and Oxidation of inorganic arsenic
const R_Tis = 32.88 #As(V) reduction in tissues (/day)
const R_K = 42 #As(V) reduction in kidney (/day)
O_Tis = 43.92 #As(III) oxidation is tissues (/day)

#Metabolism constant for methylation of arsenic
const M_Li3M = 763 #from As3 to MMA in the liver (μmol/day)
const M_Li3N = 2880 #from As3 to DMA in the liver (μmol/day)
const M_LiMN = 950.4 #from MMA to DMA (μmol/day)
const M_K3M = 305.2 #from As3 to MMA in the kidney (μmol/day)
const M_K3N = 1152 #from As3 to DMA in the kidney (μmol/day)
const M_KMN = 380.16 #from MMA to DMA in the kidney (μmol/day)
const M = 3 #another metabolism constant (μmol/L), all of them have value 3 so it can be combined into one variable
#G₅ = 1.67




arsenic_once = @ode_def begin #arsenic model for one time ingestion
#Differential equations for As(III)
    dG₃ = - (k_GB₃ + eF + O_Tis) * G₃ + k_BG * B₃ + R_Tis * G₅ #Gut
    dLi₃ = - (k_LiB₃ + eB + O_Tis) * Li₃ + k_BLi * B₃ + R_Tis * Li₅ - (M_Li3M * (Li₃ / V_liv)) / (M + Li₃ / V_liv) - (M_Li3N * (Li₃ / V_liv)) / (M + Li₃ / V_liv) #Liver
    dK₃ = - (k_KB₃ + eUᵢ + O_Tis) * K₃ + k_BK * B₃ + R_K * K₅ - (M_K3M *( K₃ / V_kid))/(M + (K₃ / V_kid)) - (M_K3N * (K₃ / V_kid)) / (M + K₃ / V_kid) #Kidney
    dS₃ = - (k_SB₃ + O_Tis) * S₃ + k_BS * B₃ + R_Tis * S₅ #skin
    dLu₃ = - (k_LuB₃ + O_Tis) * Lu₃ + k_BLu * B₃ + R_Tis * Lu₅ #lung
    dB₃ = - (k_BG + k_BLi + k_BK + k_BS + k_BLu + k_BM + k_BH + k_BBr + O_Tis) * B₃ + k_GB₃ * G₃ + k_LiB₃ * Li₃ + k_KB₃ * K₃ + k_SB₃ * S₃ + k_LuB₃ * Lu₃ + k_MB₃ * M₃ + k_HB₃ * H₃ + k_BrB₃ * Br₃ + R_Tis * B₅ #blood
    dM₃ = - (k_MB₃ + O_Tis) * M₃ + k_BM * B₃ + R_Tis * M₅ #muscle
    dH₃ = - (k_HB₃ + O_Tis) * H₃ + k_BH * B₃ + R_Tis * H₅ #heart
    dBr₃ = - (k_BrB₃ + O_Tis) * Br₃ + k_BBr * B₃ + R_Tis * Br₅ #brain
    dU₃ = eUᵢ * K₃ #Urine
    dF₃ = eF * G₃ #Faeces
    dBil₃ = eB * Li₃ #Biliary

#Differential equations for As(V)
    dG₅ = - (k_GB₅ + eF + R_Tis) * G₅ + k_BG * B₅ + O_Tis * G₃
    dLi₅ = - (k_LiB₅ + eB + R_Tis) * Li₅ + k_BLi * B₅ + O_Tis * Li₃
    dK₅ = - (k_KB₅ + eUᵢ + R_K) * K₅ + k_BK * B₅ + O_Tis * K₃
    dS₅ = - (k_SB₅ + R_Tis) * S₅ + k_BS * B₅ + O_Tis * S₃
    dLu₅ = - (k_LuB₅ + R_Tis) * Lu₅ + k_BLu * B₅ + O_Tis * Lu₃
    dB₅ = - (k_BG + k_BLi + k_BK + k_BS + k_BLu + k_BM + k_BH + k_BBr + R_Tis) * B₅ + k_GB₅ * G₅ + k_LiB₅ * Li₅ + k_KB₅ * K₅ + k_SB₅ * S₅ + k_LuB₅ * Lu₅ + k_MB₅ * M₅ + k_HB₅ * H₅ + k_BrB₅ * Br₅ + O_Tis * B₃
    dM₅ = - (k_MB₅ + R_Tis) * M₅ + k_BM * B₅ + O_Tis * M₃
    dH₅ = - (k_HB₅ + R_Tis) * H₅ + k_BH * B₅ + O_Tis * H₃
    dBr₅ = - (k_BrB₅ + R_Tis) * Br₅ + k_BBr * B₅ + O_Tis * Br₃
    dU₅ = eUᵢ * K₅
    dF₅ = eF * G₅
    dBil₅ = eB * Li₅

#Differential equations for MMA
    dGₘ = - (k_GBₘ + eF) * Gₘ + k_BG * Bₘ
    dLiₘ = - (k_LiBₘ + eB) * Liₘ + k_BLi * Bₘ + (M_Li3M * (Li₃ / V_liv)) / (M + Li₃ / V_liv) - (M_LiMN * (Liₘ / V_liv)) / (M + Liₘ / V_liv)
    dKₘ = - (k_KBₘ + eUₘ) * Kₘ + k_BK * Bₘ + (M_K3M * (K₃ / V_kid)) / (M + K₃ / V_kid) - (M_KMN * (Kₘ / V_kid)) / (M + Kₘ / V_kid)
    dSₘ = - k_SBₘ * Sₘ + k_BS * Bₘ
    dLuₘ = - k_LuBₘ * Luₘ + k_BLu * Bₘ
    dBₘ = - (k_BG + k_BLi + k_BK + k_BS + k_BLu + k_BM + k_BH + k_BBr) * Bₘ + k_GBₘ * Gₘ + k_LiBₘ * Liₘ + k_KBₘ * Kₘ + k_SBₘ * Sₘ + k_LuBₘ * Luₘ + k_MBₘ * Mₘ + k_HBₘ * Hₘ + k_BrBₘ * Brₘ
    dMₘ = - k_MBₘ * Mₘ + k_BM * Bₘ
    dHₘ = - k_HBₘ * Hₘ + k_BH * Bₘ
    dBrₘ = - k_BrBₘ * Brₘ + k_BBr * Bₘ
    dUₘ = eUₘ * Kₘ
    dFₘ = eF * Gₘ
    dBilₘ = eB * Liₘ

#Differential equations for DMA
    dGₙ = - (k_GBₙ + eF) * Gₙ + k_BG * Bₙ
    dLiₙ = - (k_LiBₙ + eB) * Liₙ + k_BLi * Bₙ + (M_Li3N * (Li₃ / V_liv)) / (M + Li₃ / V_liv) + (M_LiMN * (Liₘ / V_liv)) / (M + Liₘ / V_liv)
    dKₙ = - (k_KBₙ + eUₙ) * Kₙ + k_BK * Bₙ + (M_K3N * (K₃ / V_kid)) / (M + K₃ / V_kid) + (M_KMN * (Kₘ / V_kid)) / (M + Kₘ / V_kid)
    dSₙ = - k_SBₙ * Sₙ + k_BS * Bₙ
    dLuₙ = - k_LuBₙ * Luₙ + k_BLu * Bₙ
    dBₙ = - (k_BG + k_BLi + k_BK + k_BS + k_BLu + k_BM + k_BH + k_BBr) * Bₙ + k_GBₙ * Gₙ + k_LiBₙ * Liₙ + k_KBₙ * Kₙ + k_SBₙ * Sₙ + k_LuBₙ * Luₙ + k_MBₙ * Mₙ + k_HBₙ * Hₙ + k_BrBₙ * Brₙ
    dMₙ = - k_MBₙ * Mₙ + k_BM * Bₙ
    dHₙ = - k_HBₙ * Hₙ + k_BH * Bₙ
    dBrₙ = - k_BrBₙ * Brₙ + k_BBr * Bₙ
    dUₙ = eUₙ * Kₙ
    dFₙ = eF * Gₙ
    dBilₙ = eB * Liₙ
end eUₘ


function arsenic_plot(intake,day,organ,K) #intake refers to the daily intake, organ refers to the order number of the organ in the differential equations (from 1 to 12)
    global initial = zeros(48)
    global initial[13] = intake
    global tspan = (0.0,1.0)
    global prob = ODEProblem(arsenic_once,initial,tspan,K)
    global sol = solve(prob)
    global x = sol.t
    global y = sol[organ,:] + sol[organ + 12,:] + sol[organ + 24,:] + sol[organ + 36,:]

    for i = 1:day
        global tspan = (i,i+1.0)
        global initial = sol[end]
        global initial[13] = initial[13] + intake
        global prob = ODEProblem(arsenic_once,initial,tspan,K)
        global sol = solve(prob)
        global x = append!(x,sol.t)
        global y = append!(y,sol[organ,:] + sol[organ + 12,:] + sol[organ + 24,:] + sol[organ + 36,:])
    end
    plot(x,y)
end

arsenic_plot(1e-3,100,6,eUₘ)

#initial = zeros(48)
#initial[13] = 0.1
#tspan = (0.0,1.0)
#prob = ODEProblem(arsenic_once,initial,tspan,1)
#sol = solve(prob)
#plot(sol)
