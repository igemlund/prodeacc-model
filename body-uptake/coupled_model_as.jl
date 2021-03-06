# Based on Dede et al 2018
using Plots
using StatsPlots
using DifferentialEquations
using ParameterizedFunctions
using IterableTables
include("../append_solution.jl")
pyplot()

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
#partition coefficient is the ratio of concentrations of a compound in a mixture of two immiscible phases at equilibrium
P_giM = 2.2 #GI tract/blood partition coefficient for MMA ((GI/V_gi)/(B/V_blood))
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
eUᵢ = 100.8 #Urinary excretion rate for As3 and As5 (/day)
eUₘ = 432 #Urinary excretion rate for MMA (/day)
eUₙ = 187.2 #Urinary excretion rate for DMA (/day)


#Rate constant of As(III) (/day)
k_GB₃ = Q_gi/(V_gi*P_gi3) #from GI tract to blood
k_LiB₃ = Q_liv/(V_liv*P_liv3) #from liver to blood
k_KB₃ = Q_kid/(V_kid*P_kid3) #from kidney to blood
k_SB₃ = Q_skin/(V_skin*P_skin3) #from skin to blood
k_LuB₃ = Q_lung/(V_lung*P_lung3) #from lung to blood
k_MB₃ = Q_mus/(V_mus*P_mus3)
k_HB₃ = Q_hrt/(V_hrt*P_hrt3)
k_BrB₃ = Q_br/(V_br*P_br3)

#Rate constant of As(V) (/day)
k_GB₅ = Q_gi/(V_gi*P_gi5) #from GI tract to blood
k_LiB₅ = Q_liv/(V_liv*P_liv5) #from liver to blood
k_KB₅ = Q_kid/(V_kid*P_kid5) #from kidney to blood
k_SB₅ = Q_skin/(V_skin*P_skin5) #from skin to blood
k_LuB₅ = Q_lung/(V_lung*P_lung5) #from lung to blood
k_MB₅ = Q_mus/(V_mus*P_mus5)
k_HB₅ = Q_hrt/(V_hrt*P_hrt5)
k_BrB₅ = Q_br/(V_br*P_br5)

#Rate constant of MMA (/day)
k_GBₘ = Q_gi/(V_gi*P_giM) #from GI tract to blood
k_LiBₘ = Q_liv/(V_liv*P_livM) #from liver to blood
k_KBₘ = Q_kid/(V_kid*P_kidM) #from kidney to blood
k_SBₘ = Q_skin/(V_skin*P_skinM) #from skin to blood
k_LuBₘ = Q_lung/(V_lung*P_lungM) #from lung to blood
k_MBₘ = Q_mus/(V_mus*P_musM)
k_HBₘ = Q_hrt/(V_hrt*P_hrtM)
k_BrBₘ = Q_br/(V_br*P_brM)

#Rate constant of DMA (/day)
k_GBₙ = Q_gi/(V_gi*P_giN) #from GI tract to blood
k_LiBₙ = Q_liv/(V_liv*P_livN) #from liver to blood
k_KBₙ = Q_kid/(V_kid*P_kidN) #from kidney to blood
k_SBₙ = Q_skin/(V_skin*P_skinN) #from skin to blood
k_LuBₙ = Q_lung/(V_lung*P_lungN) #from lung to blood
k_MBₙ = Q_mus/(V_mus*P_musN)
k_HBₙ = Q_hrt/(V_hrt*P_hrtN)
k_BrBₙ = Q_br/(V_br*P_brN)

#Rate constant from blood to other organs for both organic and inorganic arsenic
k_BG = Q_gi/V_blood#from Blood to GI tract
k_BLi = Q_liv/V_blood #from blood to liver
k_BK = Q_kid/V_blood #from blood to kidney
k_BS = Q_skin/V_blood #from blood to skin
k_BLu = Q_lung/V_blood #from blood to lung
k_BM = Q_mus/V_blood
k_BH = Q_hrt/V_blood
k_BBr = Q_br/V_blood

#Reduction and Oxidation of inorganic arsenic
R_Tis = 32.88 #As(V) reduction in tissues (/day)
R_K = 42 #As(V) reduction in kidney (/day)
O_Tis = 43.92 #As(III) oxidation is tissues (/day)

Av = 6.022e23
V_cell = 8e-16
V_gut = 1.2

mcg2ions(mcg, molarmass) = mcg*1e-6/molarmass*Av
mcg2mcm(mcg,molarmass) = mcg/molarmass
ions2mcg(ions, molarmass) = ions/Av*molarmass*1e6
mcm2ions(mcm) = mcm*1e-6*Av

#Metabolism constant for methylation of arsenic
M_Li3M = 763 #from As3 to MMA in the liver (μmol/day)
M_Li3N = 2880 #from As3 to DMA in the liver (μmol/day)
M_LiMN = 950.4 #from MMA to DMA (μmol/day)
M_K3M = 305.2 #from As3 to MMA in the kidney (μmol/day)
M_K3N = 1152 #from As3 to DMA in the kidney (μmol/day)
M_KMN = 380.16 #from MMA to DMA in the kidney (μmol/day)
M = 3 #another metabolism constant (μmol/L), all of them have value 3 so it can be combined into one variable


k_growth = log(2)/40*24        # Growth rate [/day], double every 40 hours in human gut.
N_max = 1e13                 # Max number of cells
k_tsl = 0.075*1e-9*Av*V_cell*60*24   # Translation rate [/day]
tsc =  15*60*24                       # Transcription rate [/day]
k_rna_deg = log(2)/5*60*24                # RNA degradation rate [/day]
k_pro_deg = log(2)/(2/3)*24                   # Protein degradation rate [/day]
#V_max = 3.5676e-1*24                  # unknown [ion/day/TP] (there's no difference for this V_max but there is when V_max is increased by say 50 times)
V_max = 60.1*24                  # unknown [ion/day/TP] (there's no difference for this V_max but there is when V_max is increased by say 50 times)
k_m = 26e-2              # unknown [ion/m^3]
k_CA = 0.1                        # unknown Accumulation protein capture rate
k_CG = 0.01                       # unknonw Ions Cytoplasm --> Gut
N0 = 1e11                       # Number of cells at start
G0 = mcg2mcm(15,75)            # Gut start amount


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
end a


arsenic_coupled = @ode_def begin #arsenic model coupled with colony model
#Differential equations for As(III)
    dG₃ = - (k_GB₃ + eF + O_Tis) * G₃ + k_BG * B₃ + R_Tis * G₅ - (V_max/(Av*1e-6)) / (1 + (k_m/(Av*1e-6)) / (G₃/V_gut)) * TP * N #Gut
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
    dG₅ = - (k_GB₅ + eF + R_Tis) * G₅ + k_BG * B₅ + O_Tis * G₃ - (V_max/(Av*1e-6)) / (1 + (k_m/(Av*1e-6)) / (G₃/V_gut)) * TP * N
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

    # Colony
    dN = k_growth * N - k_growth * (N^2) / N_max    # Number of bacteria in gut

    # Single Bacteria
    dTP = (k_tsl * mRNA
        - dN/N * TP
        - k_pro_deg * TP)
    dAP_free = (k_tsl * mRNA           # Number of free accumulators
        - dN/N * AP_free
        - k_CA * AP_free * C
        - k_pro_deg * AP_free)
    dAP_occupied = (k_CA * AP_free * C  # Number of occupied accumulators
        - dN/N * AP_occupied
        - k_pro_deg * AP_occupied)
    dC = (V_max/(1+k_m/(G₃/V_gut))*TP + V_max/(1+k_m/(G₅/V_gut))*TP# Number of ions in cytoplasm
        - k_CA * C * AP_free
        + k_pro_deg * AP_occupied
        - dN/N * C)
    dmRNA = tsc - k_rna_deg * mRNA
end a


function arsenic_body_uptake(intake,days,init_bac,coupled=false)
    tspan = (0.0,1.0)
    if coupled
         model = arsenic_coupled
         u0 = zeros(54)
         u0[13] = intake
         u0[49] = init_bac
    else
        model = arsenic_once
        u0 = zeros(48)
        u0[13] = intake
    end
    sol = solve(ODEProblem(model, u0, tspan), dense=false)
    for t in 2:days
        tspan = tspan .+ 1
        u0 = sol.u[end]
        u0[13] += intake
        if coupled u0[49] = init_bac end
        sol1 = sol
        sol2 = solve(ODEProblem(model, u0, tspan), dense=false)
        sol = sol_append(sol1, sol2)
    end
    sol
end

sol1 = arsenic_body_uptake(G0,30,N0) #without colony
sol2 =  arsenic_body_uptake(G0,30,N0,true) #with colony
plot(sol2[13,:])

sol1_end = sol1[end]
sol2_end = sol2[end]

index_liver = [2, 14, 26, 38]
index_kidneys = [2, 14, 26, 38].+1
index_skin = [2, 14, 26, 38].+2
index_lungs = [2, 14, 26, 38].+3

levels1 = [
    sum(sol1_end[index_liver]),
    sum(sol1_end[index_kidneys]),
    sum(sol1_end[index_skin]),
    sum(sol1_end[index_lungs]),
    ]

levels2 = [
    sum(sol2_end[index_liver]),
    sum(sol2_end[index_kidneys]),
    sum(sol2_end[index_skin]),
    sum(sol2_end[index_lungs]),
    ]

1-sum(levels2./levels1)/4

plotlyjs()
data = []
append!(data, levels1)
append!(data, levels2)
data = float.(data)
g = repeat(["Without Probiotic", "With Probiotic"], inner = 4)
nam = [
    "Liver", "Kidneys", "Skin", "Lungs",
    "Liver", "Kidneys", "Skin", "Lungs"
    ]
groupedbar(nam, data,  group = g, ylabel="Micro mol of As ions")
savefig("As_compare_30d_15mcmol.png")
savefig("As_compare_30d_15mcmol.svg")


#Dietary arsenic intakes estimated from various countries range from less than 10 µg/day to 200 µg/day (WHO)
