using PyCall
using DifferentialEquations
using ParameterizedFunctions
using Plots

I_tsc = 0.1     # Transcription rate
k_tsl = 0.1     # Translation rate
k_rna_deg = 0.1 # RNA degradation rate
k_pro_deg = 0.1 # Protein degradation rate
k_tr⁺ = 0.1     # Transportation rate extracellular -> intracellular
k_tr⁻ = 0.1     # Transportation rate intracellular -> extracellular
k_ac  = 0.1     # Accumulation rate

bacteria_uptake = @ode_def begin
    dmRNAᵀ = I_tsc - k_rna_deg * mRNAᵀ      # Transport mRNA
    dTP = mRNAᵀ * k_tsl - k_pro_deg * TP    # Transport protein

    dmRNAᴬ = I_tsc - k_rna_deg * mRNAᴬ      # Accumulation mRNA
    dAP = mRNAᴬ * k_tsl - k_pro_deg * AP    # Accumulation protein

    dEC = k_tr⁻ * IC * TP - k_tr⁺ * EC * TP                     # Extracellular HM
    dIC = -k_tr⁻ * IC * TP + k_tr⁺ * EC * TP - k_ac * AP * IC   # intracellular HM
    dIC_accumulated = k_ac * AP * IC - dIC_accumulated * k_pro_deg # Accumulated HM
end a


u0 = [0;0;0;0;1.0;0;0]
tspan = (0.0, 20.0)
prob = ODEProblem(bacteria_uptake,u0,tspan,1)
sol = solve(prob)
plotly()
plot(sol)
