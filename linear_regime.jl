using RHEOS
using PyPlot
using Plots
using DataFrames

###############################################################################################################################
###################################### Fit experimental data to Fractional Solid ##############################################
###############################################################################################################################

data = importcsv("C:/Users/eleni/OneDrive/Documents/MRes_Sensor_CDT/Mini_Project/cycling/22-03-01_m4/single_cycle.csv", t_col=7, ϵ_col=6, σ_col=9)

fract_solid = modelfit(data, FractSolid, strain_imposed, lo=(η = 0.0001, cᵦ = 0.0001, β = 0.0001, k = 0.0001), hi=(η = 10000, cᵦ = 10000, β = 10000, k = 10000), verbose=true)
fractsolid_predict = modelpredict(extract(data, strain_only), fract_solid)
#sls_zener = modelfit(data, SLS_Zener, strain_imposed, lo=(η = 0.0001, kᵦ = 0.0001, kᵧ = 0.0001), hi=(η = 10000, kᵦ = 10000, kᵧ = 10000), verbose=true)
#slszener_predict = modelpredict(extract(data, strain_only), sls_zener)

pygui(true)
fig, ax = subplots()

ax.plot(data.t, data.ϵ, "--", label="Strain (Data)")
ax.plot(data.t, data.σ, ".", label="Stress (Data)")
ax.plot(fractsolid_predict.t, fractsolid_predict.σ, label="Stress (Predicted)")

ax.set_xlabel("Time (s)")
ax.set_ylabel("Stress, Strain")
ax.legend()

###############################################################################################################################
################################# Down-sample experimental data and re-fit to Fractional Solid ################################
###############################################################################################################################

data_resampled = resample(data, scale = 0.1)
println(length(data_resampled.t)) #expected 113

fract_solid_resampled = modelfit(data_resampled, FractSolid, strain_imposed, lo=(η = 0.0001, cᵦ = 0.0001, β = 0.0001, k = 0.0001), hi=(η = 10000, cᵦ = 10000, β = 10000, k = 10000), verbose=true)
fractsolid_predict_resampled = modelpredict(extract(data_resampled, strain_only), fract_solid_resampled)

pygui(true)
fig, ax = subplots()

ax.plot(data_resampled.t, data_resampled.ϵ, "--", label="Strain (Data)")
ax.plot(data_resampled.t, data_resampled.σ, ".", label="Stress (Data)")
ax.plot(fractsolid_predict_resampled.t, fractsolid_predict_resampled.σ, label="Stress (Predicted)")

ax.set_title("Fractional Solid Fit")
ax.set_xlabel("Time (s)")
ax.set_ylabel("Stress, Strain")
ax.legend()

###############################################################################################################################
####################################################### All cycles ############################################################
###############################################################################################################################

dc_data = importcsv("C:/Users/eleni/OneDrive/Documents/MRes_Sensor_CDT/Mini_Project/cycling/22-03-01_m4/all_cycles.csv", t_col=7, ϵ_col=6, σ_col=9)
dc_data_resampled = resample(dc_data, scale = 0.2)

fract_solid_dc = modelfit(dc_data_resampled, FractSolid, strain_imposed, lo=(η = 0.0001, cᵦ = 0.0001, β = 0.0001, k = 0.0001), hi=(η = 10000, cᵦ = 10000, β = 10000, k = 10000), verbose=true, rel_tol=1e-6)
fractsolid_predict_dc = modelpredict(extract(dc_data_resampled, strain_only), fract_solid_dc)

pygui(true)
fig, ax = subplots()

ax.plot(dc_data_resampled.t, dc_data_resampled.ϵ, "--", label="Strain (Data)")
ax.plot(dc_data_resampled.t, dc_data_resampled.σ, ".", label="Stress (Data)")
ax.plot(fractsolid_predict_dc.t, fractsolid_predict_dc.σ, label="Stress (Predicted)")

ax.set_title("Fractional Solid Fit")
ax.set_xlabel("Time (s)")
ax.set_ylabel("Stress, Strain")
ax.legend()

#Time: 9207.7284172 s, Why: XTOL_REACHED, Parameters: [2.125555869960186, 1.633632673867637, 0.0001, 0.2760957120720441], Error: 0.0031826757770975524


###############################################################################################################################
############################## Fit undersampled experimental data to Standard Linear Solid ####################################
###############################################################################################################################

dc_data = importcsv("C:/Users/eleni/OneDrive/Documents/MRes_Sensor_CDT/Mini_Project/cycling/22-03-01_m4/all_cycles.csv", t_col=7, ϵ_col=6, σ_col=9)
dc_data_resampled = resample(dc_data, scale = 0.2)

sls_resampled = modelfit(dc_data_resampled, SLS_Zener, strain_imposed,  verbose=true, rel_tol=1e-6)
sls_predict_resampled = modelpredict(extract(dc_data_resampled, strain_only), sls_resampled)

pygui(true)
fig, ax = subplots()

ax.plot(dc_data_resampled.t, dc_data_resampled.ϵ, "--", label="Strain (Data)")
ax.plot(dc_data_resampled.t, dc_data_resampled.σ, ".", label="Stress (Data)")
ax.plot(sls_predict_resampled.t, sls_predict_resampled.σ, label="Stress (Predicted)")

ax.set_title("Standard Linear Solid Fit")
ax.set_xlabel("Time (s)")
ax.set_ylabel("Stress, Strain")
ax.legend()

#Result: Time: 3.1867916 s, Why: XTOL_REACHED, Parameters: [2.113667257378493, 1.5586247109007703, 0.2758070503146254], Error: 0.003193660629536194

#optimiser_list = [:LN_COBYLA,:LN_SBPLX,:LN_BOBYQA]
#Can add optmethod = "" in modelfit to change optimised method

