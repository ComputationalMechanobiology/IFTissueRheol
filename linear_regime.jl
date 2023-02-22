using RHEOS
using PyPlot
using Plots

###############################################################################################################################
###################################### Fit experimental data to Fractional Solid ##############################################
###############################################################################################################################



data = importcsv("C:/Users/eleni/OneDrive/Documents/MRes_Sensor_CDT/Mini_Project/cycling/22-03-01_m4/single_cycle.csv", t_col=7, ϵ_col=6, σ_col=9)

println(data.t)

fract_solid = modelfit(data, FractSolid, strain_imposed, lo=(η = 0.0001, cᵦ = 0.0001, β = 0.0001, k = 0.0001), hi=(η = 10000, cᵦ = 10000, β = 10000, k = 10000), verbose=true)
fractsolid_predict = modelpredict(extract(data, strain_only), fract_solid)

pygui(true)
fig, ax = subplots()

ax.plot(data.t, data.ϵ, "--", label="Strain (Data)")
ax.plot(data.t, data.σ, ".", label="Stress (Data)")
ax.plot(fractsolid_predict.t, fractsolid_predict.σ, label="Stress (Predicted)")

ax.set_xlabel("Time (s)")
ax.set_ylabel("Stress, Strain")
ax.legend()


