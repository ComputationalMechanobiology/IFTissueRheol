using RHEOS
using PyPlot
using Plots

###############################################################################################################################
###################################### Fit experimental data to Fractional Solid ##############################################
###############################################################################################################################

data = importcsv("C:/Users/eleni/OneDrive/Documents/MRes_Sensor_CDT/Mini_Project/cycling/22-03-01_m4/cyclingm4300.csv", t_col=1, ϵ_col=6, σ_col=3)

fract_solid = modelfit(data, FractSolid, strain_imposed)
fractsolid_predict = modelpredict(extract(data, strain_only), fract_solid)

pygui(true)
fig, ax = subplots()

ax.plot(data.t, data.ϵ, "--", label="Strain (Data)")
ax.plot(data.t, data.σ, ".", label="Stress (Data)")
ax.plot(fractsolid_predict.t, fractsolid_predict.σ, label="Stress (Predicted)")

ax.set_xlabel("Time (s)")
ax.set_ylabel("Stress, Strain")
ax.legend()


