using PyPlot
using DelimitedFiles
using LsqFit
using ImageFiltering: imfilter, Kernel



function main(folder, file_name)

    print("CALIBRATION running... \n")

    source = string(folder, file_name);

    print("Sample = ", source, "\n")

    file = readdlm(source, '\t'; header=true)[1]
    file_calibration = readdlm(string(folder, "calibration.txt"), ' ');

    print("String values in file = ", file[findall( x -> isa(x,SubString{String}), file[:,3]),1], "\n")
    

    L_wire = file_calibration[1,1]; # flexible wire length (m)
    L = file_calibration[1,2]; # monolayer length (m)
    w = (file_calibration[2,1] + 2*file_calibration[2,2] + file_calibration[2,3])/4; # monolayer width (m)
    t = 10e-6 # thickness (m)

    print(L_wire, " ", file_calibration[3:end,2])
    alpha = calibration(L_wire, file_calibration[3:end,2], file_calibration[3:end,1]);   # Each step is 100um
    force = file[:,3].* alpha;

    strain = (file[1:end,2] .- file[1,2]).*1e-3 ./L; # strain ( not in %)
    stress = force ./ (w*t) # stress (N/m)    
    time = (file[1:end,1] .- file[1,1])

    A = [file[:,1] file[:,2] force.-force[1]] # time (s)   position (mm)    force(N)-> first point shifted to zero
    B = [time strain stress.-stress[1]]
    
    (fig, ax) = subplots(1, 2, figsize=(10,5))
    fig.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.95, wspace = 0.3)

    ax[1].plot(B[:,2], B[:,3], "blue")
    ax[1].set_xlabel("strain");
    ax[1].set_ylabel("stress (N/m^2) "); 

    ax[2].plot(B[:,1], B[:,3], "blue");
    ax[2].set_xlabel("time");
    ax[2].set_ylabel("stress (N/m^2) "); 

    writedlm(string(folder,file_name[1:end-4], "_calibrated.txt"), A, '\t')
    writedlm(string(folder,file_name[1:end-4], "_tss.txt"), B, '\t') # _time strain stress

    print("alpha = ", alpha, "\n")
    print("Wire length (m) = ", L_wire , "\n")
    print("Monolayer width - w (m) = ", w, "\n")
    print("Monolayer length - L (m) = ", L, "\n")

    print("Max force (N) = ", maximum(A[:,3]), "\n")
    print("Max stress (N/m^2) = ", maximum(B[:,3]))


end




function calibration(L,V,dist)
    # alpha times the voltage is the force the monolayer experience
        # Wire information
        # L # length of the flexible wirein meters
         E = 75e9; # Pascals
         r = 0.05e-3; # in meters
         I = (pi*r^4)/4;
         k = (3*E*I)/L^3;
    
    (fig1, ax1) = subplots(1, 1, figsize=(5,5))

    ax1.plot(V,k.*dist, "o", color = "blue")
    ax1.set_xlabel("Voltage (V)");
    ax1.set_ylabel("k*dist (Pa*m^2)"); 

    @. model(x, p) = p[1].* x + p[2]
    p0 = [0.5, 0.5]
    fit = curve_fit(model, V, k.*dist, p0)

    ax1.plot(V, fit.param[1].*V .+ fit.param[2],"red")
    savefig("Calibration.png")


    return fit.param[1]

    

end

#-------------------------------------------------------------


folder = "C:/Users/eleni/OneDrive/Documents/MRes_Sensor_CDT/Mini_Project/cycling/22-03-01_m4/"
file_name = "cycling_mod.txt"
main(folder, file_name)