#using PyPlot
using Images, FileIO
using DelimitedFiles
using Plots
using PyCall, Conda
#pyplot()

#Conda.add("opencv") #install opencv from conda
cv2 = pyimport("cv2") #import the cv2 module from conda and define as cv2

function calculate_centroid(image_process)
    # convert image to grayscale image
    gray_image = cv2.cvtColor(image_process, cv2.COLOR_BGR2GRAY);
    # convert the grayscale image to binary image
    ret,thresh = cv2.threshold(gray_image,127,255,0);
    # calculate moments of binary image
    M = cv2.moments(thresh);

    # calculate x,y coordinate of center
    cX = Integer(round(M["m10"] / M["m00"]));
    cY = Integer(round(M["m01"] / M["m00"]));

    print("c_x = ", cX, " c_y = ", cY, "\n");

    return cX, cY

end

function calculate_displacement(folder_path, folder_calibration)
    print("Calculate centroid BEFORE... \n")
    path_before = string(folder_path, "/before.tif")
    img_before = cv2.imread(path_before);

    cX_before, cY_before = calculate_centroid(img_before)
    img_plot = Images.load(path_before);
    plt_before = Plots.plot(img_plot)
    plt_before = Plots.plot!([0, cX_before], [0, cY_before], markersize = 8)
    display(plt_before)
    savefig(plt_before, string(folder_path,"/before_centroid.png"))

    print("Calculate centroid AFTER... \n")
    path_after = string(folder_path, "/after.tif")
    img_after = cv2.imread(path_after);


    cX_after, cY_after = calculate_centroid(img_after)
    img_plot = Images.load(path_after);
    plt_after = Plots.plot(img_plot)
    plt_after = Plots.plot!([0, cX_after], [0, cY_after], markersize = 8)
    display(plt_after)
    savefig(plt_after, string(folder_path,"/after_centroid.png"))

    displacement = cX_after-cX_before;

    print("Displacement X (pixels)= ", displacement, "\n")
    print("Displacement X (m)= ", displacement/215 * 1e-3, "\n")

    pre_stress = calculate_pre_stress(folder_calibration, displacement)
    writedlm(string(folder_path, "/pre_stress.txt"), pre_stress, '\t')

    print("Pre stress (Pa)= ", pre_stress)
end

function calculate_pre_stress(folder_calibration, displacement_pixels)

    
    file_calibration = readdlm(string(folder_calibration, "calibration.txt"), ' ');
    L_wire = file_calibration[1,1]; # flexible wire length (m)
    L = file_calibration[1,2]; # monolayer length (m)
    w = (file_calibration[2,1] + 2*file_calibration[2,2] + file_calibration[2,3])/4; # monolayer width (m)


    # alpha times the voltage is the force the monolayer experience
    # Wire information
    # L # length of the flexible wirein meters
    E = 75e9; # Pascals
    r = 0.05e-3; # in meters
    I = (pi*r^4)/4;
    k = (3*E*I)/L_wire^3;  #Pa * m
    t = 10e-6 #m

    pre_stress = -k * (displacement_pixels/215 * 1e-3) / (w*t)   # Pa m^2 / m^2

    return pre_stress

end

folder_root = "C:/Users/eleni/OneDrive/Documents/MRes_Sensor_CDT/Mini_Project/cycling/test"

folder_path = string(folder_root)
folder_calibration = string(folder_root, "/")

calculate_displacement(folder_path, folder_calibration)