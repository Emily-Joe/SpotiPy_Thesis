import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Define the model function
def model(phi, A, B, C):
    return A + B * np.sin(phi)**2 + C * np.sin(phi)**4

def fit_rotation_rate(phi, omega, omega_errors=None):
    # Specify initial values for A, B, and C
    initial_guess = [14.2, -2.5, 0]

    # Specify bounds for parameters
    parameter_bounds = ([-np.inf, -np.inf, -np.inf], [np.inf, 0, np.inf])

    if omega_errors is None:
        popt, pcov = curve_fit(model, phi, omega, p0=initial_guess, bounds=parameter_bounds)
    else:
        popt, pcov = curve_fit(model, phi, omega, sigma=omega_errors, absolute_sigma=True, p0=initial_guess, bounds=parameter_bounds)

    # Ensure parameter B is negative
    if popt[1] > 0:
        popt[1] *= -1

    A_fit, B_fit, C_fit = popt

    # Calculate R-squared value
    residuals = omega - model(phi, A_fit, B_fit, C_fit)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((omega - np.mean(omega))**2)
    r_squared = 1 - (ss_res / ss_tot)

    # Calculate uncertainties on parameters A, B, and C from the covariance matrix
    A_uncertainty, B_uncertainty, C_uncertainty = np.sqrt(np.diag(pcov))

    return A_fit, B_fit, C_fit, A_uncertainty, B_uncertainty, C_uncertainty, r_squared

def model2(phi, A, B): # without the C parameter
    return A + B * np.sin(phi)**2

def fit_rotation_rate2(phi, omega, omega_errors=None): # without the C parameter
    # Specify initial values for A, B, and C
    initial_guess = [14.2, -2.5]

    # Specify bounds for parameters
    parameter_bounds = ([-np.inf, -np.inf], [np.inf, 0])

    if omega_errors is None:
        popt, pcov = curve_fit(model2, phi, omega, p0=initial_guess, bounds=parameter_bounds)
    else:
        popt, pcov = curve_fit(model2, phi, omega, sigma=omega_errors, absolute_sigma=True, p0=initial_guess, bounds=parameter_bounds)

    # Ensure parameter B is negative
    if popt[1] > 0:
        popt[1] *= -1

    A_fit2, B_fit2 = popt

    # Calculate R-squared value
    residuals = omega - model2(phi, A_fit2, B_fit2)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((omega - np.mean(omega))**2)
    r_squared2 = 1 - (ss_res / ss_tot)

    # Calculate uncertainties on parameters A, B, and C from the covariance matrix
    A_uncertainty2, B_uncertainty2 = np.sqrt(np.diag(pcov))

    return A_fit2, B_fit2, A_uncertainty2, B_uncertainty2, r_squared2

# folder containing all the data:
mydir = '/work2/loessnitz/'

dt = 12 # do we still need this here?

# set up my fit model
Loessnitz = model(x*np.pi/180, A_fit, B_fit, C_fit)

NOAA_list = np.loadtxt("all_NOAA.txt", dtype=str, comments="#", unpack=False) # txt list containing  all the NOAAs of all the spots that have been downloaded and ran at least once, which should be considered in this plot
length=len(NOAA_list)

# make a file containing all the data from all the spots: (delete later, if done once):
all_data = open("full_data_sunspot_rotation.txt", "w+")
all_data.write('#NOAA,SCORE,CADANCE,LATITUDE,Rotational_Velocity,Standard_Deviation,Mean_Standard_Deviation\n')

# set up arrays to save data in:
phi_data1   = []
omega_data1 = []
omega_errors1= []
weights1    = []

for i in range(length):
    NOAA = NOAA_list[i]
    path = 'NOAA_' + str(NOAA) + '_dt_' + str(dt) + 'h'
    NOAA_result = np.loadtxt(path+'/'+str(NOAA)+"_results.txt", dtype=str, comments="#", delimiter=" ", unpack=False)   # read out results for each spot in the NOAA-list

    # read out the following informatio:
    score = int(NOAA_result[1][1])
    cadance = int(NOAA_result[2][1])
    latitude = int(NOAA_result[3][1])
    c_omega = float(NOAA_result[7][1])
    c_omega_std = float(NOAA_result[8][1])
    c_mean_omega_std = float(NOAA_result[9][1])
    e_omega = float(NOAA_result[10][1])
    e_omega_std = float(NOAA_result[10][1])
    e_mean_omega_std = float(NOAA_result[12][1])

    # save to arrays:
    phi_data1.append(latitude)
    omega_data1.append(e_omega)
    omega_errors1.append(e_mean_omega_std)
    weights1.append(score)

    #write to final output file (delete later):
    all_data.write('{0},{1},{2},{3},{4},{5},{6}\n'.format(NOAA, score, cadance, latitude, e_omega, e_omega_std, e_mean_omega_std))

#delete later:
all_data.close()

# make sure data is an array
phi_data = np.array(phi_data1)
omega_data= np.array(omega_data1)
omega_errors=np.array(omega_errors1)
weights=np.array(weights1)

# Sort the data based on phi (latitude)
sorted_indices = np.argsort(phi_data)
phi_data_sorted = phi_data[sorted_indices]
omega_data_sorted = omega_data[sorted_indices]
omega_errors_sorted = omega_errors[sorted_indices]
weights_sorted = weights[sorted_indices]
weigthed_errors = omega_errors_sorted/(weights_sorted+0.1) # weigh the errors based on the score given (0=worst; 3=best) low-score spots should have bigger errors, but need to add +0.1 so that we dont divide by zero

# Fit the sorted data with errors and obtain uncertainties
A_fit, B_fit, C_fit, A_uncertainty, B_uncertainty, C_uncertainty, r_squared = fit_rotation_rate(phi_data_sorted*np.pi/180, omega_data_sorted, weigthed_errors)
A_fit2, B_fit2, A_uncertainty2, B_uncertainty2, r_squared2 = fit_rotation_rate2(phi_data_sorted*np.pi/180, omega_data_sorted, weigthed_errors)

# Print the fitted parameters, uncertainties, and R-squared value
print("A:", A_fit)
print("B:", B_fit)
print("C:", C_fit)
print("A Uncertainty:", A_uncertainty)
print("B Uncertainty:", B_uncertainty)
print("C Uncertainty:", C_uncertainty)
print("R-squared:", r_squared)

# write fit results into output file
fit_parameters = open("rotation_fit_results.txt", "w+")
fit_parameters.write('#Fit Parameters for all alpha sunspots\n#Model:A+B*np.sin(phi)**2+C*np.sin(phi)**4 \nA: {0}\nB: {1}\nC: {2}\nA_Uncertainty: {3}\nB_Uncertainty: {4}\nC_Uncertainty: {5}\nR-squared: {6}\n\n#Model:A+B*np.sin(phi)**2 \nA: {7}\nB: {8}\nA_Uncertainty: {9}\nB_Uncertainty: {10}\nR-squared: {11}'.format(A_fit, B_fit, C_fit, A_uncertainty, B_uncertainty, C_uncertainty, r_squared, A_fit2, B_fit2, A_uncertainty2, B_uncertainty2, r_squared2))
fit_parameters.close()


# seperate the data by thier weight (this is to plot them with different colors)
# class 0:
phi_data0 = phi_data_sorted[np.where(weights_sorted == 0)]
ome_data0 = omega_data_sorted[np.where(weights_sorted == 0)]
err_data0 = omega_errors_sorted[np.where(weights_sorted == 0)]
# class 1:
phi_data1 = phi_data_sorted[np.where(weights_sorted == 1)]
ome_data1 = omega_data_sorted[np.where(weights_sorted == 1)]
err_data1 = omega_errors_sorted[np.where(weights_sorted == 1)]
# class 2:
phi_data2 = phi_data_sorted[np.where(weights_sorted == 2)]
ome_data2 = omega_data_sorted[np.where(weights_sorted == 2)]
err_data2 = omega_errors_sorted[np.where(weights_sorted == 2)]
# class 3:
phi_data3 = phi_data_sorted[np.where(weights_sorted == 3)]
ome_data3 = omega_data_sorted[np.where(weights_sorted == 3)]
err_data3 = omega_errors_sorted[np.where(weights_sorted == 3)]

# set up x-axis:
x=np.linspace(-40,40,100)

# reference models:
Balthasar = 14.551 + (-2.87)* (np.sin(x*np.pi/180)**2)
Snodgrass = 14.050 + (-1.492)* (np.sin(x*np.pi/180)**2) + (-2.606)* (np.sin(x*np.pi/180)**4)
Loessnitz = model(x*np.pi/180, A_fit, B_fit, C_fit)    # my models
Loessnitz2= model2(x*np.pi/180, A_fit2, B_fit2)        # my models

# Plot the data and the model curve using the fitted parameters
plt.figure(figsize=(12, 6))
plt.errorbar(phi_data0, ome_data0, yerr=err_data0, c='gray',  label='Data class 0', capsize=5, fmt ='^')        #class 0
plt.errorbar(phi_data1, ome_data1, yerr=err_data1, c='orange',  label='Data class 1', capsize=5, fmt ='s')      #class 1
plt.errorbar(phi_data2, ome_data2, yerr=err_data2, c='forestgreen',  label='Data class 2', capsize=5, fmt ='p') #class 2
plt.errorbar(phi_data3, ome_data3, yerr=err_data3, c='blue',  label='Data class 3', capsize=5, fmt ='o')        #class 3

# plot/ draw models:
plt.plot(x, Balthasar, 'dimgray', linestyle="dotted", label='Balthasar (all sunspots)' )
plt.plot(x, Snodgrass, 'dimgray', linestyle="dashed", label='Snodgrass (quiet sun)' )
plt.plot(x, Loessnitz, 'k', linestyle="solid", label=r'Lößnitz 1 (Alpha Spots)' )
plt.plot(x, Loessnitz2, 'k', linestyle="dashdot", label=r'Lößnitz 2 (Alpha Spots)')


plt.title('rotational velocity distribution (from 09.2013 to 01.2024)')
plt.xlabel("latitude [deg]")
plt.ylabel("rotational velocity [deg/day]")
plt.ylim([13, 15])

ax = plt.gca()
ax.grid(True)

leg = plt.legend()
plt.savefig("full_rotation_plot_color.pdf", dpi=500)
plt.show()


