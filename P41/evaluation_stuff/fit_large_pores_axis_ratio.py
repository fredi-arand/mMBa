import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
import math

home_folder = os.getenv("HOME")
file_path = home_folder + "/build/P41/anisotropy_comparison_7.0-large_pores_downsampled.txt"
file_input = open(file_path, "r")

lines = file_input.readlines()
lines.pop(0)

x_values = []
y_values = []

for line in lines:
    bin_low, bin_up, sig32, vx, vy, vz, sig12, vx1, vy1, vz1, vx2, vy2, vz2 = ( float(ele) for ele in line.split() )
    x_values.append( (math.log(2.0*bin_low)+math.log(2.0*bin_up))/2.0 )
    y_values.append( sig32 )

plt.plot(x_values, y_values, 'b.')

def fit_function(x, A, B):
    return A*x+B

popt, pcov = curve_fit(fit_function, x_values, y_values)

plt.plot(np.linspace(0.0,math.log(300.0)), fit_function(np.linspace(0.0,math.log(300.0)),*popt), 'r-')
plt.show()

print("fitted parameters:")
print("A = " + str(popt[0]) + "   B = " + str(popt[1]) )