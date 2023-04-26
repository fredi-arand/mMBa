import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
import math

def gauss_fit(x, A, mu, sigma):
    return A/np.sqrt(2*math.pi*sigma**2)*np.exp(-(x-mu)**2/(2*sigma**2))

def boundary_correction(x):
    return ((325-x)/325)**3

home_folder = os.getenv("HOME")
file_path = home_folder + "/build/P41/pores_smaller_7um.txt"
file_input = open(file_path, "r")

lines = file_input.readlines()
lines.pop(0)

x_values = []
y_values = []
y2_values = []

for line in lines:
    x_low, x_up, y = ( float(ele) for ele in line.split() )
    if x_low <1.0:
        continue
    if x_up > 6.0:
        continue

    x_center = (math.log(x_low)+math.log(x_up))/2.0
    
    x_values.append( x_center )
    y_values.append( y/boundary_correction(x_center) )
    y2_values.append(y)

plt.plot(x_values, y_values, 'b-')

popt, pcov = curve_fit(gauss_fit, x_values, y_values)

plt.plot(x_values, gauss_fit(x_values,*popt), 'r-')
plt.show()

print("fitted parameters:")
print("A = " + str(popt[0]) + "   mu = " + str(popt[1]) + "   sigma = " + str(popt[2]) )

actual_counts = 0
estimated_counts = 0

for estimated_count,count in zip(y_values,y2_values):
    estimated_counts += estimated_count
    actual_counts += count

print(actual_counts)
print(estimated_counts)