import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
import math

def normal_distribution(x, mu, sigma):
    return 1.0/np.sqrt(2*math.pi*sigma**2)*np.exp(-(x-mu)**2/(2*sigma**2))

def test_function(x, A, mu1, sigma1, B, mu2, sigma2):
    return A*normal_distribution(x,mu1,sigma1)+B*normal_distribution(x,mu2,sigma2)

def boundary_correction(x):
    return ((750-x)/750)**3

home_folder = os.getenv("HOME")
file_path = home_folder + "/build/P41/pores_larger_7um_downsampled_log.txt"
file_input = open(file_path, "r")

lines = file_input.readlines()
lines.pop(0)

x_values = []
y_values = []
y2_values = []

for line in lines:
    x_low, x_up, y = ( float(ele) for ele in line.split() )
    x_values.append( (math.log(2*x_low)+math.log(2*x_up))/2.0 )
    y_values.append( y/boundary_correction( np.exp(x_values[-1]) ) )
    y2_values.append(y)


plt.plot(x_values, y_values, 'b-')
plt.plot(x_values, y2_values, 'g-')

popt, pcov = curve_fit(test_function, x_values, y_values)
plt.plot(x_values, test_function(x_values,*popt), 'r-')
plt.show()

print("fitted parameters:")
print("A = " + str(popt[0]) + "   mu1 = " + str(popt[1]) + "   sigma1 = " + str(popt[2]) )
print("B = " + str(popt[3]) + "   mu2 = " + str(popt[4]) + "   sigma2 = " + str(popt[5]) )

# estimate total number of pores
actual_pores = 0
estimated_pores = 0
for y_value,y2_value in zip(y_values,y2_values) : 
    actual_pores += y2_value
    estimated_pores += y_value

print(actual_pores)
print(estimated_pores)
