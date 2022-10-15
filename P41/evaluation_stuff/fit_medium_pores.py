import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

def func(x, a, b):
    return a/pow(x, b)

home_folder = os.getenv("HOME")
file_path = home_folder + "/build/P41/mediumPoreRadii.txt"
file_input = open(file_path, "r")

lines = file_input.readlines()
lines.pop(0)

x_values = []
y_values = []

for line in lines:
    x_low, x_up, y = ( float(ele) for ele in line.split() )
    x_values.append( (x_low+x_up)/2.0 )
    y_values.append( y )


plt.plot(x_values, y_values, 'b-')

popt, pcov = curve_fit(func, x_values, y_values)

plt.plot(x_values, func(x_values,*popt), 'r-')
plt.show()

print "fitted parameters:"
print("a = " + str(popt[0]) + "   b = " + str(popt[1]) )

def func2(x,a):
    return a/x

popt2, pcov2 = curve_fit(func2,x_values, y_values)
print("fitted parameter 2:" + str(popt2))

