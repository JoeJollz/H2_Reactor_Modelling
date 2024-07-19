# -*- coding: utf-8 -*-
"""
Created on Sun Jul 14 17:48:31 2024

@author: jrjol

Quick plot of the Saturation Temperature (K) vs Saturation Pressure (bar)

Quick plot of the Saturation Temperature (deg C) vs Saturation Pressure (bar) 

"""

import numpy as np
import math as m
import matplotlib.pyplot as plt

temps_K = np.linspace(273.15+-99,273.15+234.95, 200) # kelvins
temps_degC = temps_K - 273.15 # deg C

# constants #
A = 4.42448
B = 1312.253
C = -32.445

log_Pressure_bar = A - (B/(temps_K+C))
Sat_Pressure_bar = 10**(log_Pressure_bar)

plt.plot(temps_K, Sat_Pressure_bar)

plt.show()

plt.plot(temps_degC, Sat_Pressure_bar)
plt.title('Acetone Phase Diagram')
plt.xlabel('Temperature °C')
plt.ylabel('Pressure (bar)')
plt.xlim(-150, 350)
plt.ylim(0, 100)
plt.axvline(x=234.95, color='red', linestyle='--')
plt.axvline(x = -94.65, color = 'green', linestyle = '--')
plt.axhline(y = 46.924, xmin = (-94.65 - -150) / (350 - -150), color = 'grey', linestyle = '--')

# Adding the text labels
plt.text(-121, 50, 'Solid\nPhase', fontsize=12, ha='center')
plt.text(50, 35, 'Liquid Phase', fontsize=12, ha='center')
plt.text(260, 10, 'Gaseous Phase', fontsize=12, ha='center')
plt.text(50, 75, 'Compressible Liquid', fontsize=12, ha='center')
plt.text(295, 80, 'Supercritical\nPhase', fontsize=12, ha='center')

# Adding the markers for Critical Point and Triple Point
plt.scatter(234.95, 46.924, color='black', label='Critical Point')
plt.text(234.95, 46.924, 'Critical Point', fontsize=12, ha = 'center', va='bottom')

plt.scatter(-99, 0, color='orange', label='Triple Point')
plt.text(-99, 0, 'Triple Point', fontsize=12, ha='left', va='bottom')


plt.show()