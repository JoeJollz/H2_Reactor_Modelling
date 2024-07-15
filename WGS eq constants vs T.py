# -*- coding: utf-8 -*-
"""
Created on Sun Jul  7 10:34:12 2024

@author: jrjol

Water Gas Shift (WGS) reaction. Simple code of the WGS reaction constants vs Temperature. Noticeably, the reaction constants always favour the formation of CO2(g)+H2(g) in contrast to H2O(g)+CO(g). 
This is crucial to prove that for the temperature range my reactor is operating at, if I can oxidise my oxygen carries in a CO2 environment, I can certainly oxidise my oxygen carriers in an H2O environment, leaving H2(g)!

"""
import numpy as np
import matplotlib.pyplot as plt

T = np.linspace(298.15, 1073.15)

#Kp = np.exp(4577.8/T-4.33)
Kp_log = 5693.5/T + 1.077*np.log(T)+5.44*10**(-4)*T-1.125*10**(-7)*T**2-49170/(T**2)-13.148

plt.title('WGS ln(Kp) vs T (K)')
plt.plot(T,Kp_log,label = 'ln(Kp)')
plt.xlabel('Temperature (Kelvin)')
plt.ylabel('ln(Kp)')



 