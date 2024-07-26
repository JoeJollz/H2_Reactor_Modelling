'''
(Gas) Mass Flow Controllers (MFCs) are calibrated for a specific species, e.g. CO.
The MFCs measure the flow rate using the Constant Temperature Anemometry 
working principle. 2 probes are inserted in a the flow channel, the upstream 
probe is a heater, the downstream probe is a temperature sensor. 

The heater energy required to maintain the delta T between the 2 probes is
proportional to the mass flow rate and is thus a measure of the mass flow of 
the gas. 

The proportionality factor depends on thermal properties of the gas such as
its thermal conductivity and specific heat capcity, dynamic viscoty, nomarlised
density.

Whilst MFCs are calibrated for a particular gas, Conversion Factors (CFs) exist
so that is you want to use a CO MFC for CO2 gas, the CFs allow the MFC to allow
for the desired mass flow for the CO2 gas. The mass flow controller thinks
it has got CO gas! The conversion factors exist for each gas, and are referenced
to a common gas, such as Air. 

This file is for F-201CV-20K MFCs operating at 20 deg C
'''

# guage pressures
CO_Air_CF = dict()

CO_Air_CF[0] = 0.9984
CO_Air_CF[1] = 0.9983
CO_Air_CF[2] = 0.9981
CO_Air_CF[3] = 0.9980
CO_Air_CF[4] = 0.9979
CO_Air_CF[5] = 0.9978

CO2_Air_CF = dict()

CO2_Air_CF[0] = 0.7481
CO2_Air_CF[1] = 0.7444
CO2_Air_CF[2] = 0.7406
CO2_Air_CF[3] = 0.7369
CO2_Air_CF[4] = 0.7331
CO2_Air_CF[5] = 0.7293


