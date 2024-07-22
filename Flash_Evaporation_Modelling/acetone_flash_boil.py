"#This is acetone flash boiling model - single equilibrium adiabatic flash expansion"

'''
required for correct modelling. 

error checker, 
1. is the lower P less than saturationP at that given temp.

'''

import math as m
import numpy as np
from CoolProp.CoolProp import PropsSI
import Acetone_Equation_of_State

class OperatingError(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class Flash_Model:
    def __init__(self, reaction_name):
        self.reaction_name = reaction_name
        
    
    def operating_error_check(self, temperature, pressure_upstream, pressure_downstream):
        # Check upstream P and downstream P are going to initate a flash boil
        # of Acetone.
        # temperature °C input!
        
        saturation_pressure_bar, sat_p_Pa = Acetone_Equation_of_State.sat_T_to_sat_P(temperature)
        
        if pressure_upstream < sat_p_Pa:
            raise OperatingError("Operating at {temperature} °C, with an upstream pressure of {pressure_upstream}Pa, the system is already vapor. (Increase the upstream pressure to be > {sat_p_Pa}Pa")
        
        if pressure_downstream > sat_p_Pa:
            raise OperatingError("NO FLASH BOIL POSSIBLE. Downstream pressure is greater than the saturation pressure at this temperature. (For flash boil, downstream pressure must be less than {sat_p_Pa}Pa")
        
        
    def H_TP(T,P,phase):
        '''
    
        Parameters
        ----------
        T : Float
            Temperature in °C.
        P : Float
            Pressure in Pa.
        phase : Int
            Phase.
            0 for gas.
            1 for liquid.
            2 for two-phase region.
            3 for Supercritical liquid.
            4 for Supercritical gas.
            5 for Supercritical. 
        Returns
        -------
        Enthalpy : Float
            Enthalpy of the working fluid.
    
        ''' 
    #TO DO: possibly add a checker in to check it is liquid or gas. 
        if phase == 0:
            # gas
            Enthalpy = PropsSI('H', 'T|gas', T, 'P', P, 'Acetone') # 
            
        elif phase == 1:
            # liquid
            Enthalpy = PropsSI('H', 'T|liquid', T, 'P', P, 'Acetone') # 
            
        elif phase == 2:
            # two phase
            Enthalpy = PropsSI('H', 'T|twophase', T, 'P', P, 'Acetone') # 
            
        elif phase == 3:
            #supercritical liquid. p>pcrit & T<Tcrit
            Enthalpy = PropsSI('H', 'T|supercritical_liquid', T, 'P', P, 'Acetone') # 
            
        elif phase ==4:
            # supercritical gas, p<pcrit & T>T crit.
            Enthalpy = PropsSI('H', 'T|supercritical_gas', T, 'P', P, 'Acetone') # 
          
        elif phase ==5:
            # supercritical gas, p>pcrit & T>T crit.
            Enthalpy = PropsSI('H', 'T|supercritical', T, 'P', P, 'Acetone') #    
        
        return Enthalpy  
        
    def vapor_content(self, temperature, pressure_upstream, pressure_downstream):
        '''
        Parameters
        ----------
        temperature : FLOAT
            Operating temperature (adiabatic) (°C).
        pressure_upstream : FLOAT
            Upstream pressure (Pa).
        pressure_downstream : FLOAT
            Downstream pressure (Pa)

        Returns
        -------
        vapor_content : FLOAT
            Mass ratio, vapor mass/total mass.

        '''
        self.operating_error_check
            
        Enthalpy_upstream_liquid = self.H_TP(temperature, pressure_upstream, int(1))
        Enthalpy_downstream_liquid = self.H_TP(temperature, pressure_downstream, int(1))
        Enthalpy_downstream_vapor = self.H_TP(temperature, pressure_downstream, int(0))
            
        vapor_content = (Enthalpy_upstream_liquid-Enthalpy_downstream_liquid)/(Enthalpy_downstream_vapor-Enthalpy_downstream_liquid)    
        
        return vapor_content
    
if __name__ == "__main__":
    # Create an instance of the Reaction class
    reaction = Flash_Model("Example Reaction")
    
    # Define test inputs
    temperature = 100  # degrees Celsius
    pressure_upstream = 2000000  # Pascals
    pressure_downstream = 100000  # Pascals
    
    # Perform the calculation
    vapor_content = Flash_Model.vapor_content(temperature, pressure_upstream, pressure_downstream)
    
    # Print the result
    print(f"Vapor Content: {vapor_content}")