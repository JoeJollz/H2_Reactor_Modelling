"MFC Class conversion main" 

import MFC_Converter


class MassFlowController:
    def __init__(self, mfc_number):
        self.mfc_number = mfc_number
        
    def CO2_CO_reading(self, desired_mass_flow_lpm, guage_pressure):
        '''
        

        Parameters
        ----------
        desired_mass_flow_lpm : float
            Litres per min of CO2 gas.
        guage_pressure : int
            Gauage pressure (bar).

        Returns
        -------
        vol_flow : flaot
            equiv volume flow of CO2 litres per min.
        rate : float
            0-32000 is equivalent to 0-100% MFC valve open.

        '''
        # assuming 32000 corresponds to 100% flow. 
        guage_pressure = int(guage_pressure)
        
        CF_co_air = MFC_Converter.CO_Air_CF[2]
        CF_co2_air = MFC_Converter.CO2_Air_CF[2]
        
        vol_flow = desired_mass_flow_lpm/CF_co2_air*CF_co_air
        
        rate = vol_flow/20*32000
        
        return vol_flow, rate
        
if __name__ == "__main__":
    
    mfc = MassFlowController(1)
    
    desired_flow_rate_CO2 = 1 #lpm
    guage_pressure = 1 #bar
    
    mass_flow, baudrate = mfc.CO2_CO_reading(desired_flow_rate_CO2,guage_pressure)
    
    print(f"Mass flow to set for MFC number 1: {mass_flow} ln/min, with rate {baudrate}")
    
    
    
        