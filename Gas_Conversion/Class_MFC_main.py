"MFC Class conversion main" 

import MFC_Converter


class MassFlowController:
    def __init__(self, mfc_number):
        self.mfc_number = mfc_number
        
    def CO2_CO_reading(self, desired_mass_flow_lpm, guage_pressure):
        
        
        return mass_flow, baudrate
        
if __name__ == "__main__":
    
    mfc = MassFlowController(1)
    
    desired_flow_rate_CO2 = 1 #lpm
    guage_pressure = 1 #bar
    
    mass_flow, baudrate = mfc.CO2_CO_reading(desired_flow_rate_CO2,guage_pressure)
    
    print(f"Mass flow to set for MFC number 1: {mass_flow} ln/min, with baudrate {baudrate}")
    
    
    
        