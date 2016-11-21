import numpy as np
from matplotlib import pyplot as plt
import time
class Temperature:
    """ Abstract base class for all temperature  models. (Test) 
    
    """
    def __init__(self, rew_id): self.rew_id = rew_id
                
    def update(self): 
        """ Update channel storage stock and volumetric discharge"""
        raise NameError('update')


class SimpleTemperature(Temperature):
    """ Simple temp model 
    
    dT/dt=(QinTin+QgwTgw-QoutTout)/V+(1-alphaw)phi/(d*rho*cp) - Lout + Lin - Kh*(T - Ta) 

    Args:
        - rew_id (int) : id of rew to which the channel belongs
        - alphaw (float): albedo of water
        - eps (float): emissivity of water
        - rho (float): density of water
        - Tgw (float): incoming groundwater temperature
        - cp (float): specific heat of water
        - kh (float): Turbulent heat exchange coefficient, Van Beek et al (2012)
        - sigma (float): 
    """
    def __init__(self, rew_id, **kwargs):
        Temperature.__init__(self, rew_id)
        
        args = ['alphaw','eps','rho','cp','kh','sigma','Tgw','temperature']
        for arg in args: setattr(self, arg, kwargs[arg])        


    def update(self, dt, **kwargs):
        """ Update channel temperature
        
        Args:
            - kwargs (dict) : dictionary of inputs required to integrate time step
                - vol_1 : upstream volumetric discharge from stream 1
                - temp_1 : upstream temperature from stream 1
                - vol_2 : upstream volumetric discharge from stream 2
                - temp_2 : upstream temperature from stream 2
                - hillslope_volumetric_discharge : discharge to stream from REW groundwaterZone
                - hillslope_volumetric_overlandFlow : overland flow from REW to stream
                - volumetric_discharge : discharge leaving stream link in current timestep
                - width : channel width
                - length : channel length
                - Ta : air temperature 
                - Lin : longwave incoming radiation
                - Sin : shortwave incoming radiation
                - ppt : incoming precipitation

        
        Returns: 
           
        """
        vol_1 = kwargs['vol_1']
        temp_1 = kwargs['temp_1']
        vol_2 = kwargs['vol_2']
        temp_2 = kwargs['temp_2']
        hillslope_volumetric_discharge = kwargs['hillslope_volumetric_discharge']
        hillslope_volumetric_overlandFlow = kwargs['hillslope_volumetric_overlandFlow']
        volumetric_discharge = kwargs['volumetric_discharge']
        width = kwargs['width']
        length = kwargs['length']
        volume = kwargs['volume']
        Ta = kwargs['Ta']
        Lin = kwargs['Lin']
        Sin = kwargs['Sin']
        ppt = kwargs['ppt']
        
        dt = dt*86400
        Ta = Ta + 273.15
        temp_1 = temp_1 + 273.15
        temp_2 = temp_2 + 273.15
        temp_curr = self.temperature + 273.15
        Lout = self.sigma*(temp_curr)**4

 
        self.temperature += dt*(ppt*length*width*Ta + vol_1*temp_1 + vol_2*temp_2 + hillslope_volumetric_discharge*(self.Tgw + 273.15) + hillslope_volumetric_overlandFlow*(Ta) - volumetric_discharge*temp_curr)/volume 
        self.temperature -= dt*self.kh*(temp_curr - Ta)/(self.rho*self.cp*volume/(length*width))
        self.temperature += dt*(1-self.alphaw)*Sin/(self.rho*self.cp*volume/(length*width))
        self.temperature += dt*Lin/(self.rho*self.cp*volume/(length*width))
        self.temperature -= dt*Lout/(self.rho*self.cp*volume/(length*width))






