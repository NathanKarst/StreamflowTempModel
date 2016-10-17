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


    def update(self, 
        dt, 
        upstream_volumetric_discharge_1, 
        upstream_temperature_1,
        upstream_volumetric_discharge_2,
        upstream_temperature_2,
        hillslope_volumetric_discharge, 
        volumetric_discharge, 
        width,
        length,
        volume,
        Ta,
        Lin, 
        Sin, 
        ppt, 
        LPI
        ):
        """ Update channel temperature
        
        Args:
            - 
        
        Returns: 
           
        """
        dt = dt*86400
        Ta = Ta + 273.15
        upstream_temperature_1 = upstream_temperature_1 + 273.15
        upstream_temperature_2 = upstream_temperature_2 + 273.15
        temp_curr = self.temperature + 273.15
        Lout = self.sigma*(temp_curr)**4


        self.temperature += dt*(ppt*length*width*Ta + upstream_volumetric_discharge_1*upstream_temperature_1 + upstream_volumetric_discharge_2*upstream_temperature_2 + hillslope_volumetric_discharge*(self.Tgw + 273.15) - volumetric_discharge*temp_curr)/volume 
        self.temperature -= dt*self.kh*(temp_curr - Ta)/(self.rho*self.cp*volume/(length*width))
        self.temperature += LPI*dt*(1-self.alphaw)*Sin/(self.rho*self.cp*volume/(length*width))
        self.temperature += LPI*dt*Lin/(self.rho*self.cp*volume/(length*width))
        self.temperature -= LPI*dt*Lout/(self.rho*self.cp*volume/(length*width))






