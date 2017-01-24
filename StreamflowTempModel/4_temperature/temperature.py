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
        # volumetric fluxes
        vol_1 = 1.15741e-11*kwargs['vol_1']
        vol_2 = 1.15741e-11*kwargs['vol_2']
        hillslope_volumetric_discharge = 1.15741e-11*kwargs['hillslope_volumetric_discharge']
        hillslope_volumetric_overlandFlow = 1.15741e-11*kwargs['hillslope_volumetric_overlandFlow']
        volumetric_discharge = 1.15741e-11*kwargs['volumetric_discharge']

        # linear fluxes
        ppt = 1.15741e-7*kwargs['ppt']

        # temperatures 
        temp_1 = kwargs['temp_1'] + 273.15
        temp_2 = kwargs['temp_2'] + 273.15
        Ta = kwargs['ta'] + 273.15
        temp_curr = self.temperature + 273.15

        # lengths
        width = 0.01*kwargs['width']
        length = 0.01*kwargs['length']

        # volumes
        volume = 1e-6*kwargs['volume']

        # atmospheric data (water vapor pressure in kPa)
        ea = kwargs['ea']
        
        # energy fluxes
        Lin = kwargs['Lin']
        Sin = kwargs['Sin']
        Lout = self.eps*self.sigma*(temp_curr)**4
        esat = 0.611*np.exp(2.5*10**6/461.0*(1/273.2 - 1/temp_curr)) # saturation vapor pressure in kPa
        u = 1.0 # windspeed in m/s; assume 0.5 m/s, Allen 1998
        

        dt = dt*86400

        #need volumetric update back in here...
        flux = dt*ppt*length*width + vol_1*dt + vol_2*dt + hillslope_volumetric_discharge*dt  \
            + hillslope_volumetric_overlandFlow*dt - volumetric_discharge*dt
        vold = volume - flux
        # Using equations from Westhoff et al. 2007, HESS
        upstream_in = vol_1*temp_1 + vol_2*temp_2
        lateral_in = ppt*length*width*Ta + hillslope_volumetric_discharge*(self.Tgw+273.15) + hillslope_volumetric_overlandFlow*(Ta)
        downstream_out = volumetric_discharge*temp_curr
        depth = volume/(length*width)
        tnew = (vold*temp_curr + dt*(upstream_in + lateral_in - downstream_out ))/volume
        tnew -= 2.0*dt*self.kh*(temp_curr - Ta)/(self.rho*self.cp*depth)
        tnew += dt*(1-self.alphaw)*Sin/(self.rho*self.cp*depth)
        tnew += dt*Lin/(self.rho*self.cp*depth)
        tnew -= dt*Lout/(self.rho*self.cp*depth)
        tnew += dt*285.9*(0.132 + 0.143*u)*(ea - esat)/(self.rho*self.cp*depth)
        self.temperature = tnew - 273.15
        
        

        # temp_next -= dt*self.kh*(temp_curr - Ta)/(self.rho*self.cp*volume/(length*width))
        # temp_next += dt*(1-self.alphaw)*Sin/(self.rho*self.cp*volume/(length*width))
        # temp_next += dt*Lin/(self.rho*self.cp*volume/(length*width))
        # temp_next -= dt*Lout/(self.rho*self.cp*volume/(length*width))
        # new_energy -= dt*self.kh*(temp_curr - Ta)*length*width
        # new_energy += dt*(1-self.alphaw)*Sin*length*width
        # new_energy -= dt*Lin*length*width
        # new_energy -= dt*Lout*length*width


        # self.temperature += dt*285.9*(0.132 + 0.143*0.5)*(ea - esat)/(self.rho*self .cp*volume/(length*width))







