import numpy as np
from matplotlib import pyplot as plt
import time
class Channel:
    """ Abstract base class for all channel  models. (Test) 
    
    """
    def __init__(self, rew_id): self.rew_id = rew_id
                
    def update(self): 
        """ Update channel storage stock and volumetric discharge"""
        raise NameError('update')     




class SimpleChannel(Channel):
    """ Simple hydrologic channel model.
    
    Args:
        - rew_id (int) : id of rew to which the channel belongs
        - volume (float): [L^3]
        - volumetric_discharge (float): [L^3/T]
        - mannings_n (int): 
        - width (float): [L] channel width
        - length (float): [L] channel length
        - slope (float): [] channel slope
        
    """
    def __init__(self, rew_id, rew_config, volume = 0, mannings_n = 0.03, e=0.01, f=0.39):
        Channel.__init__(self, rew_id)
        self.volume = volume
        self.mannings_n = mannings_n
        self.width = e*rew_config.upstream_area.loc[rew_id]**f
        self.length = rew_config.length.loc[rew_id]
        self.slope = rew_config.gradient.loc[rew_id]
        h = self.volume/(self.length*self.width)
        self.volumetric_discharge = _manning_u(h, self.mannings_n, self.slope)*self.width*h

    def update(self, dt, upstream_volumetric_discharge, hillslope_volumetric_discharge, ppt):
        """ Update channel storage stock and volumetric discharge
        
        Args:
            - dt (float): [T] time step
            - upstream_volumetric_discharge (float): [L^3/T] 
            - hillslope_volumetric_discharge (float): [L^3/T]
            - ppt (float): [L/T] precipitation flux falling channel itself
        
        Returns: 
            - Boolean variable indicating whether or not instability approximation was used to solve kinematic wave.
        
        """

        h = self.volume/(self.length*self.width)
        if _manning_u(h, self.mannings_n, self.slope)*dt>self.length:
            self.volumetric_discharge = h*self.width*self.length/dt
            self.volume = 0
            self.volume += ppt*dt*self.width*self.length
            self.volume += upstream_volumetric_discharge*dt
            self.volume += hillslope_volumetric_discharge*dt
            return 1
        else:
            self.volumetric_discharge = _manning_u(h, self.mannings_n, self.slope)*self.width*h
            self.volume += ppt*dt*self.width*self.length
            self.volume += upstream_volumetric_discharge*dt
            self.volume += hillslope_volumetric_discharge*dt
            self.volume -= _manning_u(h, self.mannings_n, self.slope)*self.width*h*dt
            return 0

def _manning_u(h, n, slope):
    #takes h in cm, uses SI mannings n values, returns u in cm/day
	return 100*86400*(h/100.)**(2/3.)*slope**(0.5)*1/n

