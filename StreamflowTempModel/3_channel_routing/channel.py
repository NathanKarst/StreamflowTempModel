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
    """ Simple hydrologic channel model using Manning's equation for open channel flow (for a reference, see doi:10.1103/PhysRevLett.88.014501.)
    
    Args:
        - rew_id (int) : id of rew to which the channel belongs
        - volume (float): [L^3]
        - volumetric_discharge (float): [L^3/T]
        - mannings_n (int): 
        - width (float): [L] channel width
        - length (float): [L] channel length
        - slope (float): [] channel slope
        
    """
    def __init__(self, rew_id, **kwargs):
        Channel.__init__(self, rew_id)
        
        args = ['volume','mannings_n','upstream_area','gradient','length','e','f']
        for arg in args: setattr(self, arg, kwargs[arg])        

        self.width = self.e*self.upstream_area**self.f
        self.volumetric_discharge = 0

    def update(self, dt, upstream_volumetric_discharge, hillslope_volumetric_discharge, ppt):
        """ Update channel storage stock and volumetric discharge
        
        Args:
            - dt (float): [T] time step
            - upstream_volumetric_discharge (float): [L^3/T] 
            - hillslope_volumetric_discharge (float): [L^3/T]
            - ppt (float): [L/T] precipitation falling channel itself
        
        Returns: 
            - Boolean variable indicating whether or not instability approximation was used to solve kinematic wave.
        
        """

        h = self.volume/(self.length*self.width)
        if _manning_u(h, self.mannings_n, self.gradient)*dt>self.length:
            self.volumetric_discharge = h*self.width*self.length/dt
            self.volume = 0
            self.volume += ppt*dt*self.width*self.length
            self.volume += upstream_volumetric_discharge*dt
            self.volume += hillslope_volumetric_discharge*dt
            return 1
        else:
            self.volumetric_discharge = _manning_u(h, self.mannings_n, self.gradient)*self.width*h
            self.volume += ppt*dt*self.width*self.length
            self.volume += upstream_volumetric_discharge*dt
            self.volume += hillslope_volumetric_discharge*dt
            self.volume -= self.volumetric_discharge*dt
            return 0


def _manning_u(h, n, slope):
    #takes h in cm, uses SI mannings n values, returns u in cm/day
	return 100*86400*(h/100.)**(2/3.)*slope**(0.5)*1/n

