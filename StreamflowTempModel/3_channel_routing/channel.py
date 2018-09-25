import numpy as np
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


class NoChannel(Channel):
    """ Channel discharge just equals upstream + hillslope, ignore PPT
    
    Args:
        - rew_id (int) : id of rew to which the channel belongs
        - volume (float): [L^3]
        - volumetric_discharge (float): [L^3/T]
        - width (float): [L] channel width
        - length (float): [L] channel length
        - slope (float): [] channel slope
        
    """
    def __init__(self, rew_id, **kwargs):
        Channel.__init__(self, rew_id)
        
        args = ['volume','gradient','length']
        for arg in args: setattr(self, arg, kwargs[arg])        

        self.volumetric_discharge = 0

        # for NoChannel, set width to NaN
        self.width = np.nan

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

        self.volumetric_discharge = upstream_volumetric_discharge + hillslope_volumetric_discharge 
        return 0

class TrapezoidalChannel(Channel):
    """ Trapezoidal hydrologic channel model using Manning's equation for open channel flow (for a reference, see doi:10.1103/PhysRevLett.88.014501.)
    	The channel has a bottom width and a bank slope, making a symmetric channel with sloped sides. If
    	the bottom width is zero, the channel is triangular. 
    Args:
        - rew_id (int) : id of rew to which the channel belongs
        - volume (float): [L^3]
        - volumetric_discharge (float): [L^3/T]
        - mannings_n (int): 
        - width (float): [L] channel width
        - length (float): [L] channel length
        - slope (float): [] channel slope
        - bank_slope (float): [] slope of left and right banks
        
    """
    def __init__(self, rew_id, **kwargs):
        Channel.__init__(self, rew_id)
        
        args = ['volume','mannings_n','upstream_area','gradient','length','e','f','g','h']
        for arg in args: setattr(self, arg, kwargs[arg])        

        # e and f are the power law scaling parameters for the bottom width
        # g and h are the power law scaling parameters for the bank slope
        self.width = self.e*self.upstream_area**self.f
        self.bank_slope = self.g*self.upstream_area**self.h
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
        xs_area = self.volume/self.length
        depth = (-(self.bank_slope*self.width) + np.sqrt(self.bank_slope)*np.sqrt(4*xs_area + self.bank_slope*self.width**2))/2.
        wetted_perimeter = self.width + 2*depth*np.sqrt(1+1/self.bank_slope**2)
        top_width = self.width + 2*depth/self.bank_slope
        h = xs_area/wetted_perimeter

        if _manning_u(h, self.mannings_n, self.gradient)*dt>self.length:
            self.volumetric_discharge = xs_area*self.length/dt
            self.volume = 0
            self.volume += ppt*dt*top_width*self.length
            self.volume += upstream_volumetric_discharge*dt
            self.volume += hillslope_volumetric_discharge*dt
            return 1
        else:
            self.volumetric_discharge = _manning_u(h, self.mannings_n, self.gradient)*xs_area
            self.volume += ppt*dt*top_width*self.length
            self.volume += upstream_volumetric_discharge*dt
            self.volume += hillslope_volumetric_discharge*dt
            self.volume -= self.volumetric_discharge*dt
            return 0


def _manning_u(h, n, slope):
    #takes h in cm (hydraulic radius), uses SI mannings n values, returns u in cm/day
	return 100*86400*(h/100.)**(2/3.)*slope**(0.5)*1/n

