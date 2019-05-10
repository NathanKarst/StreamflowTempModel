import numpy as np
import time
import scipy
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
        self.base_width = self.e*self.upstream_area**self.f
        self.base_width = 0
        self.width = self.e*self.upstream_area**self.f
        self.bank_slope = self.g*self.upstream_area**self.h
        self.bank_slope = 45.0
        self.volumetric_discharge = 0
        self.depth = 0
        self.u = 0
        self.area = 0
        self.x = 0

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
        self.area = xs_area
        depth = (-(self.bank_slope*self.base_width) + np.sqrt(self.bank_slope)*np.sqrt(4*xs_area + self.bank_slope*self.base_width**2))/2.
        x = (-self.base_width*self.bank_slope + np.sqrt(self.bank_slope)*np.sqrt(4*xs_area + self.bank_slope*self.base_width**2))/(2*self.bank_slope)
        self.depth = depth
        side = np.sqrt(depth**2 + x**2)
        wetted_perimeter = side*2+self.base_width
        h = xs_area/wetted_perimeter
        top_width = self.base_width + 2*x
        self.width = top_width

        if _manning_u(h, self.mannings_n, self.gradient)*dt>self.length:
            self.u = self.length/dt
            self.volumetric_discharge = xs_area*self.length/dt
            self.volume = 0
            self.volume += ppt*dt*top_width*self.length
            self.volume += upstream_volumetric_discharge*dt
            self.volume += hillslope_volumetric_discharge*dt
            return 1
        else:
            self.u = _manning_u(h, self.mannings_n, self.gradient)
            self.volumetric_discharge = self.u*xs_area
            self.volume += ppt*dt*top_width*self.length
            self.volume += upstream_volumetric_discharge*dt
            self.volume += hillslope_volumetric_discharge*dt
            self.volume -= self.volumetric_discharge*dt
            return 0


class AllenChannel(Channel):
    """ Channel geometry from Alen et al 2018
    Args:
        
    """
    def __init__(self, rew_id, **kwargs):
        Channel.__init__(self, rew_id)
        
        args = ['volume','mannings_n','upstream_area','gradient','length','c', 'd', 'e','f','g','h']
        for arg in args: setattr(self, arg, kwargs[arg])        

        # c d are power law params for bankful depth
        # e f for bankful width
        # g h for r
        # area must be in km to match with units on coefficient
        self.bankful_depth = self.c*(1e-10*self.upstream_area)**self.d
        self.bankful_width = self.e*(1e-10*self.upstream_area)**self.f
        self.r = self.g*(1e-10*self.upstream_area)**self.h

        # in m^2
        xs_area = 0.0001*self.volume/self.length
        h = (xs_area/((1-1/(self.r+1))*self.bankful_width*self.bankful_depth**(-1/self.r)))**(1/(1+1/self.r))

        # def func(h):
        #     return self.area - _allen_a(h, self.bankful_depth, self.bankful_width, self.r)
        # # get a depth in m
        # self.depth = scipy.optimize.newton(func, 0.1)
        # units with m
        k = (8.1*9.81**0.5*self.mannings_n)**6
        u = 8.1*(9.81*h*self.gradient)**0.5*(h/k)**(1/6.0)
        volumetric_discharge = xs_area*u
        width = volumetric_discharge**(3/(5*self.r+3))*self.bankful_width**( (self.r-1)/(self.r + 3/5) )*(
            8.1*(9.81*self.gradient)**0.5*k**(-1/6.0)*(self.bankful_width/self.bankful_depth)**(-3/5.)*(1-1/(self.r+1)))**(-3/(5*self.r+3))


        self.volumetric_discharge = volumetric_discharge*86400000000
        self.width = width*0.01
        self.u = u*8640000
        self.area = xs_area*10000
        self.depth = h*0.01

        
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
        # volume and length are in cm, 
        xs_area_cm2 = self.volume/self.length
        xs_area = xs_area_cm2*0.0001
        dt_seconds = dt*24*60*60
        length_meters = self.length*0.01

        # all in meters
        h = (xs_area/((1-1/(self.r+1))*self.bankful_width*self.bankful_depth**(-1/self.r)))**(1/(1+1/self.r))
        k = (8.1*9.81**0.5*self.mannings_n)**6
        u = 8.1*(9.81*h*self.gradient)**0.5*(h/k)**(1/6.0)
        # in meters
        volumetric_discharge = xs_area*u
        width = volumetric_discharge**(3/(5*self.r+3))*self.bankful_width**( (self.r-1)/(self.r + 3/5) )*(
            8.1*(9.81*self.gradient)**0.5*k**(-1/6.0)*(self.bankful_width/self.bankful_depth)**(-3/5.)*(1-1/(self.r+1)))**(-3/(5*self.r+3))
        width_cm = width*100.0

        # u is in units of m/s
        if u*dt_seconds>length_meters:
            # store all with units of cm and days
            self.u = self.length/dt
            self.volumetric_discharge = xs_area_cm2*self.length/dt
            self.volume = 0
            self.volume += ppt*dt*width_cm*self.length
            self.volume += upstream_volumetric_discharge*dt
            self.volume += hillslope_volumetric_discharge*dt
            self.width = width_cm
            self.depth = h*100
            return 1
        else:
            # store u in cm/day
            self.u = 8640000*u
            self.volumetric_discharge = self.u*xs_area_cm2
            self.volume += ppt*dt*width_cm*self.length
            self.volume += upstream_volumetric_discharge*dt
            self.volume += hillslope_volumetric_discharge*dt
            self.volume -= self.volumetric_discharge*dt
            self.width = width_cm
            self.area = xs_area_cm2
            self.depth = h*100
            return 0


def _manning_u(h, n, slope):
    #takes h in cm (hydraulic radius), uses SI mannings n values, returns u in cm/day
	return 100*86400*(h/100.)**(2/3.)*slope**(0.5)*1/n

