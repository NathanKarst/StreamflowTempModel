import numpy as np
import time
class ChannelTemperature:
    """ Abstract base class for all combined channel/temperature  models. 
    
    """
    def __init__(self, rew_id): self.rew_id = rew_id
                
    def update(self): 
        """ Update channel storage stock and volumetric discharge"""
        raise NameError('update')

class EulerianSimpleChannelTemperature(ChannelTemperature):
    """ 
 		
    """
    def __init__(self, rew_id, **kwargs):
        ChannelTemperature.__init__(self, rew_id)
        
        args = ['upstream_area', 'mannings_n', 
        		'gradient', 'alphaw','eps','rho','cp','kh',
        		'sigma','temperature', 'kf', 'tau0', 'ktau', 
        		'Tgw_offset','volume',
        		'upstream_area','gradient','length','e','f','g']



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
        width = self.volume/
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
            return 0        doy = kwargs['doy']

        # channel angle
        angle = np.pi/180*self.angle

        # volumetric fluxes
        vol_1 = 1.15741e-11*kwargs['vol_1']
        vol_2 = 1.15741e-11*kwargs['vol_2']
        hillslope_volumetric_discharge = 1.15741e-11*kwargs['hillslope_volumetric_discharge']
        hillslope_volumetric_overlandFlow = 1.15741e-11*kwargs['hillslope_volumetric_overlandFlow']
        q = 1.15741e-11*kwargs['volumetric_discharge']

        # temperatures 
        temp_1 = kwargs['temp_1'] + 273.15
        temp_2 = kwargs['temp_2'] + 273.15
        Ta = kwargs['ta'] + 273.15
        Ta_mean = kwargs['ta_mean'] + 273.15
        temp_curr = self.temperature + 273.15


        # lengths
        length = 0.01*kwargs['length']

        width = 2**(5/8.)*self.mannings_n**(3/8.)*q**(3/8.)*np.tan(angle)**(3/8.)/(np.cos(angle)**.25*self.gradient**(3/16.))
        depth = self.mannings_n**(3/8.)*q**(3/8.)/(2**(3/8.)*np.cos(angle)**.25*self.gradient**(3/16.)*np.tan(angle)**(5/8.))
        u = np.sqrt(np.cos(angle))*q**.25*self.gradient**(3/8.)*np.tan(angle)**.25/(2**.25*self.mannings_n**(3/4.))

        # atmospheric data (water vapor pressure in kPa)
        ea = kwargs['ea']
        
        # energy fluxes
        Lin = kwargs['Lin']
        Sin = kwargs['Sin']

        # get groundwater temperature
        tau = self.tau0*np.exp(-self.ktau*hillslope_volumetric_discharge)
        Tgw = (self.Tgw_offset+273.15)*np.exp(-self.kf*tau) + Ta_mean*(1-np.exp(-self.kf*tau))


        # timestep        
        dt = dt*86400


        #compute distance that outlet parcel travels during this timestep
        distance = np.min([u*dt, length])

        #interpolate between upstream temp and current outlet temp
        if vol_1+vol_2==0: 
            temp_up =  Tgw 
        else:
            temp_up = vol_1/(vol_1+vol_2)*temp_1 + vol_2/(vol_1+vol_2)*temp_2
        
        temp_start = (temp_curr-temp_up)/length*(length-distance) + temp_up
        tnew = temp_start

        # now add in various heat fluxes
        back_radiation = 0.96*self.sigma*(temp_start)**4

        # atmospheric data (water vapor pressure in kPa)
        ea = kwargs['ea']
        
        # energy fluxes
        Lin = kwargs['Lin']
        Sin = kwargs['Sin']

        ## NJK: why self.eps here? compare to (18) in Westhoff 07. 
        Lout = 0.96*self.sigma*temp_start**4 ## suggested change
        esat = 0.611*np.exp(2.5*10**6/461.0*(1/273.2 - 1/temp_start)) # saturation vapor pressure in kPa
   
        dt = dt*86400

        VTS = 0.9
        land_cover = 0.96*(1-VTS)*self.sigma*(Ta_mean)**4

        added_groundwater = hillslope_volumetric_discharge/length*dt
        added_overlandFlow = hillslope_volumetric_overlandFlow/length*dt
        vol_initial = width*depth

        tnew -= 2*dt*self.kh*(temp_start - Ta)/(self.rho*self.cp*depth) ## sensible heat
        tnew += 2*dt*(1-self.alphaw)*Sin/(self.rho*self.cp*depth) ## solar
        tnew += 2*dt*Lin/(self.rho*self.cp*depth) ## atmospheric
        tnew += 2*dt*land_cover/(self.rho*self.cp*depth) ## atmospheric
        tnew -= 2*dt*back_radiation/(self.rho*self.cp*depth) ## back radiation 

        ## Garner, What causes cooling water..., 2014
        # tnew += 2*dt*285.9*(0.132 + 0.143*kwargs['va'])*(ea - esat)/(self.rho*self.cp*depth) ## latent heat

        tnew = (tnew*vol_initial + (Tgw)*added_groundwater + Ta*added_overlandFlow)/(vol_initial + added_groundwater + added_overlandFlow)

        self.temperature = tnew - 273.15

        Tgw = Tgw - 273.
        Ta_mean = Ta_mean - 273.0
        return (Tgw, Ta_mean, depth)


def _manning_u(h, n, slope):
    #takes h in cm, uses SI mannings n values, returns u in cm/day
	return 100*86400*(h/100.)**(2/3.)*slope**(0.5)*1/n