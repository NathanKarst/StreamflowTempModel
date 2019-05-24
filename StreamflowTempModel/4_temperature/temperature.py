import numpy as np
import time
import copy
import warnings
import scipy.optimize
import sys
import os
parent_dir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.append(os.path.join(parent_dir,'StreamflowTempModel','lib'))
import meteolib as meteo


class Temperature:
    """ Abstract base class for all temperature  models. (Test) 
    
    """
    def __init__(self, rew_id): self.rew_id = rew_id
                
    def update(self): 
        """ Update channel storage stock and volumetric discharge"""
        raise NameError('update')


class ChengHeatedGw(Temperature):
    """ 
 
    """
    def __init__(self, rew_id, **kwargs):
        Temperature.__init__(self, rew_id)
        
        args = ['Vf','gradient','eps','rho','cp','sigma','k_f', 'temperature', 'Tgw_hillslope']
        for arg in args: setattr(self, arg, kwargs[arg])        

        self.internalCounter = 0

    def update(self, dt, **kwargs):
        """ Update channel temperature
        
        Args:
            - kwargs (dict) : dictionary of inputs required to integrate time step
                - vol_1 : upstream volumetric discharge from stream 1
                - temp_1 : upstream temperature from stream 1x
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
        # day of year
        doy = kwargs['doy']

        # volumetric fluxes
        vol_1 = 1.15741e-11*kwargs['vol_1']
        vol_2 = 1.15741e-11*kwargs['vol_2']
        hillslope_volumetric_discharge = 1.15741e-11*kwargs['hillslope_volumetric_discharge']
        hillslope_volumetric_overlandFlow = 1.15741e-11*kwargs['hillslope_volumetric_overlandFlow']
        q = 1.15741e-11*kwargs['volumetric_discharge']

        # thermal depth from model 
        # dh = np.exp(-self.a_dh)*kwargs['upstream_area']**self.b_dh*q**self.c_dh

        # lengths
        length = 0.01*kwargs['length']
        width = 0.01*kwargs['width']

        # energy fluxes
        Sin = kwargs['Sin']

        # heat transfer coefficient from Cheng 2016, where va is wind speed
        # hc = 1.5*10**6*1.74*10**(-6)*(1 + 0.72*kwargs['va'])
        hc = 1.5*10**6*1.74*10**(-6)*(1 + 0)

        # temperatures 
        temp_1 = kwargs['temp_1'] + 273.15
        temp_2 = kwargs['temp_2'] + 273.15
        Ta = kwargs['ta'] + 273.15
        Ta_mean = kwargs['ta_mean'] + 273.15
        temp_curr = self.temperature + 273.15
        tau =length*self.Vf/hillslope_volumetric_discharge
        Tgw = self.Tgw_hillslope*np.exp(-self.k_f*tau) + (Ta_mean-273.15)*(1-np.exp(-self.k_f*tau)) + 273.15


        depth = 1e-6*kwargs['volume']/(length*width)
        u = q/(width*depth)
        # heat exchange coefficient, Cheng 2016
        k = (4*273**3*self.sigma*self.eps + hc)/(depth*self.rho*self.cp)

        # timestep        
        dt = dt*86400

        #compute distance that outlet parcel travels during this timestep
        distance = np.min([u*dt, length])

        #interpolate between upstream temp and current outlet temp
        #upstream temperature is assumed to be a mixture of groundwater
        #overland flow, and upstream flow links
        if vol_1+vol_2==0: 
            temp_up =  Tgw*hillslope_volumetric_discharge/q + (Ta_mean)*hillslope_volumetric_overlandFlow/q
        else:
            total_vol_in = vol_1 + vol_2 + hillslope_volumetric_discharge + hillslope_volumetric_overlandFlow
            temp_up = 1/total_vol_in*(vol_1*temp_1 + vol_2*temp_2 + hillslope_volumetric_discharge*Tgw + Ta_mean*hillslope_volumetric_overlandFlow)
            temp_up = temp_up
        
        temp_start = (temp_curr-temp_up)/length*(length-distance) + temp_up
        tnew = Ta + Sin/(k*depth*self.rho*self.cp) + (temp_start - Ta - Sin/(k*depth*self.rho*self.cp))*np.exp(-k*dt)
        self.temperature = tnew - 273.15
        return (depth, u)



class LagrangianSimpleTemperatureChengHeatedGW(Temperature):
    """ 
 
    """
    def __init__(self, rew_id, **kwargs):
        Temperature.__init__(self, rew_id)
        
        args = ['alphaw','eps','rho','cp','kh','sigma','temperature', 'kf', 'tau0', 'ktau', 'Tgw_offset']
        for arg in args: setattr(self, arg, kwargs[arg])        

        self.internalCounter = 0

    def update(self, dt, **kwargs):
        """ Update channel temperature
        
        Args:
            - kwargs (dict) : dictionary of inputs required to integrate time step
                - upstream_area: 
                - mannings_n: 
                - gradient: 
                - alphaw: 
                - eps: 
                - rho:
                - cp:
                - kh: 
                - sigma: 
                - temperature: 
                - kf: 
                - tau0: 
                - ktau: 
                - Tgw_offset: 
                - e


        
        Returns: 
           
        """
 
        # volumetric fluxes
        vol_1 = 1.15741e-11*kwargs['vol_1']
        vol_2 = 1.15741e-11*kwargs['vol_2']
        hillslope_volumetric_discharge = 1.15741e-11*kwargs['hillslope_volumetric_discharge']
        hillslope_volumetric_overlandFlow = 1.15741e-11*kwargs['hillslope_volumetric_overlandFlow']
        q = 1.15741e-11*kwargs['volumetric_discharge']

        # geometry
        length = 0.01*kwargs['length']
        # width here is the width of the upper water surface
        width = 0.01*kwargs['width']
        u = 1.15740741e-7*kwargs['u']
        volume = 1e-6*kwargs['volume']
        # this represents an effective depth (since depth varies across channel if not rectangular)
        depth = 0.01*kwargs['depth']


        # temperatures 
        temp_1 = kwargs['temp_1'] + 273.15
        temp_2 = kwargs['temp_2'] + 273.15
        Ta = kwargs['ta'] + 273.15
        Ta_mean = kwargs['ta_mean'] + 273.15
        temp_curr = self.temperature + 273.15


        # atmospheric data (water vapor pressure in kPa)
        ea = kwargs['ea']
        
        # energy fluxes
        Lin = kwargs['Lin']
        Sin = kwargs['Sin']

        # get groundwater temperature
        tau = self.tau0*np.exp(-self.ktau*hillslope_volumetric_discharge)
        Tgw = (self.Tgw_offset+273.15)*np.exp(-self.kf*tau) + Ta_mean*(1-np.exp(-self.kf*tau))

        # timestep to seconds       
        dt = dt*86400

        # compute distance that outlet parcel travels during this timestep, limited to channel length
        distance = np.min([u*dt, length])

        #interpolate between upstream temp and current outlet temp
        if vol_1+vol_2==0: 
            temp_up = Tgw 
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

        esat = 0.611*np.exp(2.5*10**6/461.0*(1/273.2 - 1/temp_start)) # saturation vapor pressure in kPa

        VTS = 0.9
        land_cover = 0.96*(1-VTS)*self.sigma*(Ta_mean)**4

        # groundwater added per unit length channel
        added_groundwater = hillslope_volumetric_discharge/length*dt
        added_overlandFlow = hillslope_volumetric_overlandFlow/length*dt

        # initial volume per unit length channel
        vol_initial = volume

        # # temperature in the channel changes as a function of: 
        # tnew -= dt*self.kh*(temp_start - Ta)/(self.rho*self.cp*depth)## sensible ChengHeatedGw
        tnew -= dt*100*(temp_start - Ta)/(self.rho*self.cp*depth)## sensible ChengHeatedGw

        # tnew += dt*(1-self.alphaw)*Sin/(self.rho*self.cp*depth) ## solar
        tnew += dt*Lin/(self.rho*self.cp*depth) ## atmospheric
        # tnew += dt*land_cover/(self.rho*self.cp*depth) ## atmospheric
        # tnew -= dt*back_radiation/(self.rho*self.cp*depth) ## back radiation 

        # tnew cannot be less than 273.15 or greater than 40 313.15
        tnew = np.min([np.max([273.15, tnew]), 313.15])

        ## Garner, What causes cooling water..., 2014
        # tnew += 2*dt*285.9*(0.132 + 0.143*kwargs['va'])*(ea - esat)/(self.rho*self.cp*depth) ## latent heat

        tnew = (tnew*vol_initial + (Tgw)*added_groundwater + Ta*added_overlandFlow)/(vol_initial + added_groundwater + added_overlandFlow)

        self.temperature = tnew - 273.15

        Tgw = Tgw - 273.
        Ta_mean = Ta_mean - 273.
        return (Tgw, Ta_mean)

class ImplicitEulerWesthoff(Temperature):
    """ 
 
    """
    def __init__(self, rew_id, **kwargs):
        Temperature.__init__(self, rew_id)
        
        args = ['alphaw','rho','cp','kh','sigma','temperature', 'kf', 'tau0', 'ktau', 'Tgw_offset']
        for arg in args: setattr(self, arg, kwargs[arg])        


        self.internalCounter = 0

    def update(self, dt, **kwargs):
        """ Update channel temperature
        
        Args:
            - kwargs (dict) : dictionary of inputs required to integrate time step
        
        Returns: 
           
        """
 
        # volumetric fluxes
        vol_1 = 1.15741e-11*kwargs['vol_1']
        vol_2 = 1.15741e-11*kwargs['vol_2']
        hillslope_volumetric_discharge = 1.15741e-11*kwargs['hillslope_volumetric_discharge']
        hillslope_volumetric_overlandFlow = 1.15741e-11*kwargs['hillslope_volumetric_overlandFlow']
        Q = 1.15741e-11*kwargs['volumetric_discharge']
        ppt = 1.157e-7*kwargs['ppt']

        # geometry
        length = 0.01*kwargs['length']
        # width here is the width of the upper water surface
        width = 0.01*kwargs['width']
        volume = 1e-6*kwargs['volume']
        area = volume/length
        u = Q/area
        # this represents an effective depth (since depth varies across channel if not rectangular)
        depth = volume/(length*width)

        # temperatures 
        temp_1 = kwargs['temp_1'] + 273.15
        temp_2 = kwargs['temp_2'] + 273.15
        Ta = kwargs['ta'] + 273.15
        Ta_mean = kwargs['ta_mean'] + 273.15
        temp_curr = self.temperature + 273.15

        # volume weighted input temps
        Qin = vol_1 + vol_2
        if Qin>0:
            Tin = vol_1*temp_1/Qin + vol_2*temp_2/Qin
        else: 
            Tin = 11.0

        # atmospheric data (water vapor pressure in kPa)
        ea = kwargs['ea']
        
        # timestep to seconds       
        dt = dt*86400
        
        # get groundwater temperature
        tau = self.tau0*np.exp(-self.ktau*hillslope_volumetric_discharge)
        Tgw = (self.Tgw_offset+273.15)*np.exp(-self.kf*tau) + Ta_mean*(1-np.exp(-self.kf*tau))

        # atmospheric data (water vapor pressure in kPa)
        ea = kwargs['ea']
        
        # energy fluxes
        Lin = kwargs['Lin'] # W/m-2
        Sin = kwargs['Sin'] # W/m-2
        esat = lambda temp: 0.611*np.exp(2.5*10**6/461.0*(1/273.2 - 1/temp)) # saturation vapor pressure in kPa
        
        VTS = 0.9
        land_cover = 0.96*(1-VTS)*0.96*self.sigma*(Ta)**4


        # now add in various heat fluxes
        back_radiation = lambda temp: 0.96*self.sigma*(temp)**4 # W/m-2
        sensible = lambda temp: -self.kh*(temp - Ta) # W/m-2
        shortwave = (1-self.alphaw)*Sin
#         latent = lambda temp: -self.latent_coefficient*(esat(temp) - ea)

        
        phi = lambda temp: Lin - back_radiation(temp) + sensible(temp) + shortwave  + land_cover #+ latent(temp)


        # groundwater added per unit length channel m^3/s/m
        qsub = hillslope_volumetric_discharge/length
        qoverland = hillslope_volumetric_overlandFlow/length
        qppt = ppt*width
        rhs = lambda temp: 1/area*(Qin/length*(Tin - temp) + qsub*(Tgw - temp) + qoverland*(Ta - temp) + width*phi(temp)/(self.rho*self.cp))
        solve = lambda tnew: temp_curr - tnew + dt*rhs(tnew)
        try: 
            tnew = scipy.optimize.newton(solve, temp_curr)
        except:
            tnew = temp_curr

        self.temperature = np.max([0,tnew - 273.15])

        Tgw = Tgw - 273.
        Ta_mean = Ta_mean - 273.
        return (Tgw, Ta_mean)





