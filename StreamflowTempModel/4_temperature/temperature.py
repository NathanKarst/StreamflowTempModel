import numpy as np
import time
import copy
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


        # triangular channel
        # angle = 2*np.pi/60
        # width = np.sqrt(2*volume/np.tan(angle)/length)

        # atmospheric data (water vapor pressure in kPa)
        ea = kwargs['ea']
        
        # energy fluxes
        Lin = kwargs['Lin']
        Sin = kwargs['Sin']

        ## NJK: why self.eps here? compare to (18) in Westhoff 07. 
        Lout = self.eps*self.sigma*(temp_curr)**4 
        #Lout = 0.96*self.sigma*temp_curr**4 ## suggested change


        esat = 0.611*np.exp(2.5*10**6/461.0*(1/273.2 - 1/temp_curr)) # saturation vapor pressure in kPa
        u = 2.0 # windspeed in m/s; assume 0.5 m/s, Allen 1998
        

        dt = dt*86400

        ## NJK: why do we need this?
        #need volumetric update back in here...
        flux = dt*ppt*length*width + vol_1*dt + vol_2*dt + hillslope_volumetric_discharge*dt  \
            + hillslope_volumetric_overlandFlow*dt - volumetric_discharge*dt
        volume_next = volume + flux
        leftover_stream_volume = volume - volumetric_discharge*dt


        tnew = (temp_curr*leftover_stream_volume 
            + Ta*dt*ppt*length*width
            + temp_2*vol_2*dt + temp_1*vol_1*dt
            + (self.Tgw + 273.15)*hillslope_volumetric_discharge*dt
            + Ta*hillslope_volumetric_overlandFlow*dt
            )/(volume_next)
        # Using equations from Westhoff et al. 2007, HESS
        depth = volume/(length*width)
        tnew -= dt*self.kh*(temp_curr - Ta)/(self.rho*self.cp*depth) ## sensible heat
        tnew += dt*(1-self.alphaw)*Sin/(self.rho*self.cp*depth) ## solar
        tnew += dt*Lin/(self.rho*self.cp*depth) ## atmospheric
        tnew -= dt*Lout/(self.rho*self.cp*depth) ## back radiation 

        ## Garner, What causes cooling water..., 2014
        tnew += dt*285.9*(0.132 + 0.143*u)*(ea - esat)/(self.rho*self.cp*depth) ## latent heat 

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

class LaxWendroffWesthoff(Temperature):
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
        
        args = ['windspeed','alphaw','eps','rho','cp','kh','sigma','Tgw','temperature']
        for arg in args: setattr(self, arg, kwargs[arg])        
        self.T_soil = self.Tgw
        self.H = 99

        n = 100.
        self.profile = self.Tgw*np.ones(np.shape(np.linspace(0,1,n))) + 273.15
        dx = 1./n


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
        dx = length/len(self.profile)       

        # volumes
        volume = 1e-6*kwargs['volume']

        # time        
        dt = dt*86400


        # energy fluxes
        Lin = kwargs['Lin']
        Sin = kwargs['Sin']               

        # rectangular channel
        # depth = volume/(length*width)

        # power law channel from Finnegan et al 2005
        # alpha = 20 # between cobble and gravel
        # c = (alpha*(alpha+2)**(2/3))**(3/8)
        # width = c*volumetric_discharge**0.5 # eq (5) + eq(7-8)
        # depth = width/alpha

        # triangular channel
        angle = 2*np.pi/20
        width = np.sqrt(2*volume/np.tan(angle)/length)
        
        cross_sectional_area = volume/length


        try: self.H = kwargs['H']
        except: pass

        # now add in various heat fluxes
        back_radiation = -0.96*self.sigma*(temp_curr)**4

        solar = (1-self.alphaw)*Sin

        VTS = 0.9
        land_cover = 0.96*(1-VTS)*self.sigma*(Ta)**4

        K_sed = 3.4 # thermal conductivity of sediment [W/m/C] (Shi et al 1996)
        eta = 0.3 # porosity [] (Westhoff 2007)
        K_w = 0.6 # thermal conductivity of water [W/m/C] (Boyd and Kasper 2003)
        d_soil = 0.071 # thickness of substrate layer [m] (Westhoff 2007)
        rho_sed = 1600 # sediment density [kg/m^3] (Boyd and Kasper 2003)
        rho_w = 1000 # density of water [kg/m^3]
        c_sed = 2219 # specific heat capcity of sediment [J/kg/C] (Boyd and Kasper 2003)
        c_w = 4182 # specific heat capacity of water [J/kg/C] (Boyd and Kasper 2003)
        K_soil = K_sed*(1-eta) + K_w*eta
        rho_soil = rho_sed*(1-eta) + rho_w*eta
        c_soil = c_sed*(1-eta) + c_w*eta

        alluvium_conduction = -K_soil*(self.T_soil - self.Tgw)/d_soil
        conduction = -K_soil*(temp_curr - 273.15 - self.T_soil)/d_soil

        phi_net = solar*0.5/0.5 - conduction + alluvium_conduction
        self.T_soil += phi_net/d_soil/rho_soil/c_soil*dt

        ## general set up 
        rho_w = 1000 # density of water [kg/m^3]
        H = self.H # humidity []
        e_s = 0.61275*np.exp(17.27*(Ta-273.15)/(Ta-273.15+237.3)) # satureated vapor pressure [kPa]
        e_a = H/100*e_s # actual vapor pressure [kPa]

        ## set up for latent
        gamma = 0.66 # psychrometric constant [kPa/C] (Dingman 2002) 
        c_a = 1004 # specific heat capacity of air [J/kg/C] (Dingman 2002)
        rho_a = 1.2 # density of air [kg/m^3] (Williams 2006)
        s = 4100*e_s/(Ta-273.15+237)**2 # slope of saturated vapour pressure curve at Ta [kPa/C]
        L_e = 1000*(2501.4 + temp_curr - 273.15)  # latent heat of vaporation [J/kg]
        r_a = 245/(0.54*self.windspeed + 0.5)  # aerodynamic resistance [m/s]
     
        ## atmospheric
        a_1 = 0.094 # empirical constant [kPa^(-1/2)]
        B_c = .6 # Brunt coefficient []
        eps_atm = 1.1*B_c + a_1*np.sqrt(e_a) # emissivity of atmosphere []
        atmospheric = 0.96*eps_atm*self.sigma*Ta**4

        ## sensible
        elev = 450 # elevation [m] 
        P_A = 101.3 - 0.01055*elev # adiabatic atmospheric pressure [kPa]
        e_s_w = 0.61275*np.exp((17.27*(temp_curr - 273.15))/(237.3+temp_curr-273.15))
        e_a_w = H/100*e_s_w
        B_r = 6.1e-4*P_A * (temp_curr- Ta)/(e_s_w - e_a_w) 

        ## since latent depends on all long wave, we have to put this last
        phi_r = solar + atmospheric + back_radiation + land_cover
        E = s*phi_r/(rho_w*L_e*(s+gamma)) + c_a*rho_a*(e_s-e_a)/(rho_w*L_e*r_a*(s+gamma))
        latent = -rho_w*L_e*E   
        sensible = B_r*latent 

        phi_total = latent + atmospheric + sensible + back_radiation + solar+ land_cover + conduction

        v = volumetric_discharge/cross_sectional_area
        r = v*dt/dx

        # print(dx/v)
        # print(dt)
        # print('')

        upper_bound = self.Tgw + 273.15
        lower_bound = self.profile[-1] + v*dt/dx*(self.profile[-2] - self.profile[-1])

        self.profile[1:-1] = (1-r**2)*self.profile[1:-1] - r/2*(1-r)*self.profile[2:] + r/2*(1+r)*self.profile[:-2] + phi_total*width/cross_sectional_area/self.rho/self.cp*dt
        self.profile[0] = upper_bound
        self.profile[-1] = lower_bound

        self.temperature = self.profile[-1] - 273.15

        return {'conduction':conduction,'land_cover':land_cover,'sensible':sensible,'solar':solar,'atmospheric':atmospheric,'latent':latent,'back_radiation':back_radiation}



class EulerianWesthoff(Temperature):
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
        
        args = ['windspeed','alphaw','eps','rho','cp','kh','sigma','Tgw','temperature']
        for arg in args: setattr(self, arg, kwargs[arg])        
        self.T_soil = self.Tgw
        self.H = 75.
        self.windspeed = 0.

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

        # time        
        dt = dt*86400.

        # energy fluxes
        Lin = kwargs['Lin']
        Sin = kwargs['Sin']               

        # rectangular channel
        # depth = volume/(length*width)

        # power law channel from Finnegan et al 2005
        # alpha = 20 # between cobble and gravel
        # c = (alpha*(alpha+2)**(2/3))**(3/8)
        # width = c*volumetric_discharge**0.5 # eq (5) + eq(7-8)
        # depth = width/alpha

        # triangular channel
        angle = 2*np.pi/20
        width = np.sqrt(2*volume/np.tan(angle)/length)
        
        cross_sectional_area = volume/length
        

        try: self.H = kwargs['H']
        except: pass

        # now add in various heat fluxes
        back_radiation = -0.96*self.sigma*(temp_curr)**4.

        solar = (1-self.alphaw)*Sin

        VTS = 1-self.alphaw
        land_cover = 0.96**2.*(1-VTS)*self.sigma*(Ta)**4.

        K_sed = 3.4 # thermal conductivity of sediment [W/m/C] (Shi et al 1996)
        eta = 0.3 # porosity [] (Westhoff 2007)
        K_w = 0.6 # thermal conductivity of water [W/m/C] (Boyd and Kasper 2003)
        d_soil = 0.071 # thickness of substrate layer [m] (Westhoff 2007)
        rho_sed = 1600. # sediment density [kg/m^3] (Boyd and Kasper 2003)
        rho_w = 1000. # density of water [kg/m^3]
        c_sed = 2219. # specific heat capcity of sediment [J/kg/C] (Boyd and Kasper 2003)
        c_w = 4182. # specific heat capacity of water [J/kg/C] (Boyd and Kasper 2003)
        K_soil = K_sed*(1-eta) + K_w*eta
        rho_soil = rho_sed*(1-eta) + rho_w*eta
        c_soil = c_sed*(1-eta) + c_w*eta

        alluvium_conduction = -K_soil*(self.T_soil - self.Tgw)/d_soil
        conduction = -K_soil*(temp_curr - 273.15 - self.T_soil)/d_soil

        phi_net = solar*0.5/0.5 - conduction + alluvium_conduction
        self.T_soil += phi_net/d_soil/rho_soil/c_soil*dt

        ## general set up 
        rho_w = 1000 # density of water [kg/m^3]
        H = self.H # humidity []
        e_s = 0.61275*np.exp(17.27*(Ta-273.15)/(Ta-273.15+237.3)) # satureated vapor pressure [kPa]
        e_a = H/100.*e_s # actual vapor pressure [kPa]

        ## set up for latent
        gamma = 0.66 # psychrometric constant [kPa/C] (Dingman 2002) 
        c_a = 1004. # specific heat capacity of air [J/kg/C] (Dingman 2002)
        rho_a = 1.2 # density of air [kg/m^3] (Williams 2006)
        s = 4100*e_s/(Ta-273.15+237)**2. # slope of saturated vapour pressure curve at Ta [kPa/C]
        L_e = 1000*(2501.4 + temp_curr - 273.15)  # latent heat of vaporation [J/kg]
        r_a = 245/(0.54*self.windspeed + 0.5)  # aerodynamic resistance [m/s]
     
        ## atmospheric
        a_1 = 0.094 # empirical constant [kPa^(-1/2)]
        B_c = .6 # Brunt coefficient []
        eps_atm = 1.1*B_c + a_1*np.sqrt(e_a) # emissivity of atmosphere []
        atmospheric = 0.96*eps_atm*self.sigma*Ta**4.


        ## sensible
        elev = 450. # elevation [m] 
        P_A = 101.3 - 0.01055*elev # adiabatic atmospheric pressure [kPa]
        e_s_w = 0.61275*np.exp((17.27*(temp_curr - 273.15))/(237.3+temp_curr-273.15))
        e_a_w = H/100.*e_s_w
        B_r = 6.1e-4*P_A * (temp_curr- Ta)/(e_s_w - e_a_w) 

        ## since latent depends on all long wave, we have to put this last
        phi_r = solar + atmospheric + back_radiation + land_cover
        E = s*phi_r/(rho_w*L_e*(s+gamma)) + c_a*rho_a*(e_s-e_a)/(rho_w*L_e*r_a*(s+gamma))
        latent = -rho_w*L_e*E   
        sensible = B_r*latent 

        #need volumetric update back in here...
        flux = ppt*length*width + vol_1 + vol_2 + hillslope_volumetric_discharge 
        flux += hillslope_volumetric_overlandFlow -  volumetric_discharge
        volume_next = volume + flux*dt
        leftover_stream_volume = volume - volumetric_discharge*dt


        tnew = (temp_curr*leftover_stream_volume
            + Ta*dt*ppt*length*width
            + temp_2*vol_2*dt + temp_1*vol_1*dt
            + (self.Tgw + 273.15)*hillslope_volumetric_discharge*dt
            + Ta*hillslope_volumetric_overlandFlow*dt
            )/(volume_next)

        # tnew = temp_curr

        phi_total = latent + atmospheric + sensible + back_radiation + solar+ land_cover + conduction
        tnew += phi_total*length*width/volume/self.rho/self.cp*dt


        self.temperature = tnew - 273.15

        return {'conduction':conduction,'land_cover':land_cover,'sensible':sensible,'solar':solar,'atmospheric':atmospheric,'latent':latent,'back_radiation':back_radiation}

        


class LagrangianSimpleTemperature(Temperature):
    """ Very basic Lagrangian temp model proposed by Yearsley (2009); doi:10.1029/2008WR007629
 
    """
    def __init__(self, rew_id, **kwargs):
        Temperature.__init__(self, rew_id)
        
        args = ['windspeed', 'c1', 'c2', 'alphaw','eps','rho','cp','kh','sigma','Tgw','temperature']
        for arg in args: setattr(self, arg, kwargs[arg])        

        self.T_soil = self.Tgw
        self.H = 20


        self.internalCounter = 0
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

        # volumes
        volume = 1e-6*kwargs['volume']

        # lengths
        length = 0.01*kwargs['length']

        # rectangular channel
        width = 0.01*kwargs['width']
        depth = volume/(length*width)

        # triangular channel
        angle = 2*np.pi/20
        width = np.sqrt(2*volume/np.tan(angle)/length)
        depth = np.tan(angle)*width

        # atmospheric data (water vapor pressure in kPa)
        ea = kwargs['ea']
        
        # energy fluxes
        Lin = kwargs['Lin']
        Sin = kwargs['Sin']

        # try: self.H = kwargs['H']
        # except: pass

        # timestep        
        dt = dt*86400


        #compute distance that outlet parcel travels during this timestep

        u = volumetric_discharge/(volume/length)
        distance = np.min([u*dt, length])

        #interpolate between upstream temp and current outlet temp
        if vol_1+vol_2==0: 
            temp_up =  self.Tgw + 273.15
        else:
            temp_up = vol_1/(vol_1+vol_2)*temp_1 + vol_2/(vol_1+vol_2)*temp_2
        
        temp_start = (temp_curr-temp_up)/length*(length-distance) + temp_up


        # now add in various heat fluxes
        back_radiation = -0.96*self.sigma*(temp_start)**4

        ## sensible and latent from Westhoff 2007
        # D_f = self.alphaw
        # D_diffuse = 0.3 # Westhoff 2007
        # C_s = 0.5 # shadow factor, [0,1]
        # phi_direct = C_s*(1-D_diffuse)*Sin
        # phi_diffuse = D_diffuse*Sin
        # solar = (1-D_f)*(phi_direct + phi_diffuse)
        solar = (1-self.alphaw)*Sin

        VTS = 0.9
        land_cover = 0.96*(1-VTS)*self.sigma*(Ta)**4

        K_sed = 3.4 # thermal conductivity of sediment [W/m/C] (Shi et al 1996)
        eta = 0.3 # porosity [] (Westhoff 2007)
        K_w = 0.6 # thermal conductivity of water [W/m/C] (Boyd and Kasper 2003)
        d_soil = 0.071 # thickness of substrate layer [m] (Westhoff 2007)
        rho_sed = 1600 # sediment density [kg/m^3] (Boyd and Kasper 2003)
        rho_w = 1000 # density of water [kg/m^3]
        c_sed = 2219 # specific heat capcity of sediment [J/kg/C] (Boyd and Kasper 2003)
        c_w = 4182 # specific heat capacity of water [J/kg/C] (Boyd and Kasper 2003)
        K_soil = K_sed*(1-eta) + K_w*eta
        rho_soil = rho_sed*(1-eta) + rho_w*eta
        c_soil = c_sed*(1-eta) + c_w*eta

        alluvium_conduction = -K_soil*(self.T_soil - self.Tgw)/d_soil
        conduction = 0 #-K_soil*(temp_start - 273.15 - self.T_soil)/d_soil

        phi_net = solar*0.5/0.5 - conduction + alluvium_conduction
        self.T_soil += phi_net/d_soil/rho_soil/c_soil*dt


        ## general set up 
        rho_w = 1000 # density of water [kg/m^3]
        H = self.H # humidity []
        e_s = 0.61275*np.exp(17.27*(Ta-273.15)/(Ta-273.15+237.3)) # satureated vapor pressure [kPa]
        e_a = H/100*e_s # actual vapor pressure [kPa]

        ## set up for latent
        gamma = 0.66 # psychrometric constant [kPa/C] (Dingman 2002) 
        c_a = 1004 # specific heat capacity of air [J/kg/C] (Dingman 2002)
        rho_a = 1.2 # density of air [kg/m^3] (Williams 2006)
        s = 4100*e_s/(Ta-273.15+237)**2 # slope of saturated vapour pressure curve at Ta [kPa/C]
        L_e = 1000*(2501.4 + temp_start - 273.15)  # latent heat of vaporation [J/kg]
        r_a = 245/(0.54*self.windspeed + 0.5)  # aerodynamic resistance [m/s]
     
        ## atmospheric
        a_1 = 0.094 # empirical constant [kPa^(-1/2)]
        B_c = .6 # Brunt coefficient []
        eps_atm = 1.1*B_c + a_1*np.sqrt(e_a) # emissivity of atmosphere []
        atmospheric = 0.96*eps_atm*self.sigma*Ta**4


        ## sensible
        elev = 450 # elevation [m] 
        P_A = 101.3 - 0.01055*elev # adiabatic atmospheric pressure [kPa]
        e_s_w = 0.61275*np.exp((17.27*(temp_start - 273.15))/(237.3+temp_start-273.15))
        e_a_w = H/100*e_s_w
        B_r = 6.1e-4*P_A * (temp_start- Ta)/(e_s_w - e_a_w) 

        ## since latent depends on all long wave, we have to put this last
        phi_r = solar + atmospheric + back_radiation + land_cover
        E = s*phi_r/(rho_w*L_e*(s+gamma)) + c_a*rho_a*(e_s-e_a)/(rho_w*L_e*r_a*(s+gamma))
        latent = -rho_w*L_e*E   
        sensible = B_r*latent 

        #sensible = -(self.kh*(temp_start - Ta))
        #esat = 0.611*np.exp(2.5*10**6/461.0*(1/273.2 - 1/temp_start))        
        #latent = (285.9*(0.132 + 0.143*self.windspeed)*(ea - esat))

        temp_start += (conduction + land_cover + sensible + solar + atmospheric + latent + back_radiation)/(depth*self.rho*self.cp)*dt

        # combine with incoming water fluxes
        volume_end = volume/length + width*ppt*dt + 1/length*hillslope_volumetric_discharge*dt
        tnew = (volume/length*temp_start + Ta*width*ppt*dt + (self.Tgw + 273.15)/length*hillslope_volumetric_discharge*dt)/volume_end

        self.startContrib = volume/length*temp_start/volume_end/tnew
        self.hsdContrib = (self.Tgw + 273.15)/length*hillslope_volumetric_discharge*dt/volume_end/tnew
        self.pptContrib = Ta*width*ppt*dt/volume_end/tnew

        self.temperature = tnew - 273.15

        return {'conduction':conduction,'land_cover':land_cover,'sensible':sensible,'solar':solar,'atmospheric':atmospheric,'latent':latent,'back_radiation':back_radiation}


class LagrangianSimpleTemperatureTriangular(Temperature):
    """ 
 
    """
    def __init__(self, rew_id, **kwargs):
        Temperature.__init__(self, rew_id)
        
        args = ['angle','upstream_area', 'mannings_n', 'gradient', 'alphaw','eps','rho','cp','kh','sigma','Tgw_amplitude', 'Tgw_phase', 'Tgw_offset', 'Tgw_sd', 'temperature']
        for arg in args: setattr(self, arg, kwargs[arg])        

        self.internalCounter = 0

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
        # day of year
        doy = kwargs['doy']

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
        Tgw = self.Tgw_offset + self.Tgw_amplitude*np.exp(-(doy-self.Tgw_phase)**2/(2*self.Tgw_sd**2)) 


        # timestep        
        dt = dt*86400


        #compute distance that outlet parcel travels during this timestep
        distance = np.min([u*dt, length])

        #interpolate between upstream temp and current outlet temp
        if vol_1+vol_2==0: 
            temp_up =  Tgw + 273.15
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

        added_groundwater = hillslope_volumetric_discharge/length*dt
        added_overlandFlow = hillslope_volumetric_overlandFlow/length*dt
        vol_initial = width*depth

        tnew -= dt*self.kh*(temp_start - Ta)/(self.rho*self.cp*depth) ## sensible heat
        tnew += dt*(1-self.alphaw)*Sin/(self.rho*self.cp*depth) ## solar
        tnew += dt*Lin/(self.rho*self.cp*depth) ## atmospheric
        tnew -= dt*back_radiation/(self.rho*self.cp*depth) ## back radiation 

        ## Garner, What causes cooling water..., 2014
        tnew += dt*285.9*(0.132 + 0.143*0.5)*(ea - esat)/(self.rho*self.cp*depth) ## latent heat

        tnew = (tnew*vol_initial + (Tgw+273.15)*added_groundwater + Ta*added_overlandFlow)/(vol_initial + added_groundwater + added_overlandFlow)

        self.temperature = tnew - 273.15

        return (depth, u)

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
        hc = 1.5*10**6*1.74*10**(-6)*(1 + 0.72*kwargs['va'])


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


class LagrangianSimpleTemperatureTriangularHeatedGW(Temperature):
    """ 
 
    """
    def __init__(self, rew_id, **kwargs):
        Temperature.__init__(self, rew_id)
        
        args = ['angle','upstream_area', 'mannings_n', 'gradient', 'alphaw','eps','rho','cp','kh','sigma','temperature', 'kf', 'tau0', 'ktau', 'Tgw_offset']
        for arg in args: setattr(self, arg, kwargs[arg])        

        self.internalCounter = 0

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
        # day of year
        doy = kwargs['doy']

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
        Tgw = (self.Tgw_offset+273.15)*np.exp(-self.kf*tau) + Ta_mean*(1-np.exp(-self.kf*tau)) - 273.15


        # timestep        
        dt = dt*86400


        #compute distance that outlet parcel travels during this timestep
        distance = np.min([u*dt, length])

        #interpolate between upstream temp and current outlet temp
        if vol_1+vol_2==0: 
            temp_up =  Tgw + 273.15
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

        tnew -= dt*self.kh*(temp_start - Ta)/(self.rho*self.cp*depth) ## sensible heat
        tnew += dt*(1-self.alphaw)*Sin/(self.rho*self.cp*depth) ## solar
        tnew += dt*Lin/(self.rho*self.cp*depth) ## atmospheric
        tnew += dt*land_cover/(self.rho*self.cp*depth) ## atmospheric
        tnew -= dt*back_radiation/(self.rho*self.cp*depth) ## back radiation 

        ## Garner, What causes cooling water..., 2014
        # tnew += dt*285.9*(0.132 + 0.143*kwargs['va'])*(ea - esat)/(self.rho*self.cp*depth) ## latent heat

        tnew = (tnew*vol_initial + (Tgw+273.15)*added_groundwater + Ta*added_overlandFlow)/(vol_initial + added_groundwater + added_overlandFlow)

        self.temperature = tnew - 273.15

        return (depth, u, Tgw, hillslope_volumetric_discharge)
