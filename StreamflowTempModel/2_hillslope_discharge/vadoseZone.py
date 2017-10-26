import numpy as np

class VadoseZone:
    """ Abstract base class for all vadose zone models. 
    """

    def __init__(self): pass

    def update(self,**kwargs): 
        """ Update vadose zone stocks and compute fluxes.
    
        """
        raise NameError('update')

class LaioVadoseZone(VadoseZone):
    """ Vadose zone model based on Laio et al (doi:10.1016/S0309-1708(01)00005-7)

    Public attributes:
         - storageVZ [cm]: vadose zone storage stock
         - vi [cm/d]: rainfall intercepted by vegetation
         - zr [cm]: root zone depth
         - n []: porosity
         - sh []: hygroscopic point
         - sw []: wilting point
         - sfc []: field capacity
         - sstar []:
         - ew [cm/d]: evapotranspiration at wiling point
         - emax [cm/d]: maximum evapotranspiration
         - ks [cm/d]: saturated hydraulic conducitvity
         - b []: 
    """

    #NATE DAWG, I THINK WE SHOULD MAKE VEG INTERCEPTION AT THE REW LEVEL. OTHERWISE, WE'D NEED TO SPREAD vi 
    #ACROSS EACH DAY; right now, if dt is very small, vi will always be larger than rain. Either that, 
    #or we need to think of vi and ppt as

    #also, in your groundwater zone, i think your output discharge is in units of depth; to get in cm/day (like in your plot)
    #you'd need to divide those outputs by dt? 
    def __init__(self, **kwargs):
        args = ['storageVZ','vi','zr','n','sh','sw','sfc','sstar','ew','emax','ks','b']
        for arg in args: setattr(self, arg, kwargs[arg])

        # main external variables
        self.leakage        = 0         # [cm/d]
        self.ET             = 0         # [cm/d]
        self.overlandFlow   = 0         # [cm/d]
        
        # auxiliary external variables
        self.asd    = self.zr*self.n
        
        # internal variables
        self._etaw  = self.ew/self.asd
        self._eta   = self.emax/self.asd
        self._m     = self.ks/(self.asd*(np.exp(self.b*(1-self.sfc))-1))
    
    def update(self,dt,**kwargs):
        """ Update vadose zone stocks and compute fluxes.
    
        Args:
            - dt (float): time step
            
        Kwargs (dict), with keys:
            - ppt (float): precipitation flux [L / T]
            - pet (float): potential evapotranspiration flux [L / T] -- if not specified, set to LaioVadoseZone::emax
        Returns: 
            - fluxes (dict): dictionary of fluxes, with keys [ET, leakage, overlandFlow, interceptedRainfall]
    
        """
        
        ppt = kwargs['ppt']
        pet = kwargs.get('pet',self._eta)
        
        self.interceptedRainfall = np.min((self.vi,ppt))    # [cm/d]    
        R = np.max((ppt - self.vi,0))                       # [cm/d]
        
        s = self.storageVZ/self.asd                           # []
        
        # et and lkg: [1/d]
        if s < self.sh: 
            et = 0
        elif s <= self.sw: 
            et = self._etaw*(s - self.sh)/(self.sw - self.sh)
        elif s <= self.sstar: 
            et = self._etaw + (pet - self._etaw)*(s - self.sw)/(self.sstar - self.sw)
        else: 
            et = pet

        if s < self.sfc: lkg = 0
        else: lkg = self._m*(np.exp(self.b*(s - self.sfc))-1)
  
        self.ET             = et*self.asd                               # [cm/d]
        self.leakage        = lkg*self.asd                              # [cm/d]
        self.storageVZ        += (R - self.ET - self.leakage)*dt          # [cm]
        self.overlandFlow   = np.max((self.storageVZ - self.asd,0))/dt    # [cm/d]
        self.storageVZ        = np.min((self.storageVZ,self.asd))           # [cm]

        return {'ET':self.ET, 'leakage':self.leakage, 'overlandFlow':self.overlandFlow,'interceptedRainfall':self.interceptedRainfall}
        
    def _closed_form(self,t,s0):
        """ Validation method.
        
        """
        tsfc    = 1/(self.b*(self._m-self._eta))*(self.b*(self.sfc-s0) + np.log((self._eta - self._m + self._m*np.exp(self.b*(s0-self.sfc)))/self._eta))
        tsstar  = (self.sfc - self.sstar)/self._eta + tsfc
        tsw     = (self.sstar-self.sw)/(self._eta - self._etaw)*np.log(self.emax/self._etaw) + tsstar


        if t < tsfc: return s0 - (1/self.b)*np.log(((self._eta - self._m + self._m*np.exp(self.b*(s0 - self.sfc)))*np.exp(self.b*(self._eta-self._m)*t) - self._m*np.exp(self.b*(s0-self.sfc )))/(self._eta-self._m))
        if t < tsstar: return self.sfc - self._eta*(t - tsfc)
        if t < tsw: return self.sw + (self.sstar - self.sw)*(self._eta/(self._eta - self._etaw)*np.exp(-(self._eta-self._etaw)/(self.sstar - self.sw)*(t - tsstar)) - self._etaw/(self._eta-self._etaw))
        else: return self.sh + (self.sw - self.sh)*np.exp(-self._etaw/(self.sw-self.sh)*(t-tsw))

class PorporatoPreferentialVadoseZone(VadoseZone): 
    """ Vadose zone model based on Porporato 
    
    Public attributes:
        - storageVZ (float): [L] current vadose zone storage stock
        - zr (float): [L] rooting zone depth
        - s0 (float): [] wilting point
        - st (float): [] field capacity
        - n (float): [] porosity 
        - alpha (float): [] fraction of incoming rainfall as preferential flow 
    
    """
    def __init__(self, **kwargs):        
    
        args = ['storageVZ','zr','s0','st','n', 'alpha', 'eta']
        for arg in args: setattr(self, arg, kwargs[arg])

        # main external variables
        self.leakage        = 0             # [cm/day]
        self.ET             = 0             # [cm/day]
        self.overlandFlow   = 0

    def update(self,dt,**kwargs):
        """ Update vadose zone stocks and compute fluxes.
    
        Args:
            - dt (float): time step
            - ppt (float): precipitation flux [L / T]
            - pet (float): potential evapotranspiration flux [L / T] -- if not specified, set to LaioVadoseZone::emax

        Returns: 
            - fluxes (dict): dictionary of fluxes, with keys [ET, leakage]
         """
        
        ppt = kwargs['ppt']
        pet = kwargs['pet']
        s = self.storageVZ/(self.n*self.zr)
        x = (s-self.s0)/(self.st - self.s0)
        bypass = ppt*self.alpha*x
        ppt = ppt - bypass

        self.ET = self.eta*pet*x
        s += ppt*dt/(self.n*self.zr) - self.ET*dt/(self.n*self.zr)
        leakage = np.max([s - self.st, 0])*self.n*self.zr/dt
        s = np.max([np.min([s, self.st]), self.s0])
        self.storageVZ = s*self.n*self.zr
        self.leakage = leakage + bypass
        return {'ET':self.ET, 'leakage':self.leakage, 'overlandFlow':self.overlandFlow}



class PorporatoVadoseZone(VadoseZone): 
    """ Vadose zone model based on Porporato 
    
    Public attributes:
        - storageVZ (float): [L] current vadose zone storage stock
        - zr (float): [L] rooting zone depth
        - s0 (float): [] wilting point
        - st (float): [] field capacity
        - n (float): [] porosity 
    
    """
    def __init__(self, **kwargs):        
    
        args = ['storageVZ','zr','s0','st','n', 'eta']
        for arg in args: setattr(self, arg, kwargs[arg])

        # main external variables
        self.leakage        = 0             # [cm/day]
        self.ET             = 0             # [cm/day]
        self.overlandFlow   = 0

    def update(self,dt,**kwargs):
        """ Update vadose zone stocks and compute fluxes.
    
        Args:
            - dt (float): time step
            - ppt (float): precipitation flux [L / T]
            - pet (float): potential evapotranspiration flux [L / T] -- if not specified, set to LaioVadoseZone::emax

        Returns: 
            - fluxes (dict): dictionary of fluxes, with keys [ET, leakage]
         """
        
        ppt = kwargs['ppt']
        pet = kwargs['pet']
        s = self.storageVZ/(self.n*self.zr)
        x = (s-self.s0)/(self.st - self.s0)

        self.ET = self.eta*pet*x
        s += ppt*dt/(self.n*self.zr) - self.ET*dt/(self.n*self.zr)
        leakage = np.max([s - self.st, 0])*self.n*self.zr/dt
        s = np.max([np.min([s, self.st]), self.s0])
        self.storageVZ = s*self.n*self.zr
        self.leakage = leakage
        return {'ET':self.ET, 'leakage':self.leakage, 'overlandFlow':self.overlandFlow}



# class MelangeVadoseZone(VadoseZone): 
#     """ Vadose zone model based on Rodriguez Iturbe (1999)
    
#     Public attributes:
#         - eta
#         - s1
#         - sstar 
#         - k1 
#         - k12
#         - n 
#         - zr 
    
#     """
#     def __init__(self, **kwargs):        
    
#         args = ['storageVZ', 'eta', 's1', 'sstar', 'k1', 'k12', 'n', 'zr']
#         for arg in args: setattr(self, arg, kwargs[arg])

#         # main external variables
#         self.leakage        = 0             # [cm/day]
#         self.ET             = 0             # [cm/day]
#         self.overlandFlow   = 0
#         self.overlandRes    = 0 
#         self.soilRes        = self.storageVZ
    
#     def update(self,dt,**kwargs):
#         """ Update vadose zone stocks and compute fluxes.
    
#         Args:
#             - dt (float): time step
#             - ppt (float): precipitation flux [L / T]
#             - pet (float): potential evapotranspiration flux [L / T] -- if not specified, set to LaioVadoseZone::emax

#         Returns: 
#             - fluxes (dict): dictionary of fluxes, with keys [ET, leakage]
#          """
        
#         ppt = kwargs['ppt']
#         pet = kwargs['pet']

#         #convert current storage to normalized relative soil moisture
#         s = self.soilRes/(self.n*self.zr)
#         self.ET = self.eta*pet if (s > self.sstar) else s*(pet*self.eta)/self.sstar

#         s += ppt*dt/(self.n*self.zr) - self.ET*dt/(self.n*self.zr)

#         #anything in excess of field capacity is drained to groundwater (slow zone) 
#         #or to surface overland flow; k12 splits leakage between these reservoirs
#         leakage = np.max([s - self.s1, 0])*self.n*self.zr/dt
#         s = np.min([s, self.s1])

#         self.soilRes = s*self.n*self.zr

#         # get output from overland flow reservoir
#         self.overlandFlow = self.overlandRes*self.k1

#         # get input to overland flow reservoir
#         overlandResInput = leakage*(1-self.k12)

#         # update overland flow reservoir
#         self.overlandRes += overlandResInput*dt - self.overlandFlow*dt 

#         self.leakage = leakage*self.k12

#         self.storageVZ = self.overlandRes + self.soilRes

#         return {'ET':self.ET, 'leakage':self.leakage, 'overlandFlow':self.overlandFlow}


class SimpleRockMoistureZone(VadoseZone): 
    """ Two layer vadose zone model. Each layer is treated as a porporato type vadose zone. 
    
    Public attributes:
        - nS (float): [] porosity of soil layer
        - nR (float): [] porosity of rock moisture layer
        - s0R (float): [] water content at which ET = 0 in rock moisture zone
        - s0S (float): [] water content at which ET = 0 in soil moisture zone
        - stR (float): [] water content at which ET = PET in rock moisture zone
        - stS (float): [] water content at which ET = PET in soil moisture zone
        - zrR (float): [L] thickness of rock moisture zone
        - zrS (float): [L] thickness of soil moisture zone
        - f (float): [] fraction of roots in soil layer (i.e., fraction of PET apportioned to soil layer)
        - storageR (float): [L] storage in rock moisture zone
        - storageS (float): [L] storage in soil moisture zone
    
    """
    def __init__(self, **kwargs):        
    
        args = ['nS','nR','s0R','s0S','stR','stS','zrR','zrS','f','storageR','storageS']
        for arg in args: setattr(self, arg, kwargs[arg])

        # main external variables
        self.leakage         = 0            # [cm/day]
        self.ETR             = 0            # [cm/day]
        self.ETS             = 0            # [cm/day]
        self.ET              = 0            # [cm/day]
        self.overlandFlow    = 0 
        self.storageVZ       = self.storageS + self.storageR

    def update(self,dt,**kwargs):
        """ Update vadose zone stocks and compute fluxes.
    
        Args:
            - dt (float): time step
            - ppt (float): precipitation flux [L / T]
            - pet (float): potential evapotranspiration flux [L / T] -- if not specified, set to ::emax

        Returns: 
            - fluxes (dict): dictionary of fluxes, with keys [ET, leakage]
         """
        
        ppt = kwargs['ppt']
        pet = kwargs['pet']

        #first compute soil moisture zone
        sS = self.storageS/(self.nS*self.zrS)

        if (sS <= self.s0S):
            self.ETS = 0 
        elif (sS <= self.stS):
            self.ETS = (self.f)*pet*(sS-self.s0S)/(self.stS-self.s0S)
        else: 
            self.ETS = (self.f)*pet    

        sS += ppt*dt/(self.nS*self.zrS) - self.ETS*dt/(self.nS*self.zrS)

        #anything in excess of stS is drained to rock moisture
        soilLeakage = np.max([sS - self.stS, 0])*self.nS*self.zrS/dt
        sS = np.min([sS, self.stS])
        self.storageS = sS*self.nS*self.zrS


        #convert current storage to normalized relative value
        sR = self.storageR/(self.nR*self.zrR)
        if (sR <= self.s0R):
            self.ETR = 0
        elif (sR <= self.stR):
            self.ETR = (1-self.f)*pet*(sR-self.s0R)/(self.stR-self.s0R)
        else: 
            self.ETR = (1-self.f)*pet

        sR += soilLeakage*dt/(self.nR*self.zrR) - self.ETR*dt/(self.nR*self.zrR)

        #anything in excess of stS is drained to rock moisture
        self.leakage = np.max([sR - self.stR, 0])*self.nR*self.zrR/dt
        sR = np.min([sR, self.stR])
        self.storageR = sR*self.nR*self.zrR

        self.ET = self.ETS + self.ETR

        self.storageVZ = self.storageS + self.storageR
        
        return {'ET':self.ET, 'leakage':self.leakage, 'overlandFlow':self.overlandFlow}


class PreferentialRockMoistureZone(VadoseZone): 
    """ Two layer vadose zone model. Each layer is treated as a porporato type vadose zone. 
    
    Public attributes:
        - nS (float): [] porosity of soil layer
        - nR (float): [] porosity of rock moisture layer
        - s0R (float): [] water content at which ET = 0 in rock moisture zone
        - s0S (float): [] water content at which ET = 0 in soil moisture zone
        - stR (float): [] water content at which ET = PET in rock moisture zone
        - stS (float): [] water content at which ET = PET in soil moisture zone
        - zrR (float): [L] thickness of rock moisture zone
        - zrS (float): [L] thickness of soil moisture zone
        - alpha (float): [] fraction of water that exits soil zone and is preferentially routed to groundwater
        - f (float): [] fraction of roots in soil layer (i.e., fraction of PET apportioned to soil layer)
        - storageR (float): [L] storage in rock moisture zone
        - storageS (float): [L] storage in soil moisture zone
    
    """
    def __init__(self, **kwargs):        
    
        args = ['nS','nR','s0R','s0S','stR','stS','zrR','zrS','f','alpha','eta','storageR','storageS']
        for arg in args: setattr(self, arg, kwargs[arg])

        # main external variables
        self.leakage         = 0            # [cm/day]
        self.ETR             = 0            # [cm/day]
        self.ETS             = 0            # [cm/day]
        self.ET              = 0            # [cm/day]
        self.overlandFlow    = 0 
        self.storageVZ       = self.storageS + self.storageR

    def update(self,dt,**kwargs):
        """ Update vadose zone stocks and compute fluxes.
    
        Args:
            - dt (float): time step
            - ppt (float): precipitation flux [L / T]
            - pet (float): potential evapotranspiration flux [L / T] -- if not specified, set to ::emax

        Returns: 
            - fluxes (dict): dictionary of fluxes, with keys [ET, leakage]
         """
        
        ppt = kwargs['ppt']
        pet = self.eta*kwargs['pet']

        #convert current storage to normalized relative value
        sR = self.storageR/(self.nR*self.zrR)
        xR = (sR - self.s0R)/(self.stR - self.s0R)
        sS = self.storageS/(self.nS*self.zrS)

        # frac = self.storageVZ/(self.nS*self.zrS + self.nR*self.zrR)

        if (sS <= self.s0S):
            self.ETS = 0 
        elif (sS <= self.stS):
            self.ETS = (self.f)*pet*(sS-self.s0S)/(self.stS-self.s0S)
        else: 
            self.ETS = (self.f)*pet    

        sS += ppt*dt/(self.nS*self.zrS) - self.ETS*dt/(self.nS*self.zrS)

        #anything in excess of stS is drained to rock moisture
        soilLeakage = np.max([sS - self.stS, 0])*self.nS*self.zrS/dt
        sS = np.min([sS, self.stS])
        self.storageS = sS*self.nS*self.zrS


        
        if (sR <= self.s0R):
            self.ETR = 0
        elif (sR <= self.stR):
            self.ETR = (1-self.f)*pet*(sR-self.s0R)/(self.stR-self.s0R)
        else: 
            self.ETR = (1-self.f)*pet

        # soil leakage is partitioned into bypass flow, and storage in rock moisture
        # bypass flow is a linear function of rock moisture storage
        bypass = soilLeakage*self.alpha*xR
        leakageToMatrix = soilLeakage - bypass

        sR += leakageToMatrix*dt/(self.nR*self.zrR) - self.ETR*dt/(self.nR*self.zrR)

        #anything in excess of stR is drained
        self.leakage = np.max([sR - self.stR, 0])*self.nR*self.zrR/dt
        # add water preferentially routed through rock moisture zone
        self.leakage += bypass

        sR = np.min([sR, self.stR])
        self.storageR = sR*self.nR*self.zrR
        self.ET = self.ETS + self.ETR
        self.storageVZ = self.storageS + self.storageR
        
        return {'ET':self.ET, 'leakage':self.leakage, 'overlandFlow':self.overlandFlow}


