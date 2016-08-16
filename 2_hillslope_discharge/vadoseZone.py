import numpy as np

class VadoseZone:
    """ Abstract base class for all vadose zone models. 

    Public attributes: 
        - rew: :class:`REW` instance to which the vadose belongs
    """

    def __init__(self,rew): self.rew = rew     

    def update(self,**kwargs): 
        """ Update vadose zone stocks and compute fluxes.
    
        """
        raise NameError('update')

class LaioVadoseZone(VadoseZone):
    """ Vadose zone model based on Laio et al (doi:10.1016/S0309-1708(01)00005-7)

    Public attributes:
         - rew: :class:`REW` instance to which the vadose zone belongs
         - storage [cm]: vadose zone storage stock
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
    def __init__(self, rew, storage = 0, vi = 0.2, zr = 75, n = 0.45, sh = 0.19, sw = 0.24,sfc = 0.65,sstar = 0.57,ew = 0.01,emax = 0.5,ks = 20,b = 14.8):
        # make sure we get the attrs that ANY vadose zone model must have
        VadoseZone.__init__(self,rew)
        
        # right now these are all set to the defaults, 
        # but we could easily set these values based on REW properties 

        # main external variables
        self.storage        = storage   # [cm]
        self.leakage        = 0         # [cm/d]
        self.ET             = 0         # [cm/d]
        self.overlandFlow   = 0         # [cm/d]
                
        # auxiliary external variables
        self.vi     = vi
        self.sh     = sh    
        self.sw     = sw    
        self.ew     = ew 
        self.emax   = emax
        self.ks     = ks
        self.sfc    = sfc
        self.b      = b 
        self.sstar  = sstar
        self.asd    = zr*n
        
        # internal variables
        self._etaw  = ew/self.asd
        self._eta   = emax/self.asd
        self._m     = ks/(self.asd*(np.exp(b*(1-sfc))-1))
        self._s     = self.storage/self.asd
    
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
        
        s = self.storage/self.asd                           # []
        
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
        self.storage        += (R - self.ET - self.leakage)*dt          # [cm]
        self.overlandFlow   = np.max((self.storage - self.asd,0))/dt    # [cm/d]
        self.storage        = np.min((self.storage,self.asd))           # [cm]

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

class PorporatoVadoseZone(VadoseZone): 
    """ Vadose zone model based on Porporato 
    
    Public attributes:
    	- rew: :class:`REW` instance to which the vadose belongs
    	- storage (float): [L] current vadose zone storage stock
    	- zr (float): [L] rooting zone depth
    	- sw (float): [] wilting point
    	- emax (float): [L/T] maximum ET
    	- sfc (float): [] field capacity
    	- n (float): [] porosity 
    
    """
    def __init__(self, rew, storage = 0, zr = 75, sw = 0.19, emax = 0.5, sfc = 0.65, n = 0.45):
        
        # make sure we get the attrs that ANY vadose zone model must have
        VadoseZone.__init__(self,rew)
        
        # right now these are all set to the defaults, 
        # but we could easily set these values based on REW properties 

        # main external variables
        self.storage        = storage       # [cm]
        self.leakage        = 0             # [cm/day]
        self.ET             = 0             # [cm/day]
                
        # auxiliary external variables
        self.zr     = zr 
        self.sw     = sw    
        self.emax   = emax
        self.sfc    = sfc
        self.n      = n
        
    
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
        pet = kwargs.get('pet',self.emax)

        #convert current storage to normalized relative soil moisture
        s = self.storage/(self.n*self.zr)
        self.ET = 0 if (s <= self.sw) else pet*(s-self.sw)/(self.sfc-self.sw)

        s += ppt*dt/(self.n*self.zr) - self.ET*dt/(self.n*self.zr)

        #anything in excess of field capacity is drained
        self.leakage = np.max([s - self.sfc, 0])*self.n*self.zr/dt
        s = np.min([s, self.sfc])
        self.storage = s*self.n*self.zr

        return {'ET':self.ET, 'leakage':self.leakage}



