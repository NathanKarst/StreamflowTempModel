import numpy as np

class VadoseZone:
    def __init__(self): pass
                
    def computeFluxes(self): raise NameError('computeFluxes')
          
class LaioVadoseZone(VadoseZone):
    
    def __init__(self,
                    rew,
                    s   = 0,    # [] soil moisture 
                    vi  = 0.2,  # [cm] rainfall intercepted by veg. (default: trees; 0.05 cm for grass)
                    asd = 1,    # [cm] active soil depth
                    sh  = 0.19, # [] hygroscopic point (default: loam)
                    sw  = 0.24, # [] wilting point (default: loam)
                    ew  = 0.01, # [cm/d] ET at wilting point
                    emax = 0.5, # [cm/d] maximum ET (default: tree-dominated; hot growing season)
                    ks  = 2,    # [cm/d] saturated hydraulic conductivity (default: loam)
                    sfc = 0.65, # [] field capacity (default: loam)
                    sstar =0.57,# []  
                    b   = 14.8  # [] (default: loam)
                    ):
        
        # make sure we get the attrs that ANY vadose zone model must have
        VadoseZone.__init__(self)  
        self.rew = rew
        self.vi     = vi
        self.asd    = asd
        self.sh     = sh    
        self.sw     = sw    
        self.s      = s     
        self.ew     = ew/self.asd    
        self.emax   = emax/self.asd
        self.ks     = ks
        self.sfc    = sfc
        self.b      = b 
        self.sstar  = sstar
        self.m      = self.ks/(self.asd*(np.exp(self.b*(1-self.sfc))-1)) # used in leakage eqn.
    
    def computeFluxes(self,time,**inputs):
        ## dimensionless!!            
       #  inputs['temp']   # this is how you access by default 
        interceptedRainfall = np.max((self.rew.ppt[time] - self.vi,0))/self.asd
        
        if self.s > self.sh and self.s <= self.sw: self.et = self.ew*(self.s - self.sh)/(self.sw - self.sh)
        elif self.s > self.sw and self.s <= self.sstar: self.et = self.ew + (self.emax - self.ew)*(self.s-self.sw)/(self.sstar - self.sw)
        else: self.et = self.emax

        if self.s < self.sfc: self.lkg = 0
        else: self.lkg = self.m*(np.exp(self.b*(self.s - self.sfc))-1)
        
        return((interceptedRainfall,self.et,self.lkg))
        