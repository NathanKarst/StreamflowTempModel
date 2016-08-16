class GroundwaterZone:
    """ Abstract base class for all groundwater zone models. 

    Public attributes: 
        - rew: :class:`REW` instance

    Public methods:
        - update: update state of groundwater zone stocks and return dictionary of fluxes
    """

    def __init__(self,rew): self.rew = rew

    def update(self): 
        """ Update state of groundwater zone stocks.
        """
        raise NameError('update')
        
class LinearReservoir(GroundwaterZone):
    """ Simple linear reservoir model.
    
    Public attributes: 
        - groundwater (float): [L] initial condition for groundwater stock
        - k (float): [1/T] recession parameter
        - discharge (float): [L/T] current discharge flux
    """
    def __init__(self,rew,groundwater=0,k=.2):
        GroundwaterZone.__init__(self,rew)
        self.groundwater = groundwater      # [cm]
        self.k = k                          # [1/d]
        self.discharge = 0                  # [cm/d]
        
    def update(self,dt,**kwargs):
        """ Update state of groundwater zone stocks.
        
        Args:
            - dt (float): [T] time step
            
        Kwargs (dict) has keys:
            - leakage (float): i[L/T] incoming leakage flux from vadose zone
        
        Returns:
            - fluxes (dict): dictionary of fluxes [L/T], with keys [discharge]
        """
        leakage = kwargs['leakage']
        
        self.discharge = self.k*self.groundwater                # [cm/d]
        self.groundwater += (-self.discharge + leakage)*dt      # [cm]
        
        return {'discharge':self.discharge}
        
class NonlinearReservoir(GroundwaterZone):
    """ Nonlinear reservoir model.
    
    Public attributes: 
        - a (float): powerlaw scale parameter
        - b (float): powerlaw exponent
        - discharge: [L/T] current discharge flux
        - groundwater: [L] current groundwater stock
    
    Public methods:
        - update:  (here, discharge only)
    """
    def __init__(self,rew,groundwater=0,a=.5,b=1):
        GroundwaterZone.__init__(self,rew)
        self.groundwater = groundwater
        self.a = a
        self.b = b
        self.discharge = self.a*self.groundwater**self.b
        
    def update(self,dt,**kwargs):
        """ Update state of groundwater zone stocks.
        
        Args:
            - dt (float): [T] time step
            
        Kwargs (dict) has keys:
            - leakage (float): i[L/T] incoming leakage flux from vadose zone
        
        Returns:
            - fluxes (dict): dictionary of fluxes [L/T], with keys [discharge]
        """
        # if isinstance(leakage,np.ndarray):
        #     discharge = np.zeros(np.size(leakage))
        #     groundwater = np.zeros(np.size(leakage))
        #     discharge[0] = self.a*self.groundwater**self.b
        #     groundwater[0] = (-self.discharge+leakage[0])*dt+self.groundwater
        #     for i in xrange(1,np.size(leakage)):
        #         discharge[i] = self.a*groundwater[i]**self.b
        #         groundwater[i] = (-discharge[i]+leakage[i])*dt+groundwater[i-1]

        #     self.groundwater = groundwater[-1]
        #     self.discharge = discharge[-1]
        #     return discharge, groundwater

        leakage = kwargs['leakage']

        # else:
        self.discharge = self.a*self.groundwater**self.b
        self.groundwater += (-self.discharge + leakage)*dt
#         self.overlandFlow = np.max((self.groundwater - self.groundwaterMax,0))
#         self.groundwater = np.min((self.groundwater,self.groundwaterMax))
        return {'discharge':self.discharge}

class TwoLinearReservoir(GroundwaterZone):
    """ Two linear reservoir model.
    
    Res. 1 and res. 2 are in series. Res. 1 drains to both the stream and res. 2; res. 2 drains only to stream.
    
    Public attributes: 
        - k1 (float): [1/T] rate at which res. 1 drains to stream
        - res1 (float): [L] current groundwater stock in res. 1
        - k2 (float): [1/T] rate at which res. 2 drains to stream
        - res2 (float): [L] current groundwater stock in res. 2
        - k12 (float): [1/T] rate at which res. 1 drains to res. 2
    """
    def __init__(self, rew, groundwater = 0, res1 = 0, res2 = 0, k1 = 0.5, k2 = 0.2, k12 = 0.1):
        GroundwaterZone.__init__(self,rew)
        self.res1 = res1
        self.res2 = res2
        self.groundwater = res1 + res2
        self.k1 = k1
        self.k2 = k2
        self.k12 = k12
        self.discharge = 0
        
    def update(self,dt, **kwargs):
        """ Update state of groundwater zone stocks.
        
        Res. 1 drains to both the stream and res. 2; res. 2 drains only to stream.
        
        Args:
            - dt (float): [T] time step
            
        Kwargs (dict) has keys:
            - leakage (float): i[L/T] incoming leakage flux from vadose zone
        
        Returns:
            - fluxes (dict): dictionary of fluxes [L/T], with keys [discharge]
        """
        leakage = kwargs['leakage']
        
        self.discharge = self.k2*self.res2 + self.k1*self.res1
        self.groundwater += (leakage - self.discharge)*dt
        self.res2 += self.k12*self.res1*dt - self.k2*self.res2*dt
        self.res1 += leakage*dt - self.k12*self.res1*dt - self.k1*self.res1*dt
        
        return {'discharge':self.discharge}

class TwoParallelLinearReservoir(GroundwaterZone):
    """ Two linear reservoir model.
    
    Res. 1 and res. 2 are in parallel; both drain directly to stream.
    
    Public attributes: 
        - k1 (float): [1/T] rate at which res. 1 drains to stream
        - res1 (float): [L] current groundwater stock in res. 1
        - k2 (float): [1/T] rate at which res. 2 drains to stream
        - res2 (float): [L] current groundwater stock in res. 2
        - f1 (float): [] fraction of incoming flux apportioned to res. 1
    """
    def __init__(self, rew, groundwater = 0, res1 = 0, res2 = 0, k1 = 0.5, k2 = 0.2, f1 = 0.5):
        GroundwaterZone.__init__(self,rew)
        self.res1 = res1
        self.res2 = res2
        self.groundwater = res1 + res2
        self.k1 = k1
        self.k2 = k2
        self.f1 = f1
        self.discharge = 0
        
    def update(self,dt,**kwargs):
        """ Update state of groundwater zone stocks.
        
        Res. 1 drains to both the stream and res. 2; res. 2 drains only to stream.
        
        Args:
            - dt (float): [T] time step
            
        Kwargs (dict) has keys:
            - leakage (float): i[L/T] incoming leakage flux from vadose zone
        
        Returns:
            - fluxes (dict): dictionary of fluxes [L/T], with keys [discharge]
        """
        leakage = kwargs['leakage']
        
        self.discharge = self.k2*self.res2 + self.k1*self.res1
        self.groundwater += (leakage - self.discharge)*dt
        self.res2 += (1-self.f1)*leakage*dt-self.k2*self.res2*dt
        self.res1 += self.f1*leakage*dt - self.k1*self.res1*dt
        
        return {'discharge':self.discharge}

class ZanardoModel(GroundwaterZone):
    """ Two linear reservoir in series model with constant leakage from first to second
    
    Res. 1 and res. 2 are in series; both drain directly to stream.
    
    Public attributes: 
        - k1 (float): [1/T] rate at which res. 1 drains to stream
        - res1 (float): [L] current groundwater stock in res. 1
        - k2 (float): [1/T] rate at which res. 2 drains to stream
        - res2 (float): [L] current groundwater stock in res. 2
        - d (float): [L/T] rate of constant leakage from res. 1 to res. 2
    """
    def __init__(self, rew, groundwater = 0, res1 = 0, res2 = 0, k1 = 0.5, k2 = 0.2, d = 0.5):
        GroundwaterZone.__init__(self,rew)
        self.res1 = res1
        self.res2 = res2
        self.groundwater = res1 + res2
        self.k1 = k1
        self.k2 = k2
        self.d = d
        self.discharge = 0
        
    def update(self,dt,**kwargs):
        """ Update state of groundwater zone stocks.
        
        Res. 1 drains to both the stream and res. 2; res. 2 drains only to stream.
        
        Args:
            - dt (float): [T] time step
            
        Kwargs (dict) has keys:
            - leakage (float): i[L/T] incoming leakage flux from vadose zone
        
        Returns:
            - fluxes (dict): dictionary of fluxes [L/T], with keys [discharge]
        """
        leakage = kwargs['leakage']
        
        self.discharge = self.k2*self.res2 + self.k1*self.res1
        self.groundwater += (leakage - self.discharge)*dt
        self.res2 += -self.k2*self.res2*dt + self.d*dt*(self.res1>0)
        self.res1 += np.max([leakage*dt - self.k1*self.res1*dt - self.d*dt,0])

        
        return {'discharge':self.discharge}
        