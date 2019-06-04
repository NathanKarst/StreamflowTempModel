class GroundwaterZone:
    """ Abstract base class for all groundwater zone models. 

    Public methods:
        - update: update state of groundwater zone stocks and return dictionary of fluxes
    """

    def __init__(self): pass

    def update(self): 
        """ Update state of groundwater zone stocks.
        """
        raise NameError('update')
        
class LinearReservoir(GroundwaterZone):
    """ Simple linear reservoir model.
    
    Public attributes: 
        - storageGZ (float): [L] initial condition for groundwater stock
        - k (float): [1/T] recession parameter
        - discharge (float): [L/T] current discharge flux
    """
    def __init__(self,**kwargs):
        args = ['storageGZ','k']
        for arg in args: setattr(self, arg, kwargs[arg])
        self.overlandFlow = 0   
        self.discharge = self.k*self.storageGZ      
        
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
        
        self.discharge = self.k*self.storageGZ                # [cm/d]
        self.storageGZ += (-self.discharge + leakage)*dt      # [cm]
        
        return {'discharge':self.discharge, 'overlandFlow':self.overlandFlow}
        
class NonlinearReservoir(GroundwaterZone):
    """ Nonlinear reservoir model.
    
    Public attributes: 
        - storageGZ: [L] current groundwater stock    
        - a (float): powerlaw scale parameter
        - b (float): powerlaw exponent
        - discharge: [L/T] current discharge flux
    
    Public methods:
        - update:  (here, discharge only)
    """
    def __init__(self,**kwargs):
        args = ['storageGZ','a','b']
        for arg in args: setattr(self, arg, kwargs[arg])

        self.discharge = 0
        self.overlandFlow = 0        
        
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
        #     storageGZ = np.zeros(np.size(leakage))
        #     discharge[0] = self.a*self.groundwater**self.b
        #     groundwater[0] = (-self.discharge+leakage[0])*dt+self.groundwater
        #     for i in xrange(1,np.size(leakage)):
        #         discharge[i] = self.a*groundwater[i]**self.b
        #         groundwater[i] = (-discharge[i]+leakage[i])*dt+groundwater[i-1]

        #     self.groundwater = groundwater[-1]
        #     self.discharge = discharge[-1]
        #     return discharge, groundwater

        leakage = kwargs['leakage']

        self.discharge = self.a*self.storageGZ**self.b
        self.storageGZ += (-self.discharge + leakage)*dt

        return {'discharge':self.discharge, 'overlandFlow':self.overlandFlow}

class TwoLinearReservoir(GroundwaterZone):
    """ Two linear reservoir model.
    
    Res. 1 and res. 2 are in series. Res. 1 drains to both the stream and res. 2; res. 2 drains only to stream.
    
    Public attributes: 
        - storageGZ (float): [L] current groundwater stock     
        - k1 (float): [1/T] rate at which res. 1 drains to stream
        - res1 (float): [L] current groundwater stock in res. 1
        - k2 (float): [1/T] rate at which res. 2 drains to stream
        - res2 (float): [L] current groundwater stock in res. 2
        - k12 (float): [1/T] rate at which res. 1 drains to res. 2
        - discharge (float): [L/T] rate of discharge from groundwater zone to stream         
    """
    def __init__(self, **kwargs):
        args = ['storageGZ','k1','k2','res1','res2','k12']
        for arg in args: setattr(self, arg, kwargs[arg])
        self.storageGZ = self.res1 + self.res2
        self.overlandFlow = 0      
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
        self.storageGZ += (leakage - self.discharge)*dt
        self.res2 += self.k12*self.res1*dt - self.k2*self.res2*dt
        self.res1 += leakage*dt - self.k12*self.res1*dt - self.k1*self.res1*dt
        
        return {'discharge':self.discharge, 'overlandFlow':self.overlandFlow}

class TwoParallelLinearReservoir(GroundwaterZone):
    """ Two linear reservoir model.
    
    Res. 1 and res. 2 are in parallel; both drain directly to stream.
    
    Public attributes: 
        - storageGZ (float): [L] current groundwater stock     
        - k1 (float): [1/T] rate at which res. 1 drains to stream
        - res1 (float): [L] current groundwater stock in res. 1
        - k2 (float): [1/T] rate at which res. 2 drains to stream
        - res2 (float): [L] current groundwater stock in res. 2
        - f1 (float): [] fraction of incoming flux apportioned to res. 1
        - discharge (float): [L/T] rate of discharge from groundwater zone to stream        
    """
    def __init__(self, **kwargs):
        args = ['storageGZ','k1','k2','res1','res2','f1']
        for arg in args: setattr(self, arg, kwargs[arg])
        self.storageGZ = self.res1 + self.res2
        self.discharge = 0
        self.overlandFlow = 0
        
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
        self.storageGZ += (leakage - self.discharge)*dt
        self.res2 += (1-self.f1)*leakage*dt-self.k2*self.res2*dt
        self.res1 += self.f1*leakage*dt - self.k1*self.res1*dt
        
        return {'discharge':self.discharge, 'overlandFlow':self.overlandFlow}

class LinearToNonlinearReservoir(GroundwaterZone):
    """ Two linear reservoir model.
    
        a, b, k12, k1, res1, res2      
    """
    def __init__(self, **kwargs):
        args = ['k1','k12','res1','res2','a','b']
        for arg in args: setattr(self, arg, kwargs[arg])
        self.storageGZ = self.res1 + self.res2
        self.discharge = 0
        self.overlandFlow = 0
        
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
        discharge = self.k1*self.res1 + self.a*(self.res2)**self.b

        self.res2 += self.k12*self.res1*dt - self.a*(self.res2)**self.b*dt
        self.res1 += leakage*dt - self.k12*self.res1*dt - self.k1*self.res1*dt

        self.discharge = discharge
        self.storageGZ = self.res1 + self.res2

        return {'discharge':self.discharge, 'overlandFlow':self.overlandFlow}
        
        
        
class Melange(GroundwaterZone):
    def __init__(self, **kwargs):
        args = ['storageGZ', 'storageFAST', 'k', 'a', 'b', 'capacity']
        for arg in args: setattr(self, arg, kwargs[arg])        
        
        self.discharge = 0
        self.overlandFlow = 0
        
    def update(self, dt, **kwargs):
        leakage = kwargs['leakage']
        
        # discharge from two reservoirs
        dischargeFAST = self.k*self.storageFAST
        dischargeGZ = self.a*self.storageGZ**self.b
        self.discharge = dischargeFAST + dischargeGZ
        self.overlandFlow = dischargeFAST
        
        # update storages
        self.storageGZ += (leakage - dischargeGZ)*dt
        self.storageFAST -= dischargeFAST*dt
        
        # if GZ is > capacity, add extra to fast bucket
        if self.storageGZ > self.capacity:
            self.storageFAST = (self.storageGZ - self.capacity)
            self.storageGZ = self.capacity
            
        return {'discharge':self.discharge, 'overlandFlow':self.overlandFlow}
         
        