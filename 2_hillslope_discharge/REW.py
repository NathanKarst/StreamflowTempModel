class REW:
    def __init__(self,temperature,potentialEvapotranspiration,precipitation,vadoseZone,groundwaterZone,**kwargs):
        for key, value in kwargs.items(): setattr(self, key, value)
        self.tmp        = temperature
        self.pet        = potentialEvapotranspiration        
        self.ppt        = precipitation
        self.vz         = vadoseZone(self)
        self.gz         = groundwaterZone(self)
        
        