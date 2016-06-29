class REW:
    def __init__(self,temp,pet,precip,vadoseZone,groundwaterZone,**kwargs):
        for key, value in kwargs.items(): setattr(self, key, value)
        self.temp               = temp
        self.precip             = precip
        self.pet                = pet
        self.vadoseZone         = vadoseZone
        self.vadoseZoneArgs     = self.asssembleInputs(self,self.vadoseZone.fluxInputArgs)
        self.groundwaterZone    = groundwaterZone
        
    def assembleInputs(self,obj,inputs): return dict([(input,getattr(obj,input,None)) for input in inputs])