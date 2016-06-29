from vadoseZone import LaioVadoseZone
from groundwaterZone import GroundwaterZone

class REW:
    def __init__(self,hillslopes=[],channel=None):
        self.hillslopes = hillslopes
        self.channel = channel      
        
class Hillslope:
    def __init__(self,temp,pet,precip,vadoseZone,groundwaterZone,**kwargs):
        for key, value in kwargs.items(): setattr(self, key, value)
        self.temp               = temp
        self.precip             = precip
        self.pet                = pet
        self.vadoseZone         = vadoseZone
        self.vadoseZoneArgs     = dict([(input,getattr(self,input,None)) for input in self.vadoseZone.fluxInputArgs])
        self.groundwaterZone    = groundwaterZone
        
def main():
    memory_address = lambda input: hex(id(input))
    
    rew = REW()
    hs = Hillslope(100,200,300,LaioVadoseZone(),GroundwaterZone(),aspect=90)
    rew.hillslopes = [hs]
    
    print(memory_address(hs.temp))
    print(memory_address(hs.precip))
    print(memory_address(hs.pet))
    print(memory_address(hs.aspect))
        
    ## option 1: pass arguments in as a dictionary
    ## pros: 
    ##  always the same, even if different models require different arguments;
    ##  dictionary values are just pointers to data -- no copying big timeseries     
    ## cons: 
    ##  have to keep track of input arguments as attr
    ##  awkward input dict construction
    # hs.vadoseZoneArgs = dict([(input,getattr(hs,input,None)) for input in hs.vadoseZone.fluxInputArgs])
    print([hs.vadoseZone.computeFluxes(**hs.vadoseZoneArgs) for hs in rew.hillslopes])
        
    ## option 2: pass arguments in individually
    ## pros: simple; also appears to not copy data
    ## cons: argument list could change for different models
    #print([hs.vadoseZone.computeFluxes(temp=hs.temp,precip=hs.precip,pet=hs.pet) for hs in rew.hillslopes])
    
    
if __name__ == '__main__': main()