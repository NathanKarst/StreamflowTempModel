from vadoseZone import LaioVadoseZone
from groundwaterZone import GroundwaterZone
from REW import REW

        
def main():
    memory_address = lambda input: hex(id(input))
    
    rew = REW()(100,200,300,LaioVadoseZone(),GroundwaterZone(),aspect=90)
    
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
    print(rew.vadoseZone.computeFluxes(**rew.vadoseZoneArgs))
        
    ## option 2: pass arguments in individually
    ## pros: simple; also appears to not copy data
    ## cons: argument list could change for different models
    #print([hs.vadoseZone.computeFluxes(temp=hs.temp,precip=hs.precip,pet=hs.pet) for hs in rew.hillslopes])
    
    
if __name__ == '__main__': main()