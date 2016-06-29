class VadoseZone:
    def __init__(self):
        self.fluxInputArgs = []
                
    def computeFluxes(self):
        raise NameError('computeFluxes')
          
class LaioVadoseZone(VadoseZone):
    
    def __init__(self,wiltingPoint = 80):
         # make sure we get the attrs that ANY vadose zone model must have
        VadoseZone.__init__(self)  
        self.fluxInputArgs = ['tmp','pet','ppt']    
        self.wp = wiltingPoint
    
    def computeFluxes(self,**inputs):
        for k,v in inputs.items(): 
            print('%s = %s at memory address %s'%(k, v,hex(id(v))))
            
        print(inputs['temp'])   # this is how you access by default 
        ## could also do: for k,v in inputs.items(): setattr(self,k,v); handier?
        return -1