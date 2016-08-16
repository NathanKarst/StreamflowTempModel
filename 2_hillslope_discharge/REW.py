class REW:
    ''' Representative watershed element class.
    
    Arguments:
        - vadoseZone: :class:`vadoseZone` class object, not instance
        - grounwaterZone: :class:`groundwaterZone` class object, not instance
    
    Public attributes:
        - vz: :class:`vadoseZone` instance
        - gz: :class:`groundwaterZone` instance
    
    '''
    def __init__(self,vadoseZone,groundwaterZone,**kwargs):

        for key, value in kwargs.items(): setattr(self, key, value)

        ## these inits must go last, b/c for instance vz might need to reference 
        ## some of the forcing data in the containing REW instance
        ## NOTE: there will still be a problem here if vz needs to refer to gz and vice versa
        self.vz         = vadoseZone(self)
        self.gz         = groundwaterZone(self)
