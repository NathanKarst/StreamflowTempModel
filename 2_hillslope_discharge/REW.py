class REW:
    ''' Representative watershed element class.
    
    Arguments:
        - vadoseZone: :class:`vadoseZone` instance
        - grounwaterZone: :class:`groundwaterZone` instance
    
    Public attributes:
        - vz: :class:`vadoseZone` instance
        - gz: :class:`groundwaterZone` instance
    
    '''
    def __init__(self,vadoseZone,groundwaterZone,**kwargs):

        for key, value in kwargs.items(): setattr(self, key, value)
        self.vz = vadoseZone
        self.gz = groundwaterZone
