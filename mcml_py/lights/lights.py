import abc

class LightSource(metaclass=abc.ABCMeta):
    """Abstract class to give some coherence to the different light sources offered
    """
    @abc.abstractmethod
    def __init__(self,x,y,z,theta,phi):
        """Inits the light source
        Args:
            x : position x
            y : position y
            z : position z
            theta : theta angle in the XY plane
            phi : phi angle (azimuthal) in the XZ plane
        """
        self.x=x
        self.y=y
        self.z=z
        self.theta=theta
        self.phi=phi        
    
    @abc.abstractmethod
    def getType(self):
        pass
    
    
    
class LaserSource(LightSource):
    def __init__(self,**kwargs):
        super(LaserSource,self).__init__(**kwargs)
    def getType(self):
        return "laser"

class LedSource(LightSource):
    def __init__(self,ledRadius,**kwargs):
        """Inits the light source
        Args:
            ledRadius : Radius (cm) of the led bulb
            x : position x
            y : position y
            z : position z
            theta : theta angle in the XY plane
            phi : phi angle (azimuthal) in the XZ plane
        """
        self.ledRadius=ledRadius
        super(LaserSource,self).__init__(**kwargs)
    def getType(self):
        return "led"