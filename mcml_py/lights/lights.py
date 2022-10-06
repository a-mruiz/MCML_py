import abc

class LightSource(metaclass=abc.ABCMeta):
    """Abstract class to give some coherence to the different light sources offered
    """
    @abc.abstractmethod
    def __init__(self,x,y,z,theta,phi):
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
        """Inits the laser source
        Args:
            x : position x
            y : position y
            z : position
        """
        super(LaserSource,self).__init__(**kwargs)
    def getType(self):
        return "laser"

class LedSource(LightSource):
    def __init__(self):
        pass
    def getType(self):
        return "led"