import abc
from ..helpers import Position,Direction
import numpy as np

class LightSource(metaclass=abc.ABCMeta):
    """Abstract class to give some coherence to the different light sources offered
    """
    @abc.abstractmethod
    def __init__(self,position,theta,phi):
        """Inits the light source
        Args:
            position : position object with x,y,z coordinates
            theta : theta angle in the XY plane
            phi : phi angle (azimuthal) in the XZ plane
        """
        self.position=position
        self.theta=theta
        self.phi=phi
        ux,uy,uz=self._computeDirection()        
        self.direction=Direction(ux,uy,uz)
    
    def _computeDirection(self):
        #ux=np.sin(self.phi)*np.cos(self.theta)+np.cos(self.phi)*np.cos(self.theta)-np.sin(self.phi)
        #uy=np.sin(self.phi)*np.sin(self.theta)+np.cos(self.phi)*np.sin(self.theta)+np.cos(self.phi)
        #uz=np.cos(self.phi)-np.sin(self.phi)
        ux=np.sin(self.phi)*np.cos(self.theta)
        uy=np.sin(self.phi)*np.sin(self.theta)
        uz=np.cos(self.phi)
        return round(ux, 3),round(uy, 3),round(uz, 3)
    
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
            position : position object with x,y,z coordinates
            theta : theta angle in the XY plane
            phi : phi angle (azimuthal) in the XZ plane
        """
        self.ledRadius=ledRadius
        super(LaserSource,self).__init__(**kwargs)
    def getType(self):
        return "led"