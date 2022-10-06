
class Position():
    """Base class to localize something in the 3D space
    """
    def __init__(self,x,y,z):
        self.x=x
        self.y=y
        self.z=z

class Direction():
    """Base class to denote the movement of something in the 3D space
    """
    def __init__(self,ux,uy,uz):
        self.ux=ux
        self.uy=uy
        self.uz=uz
    