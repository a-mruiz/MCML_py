
class Layer():
    """
    Class that represents a Layer in MC simulation
    """
    def __init__(self,id,init_depth, end_depth,n,mua,mus,g,isOutLayer:bool):
        """Inits the Layer with the given properties

        Args:
            id (int): Descriptive ID of the layer (FIRST LAYER MUST HAVE ID=0, SECOND ID=1...)
            init_depth (float): Position of the z axis where the layer starts (cm)
            end_depth (float): Position of the z axis where the layer ends (cm)
            n (float): Refractive index of the layer
            mua (float): Attenuation coefficient of the layer (1/cm)
            mus (float): Scattering coefficient of the layer (1/cm)
            mut (float): Transport coefficient of the layer (mua+mus) (1/cm)
            g (float): Anisotropy factor for scattering
            isOutLayer (bool): True for outside layers (from light to tissue and when light exits tissue (transmitance))
        """
        self.id=id
        self.init_depth=init_depth
        self.end_depth=end_depth
        self.n=n
        self.mua=mua
        self.mus=mus
        self.mut=self.mua+self.mus 
        self.g=g
        self.isOutLayer=isOutLayer
        
    def __dir__(self):
        return super().__dir__() + [str(k) for k in self.keys()]
    