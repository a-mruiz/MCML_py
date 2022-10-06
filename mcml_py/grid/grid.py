import numpy as np
from ..layer import Layer
#from ..photon import Photon

class WorldGrid():
    """
    In order to score physical quantities in the MC simulation, we need to set up grid systems.
    
        - For scoring absorption, a two-dimensional homogeneous grid in the 'r' and 'z' directions, with grid separations of 'delta_r' and 'delta_z'. The total number of grid elements
            in the 'r' and 'z' direction are 'Nr' and 'Nz'.
        - For scoring diffuse reflectance (Rd) and transmittance, a two-dimensional homogeneous grid sustem is set up in the 'r' and 'alpha' directions. This grid system shares the 'r'
            direction with the grid system for photon absorption, therefore we need to set up one extra one-dimensional grid of elements in Nalpha. 
            The grid separation being delta_alpha=pi/(2*Nalpha), since alpha is in range [0,pi/2].
            Alpha is the angle between the photon exiting direction and the normal (-z axis for reflectance
            and z axis for transmittance).
    
    The simulation propagates photons in three dimensions and records photon deposition A_rz in the neighborhood of (r,z).
    """
    def __init__(self,delta_r,delta_z,Nr,Nz,Na,layers:list[Layer]):
        self.Nr=Nr
        self.Nz=Nz
        self.Na=Na
        self.delta_r=delta_r
        self.delta_z=delta_z
        self.delta_a=np.pi/(2*Na)
        self.layers=layers
        
        
        #Initialize the absorption 2D array, which is in function of 'r' and 'z'
        self.A_rz=np.zeros([self.Nr,self.Nz])
        self.Rd_ra=np.zeros([self.Nr,self.Na])
        self.Td_ra=np.zeros([self.Nr,self.Na])
        #Total Specular Reflectance
        self.total_Rsp=0
        #Total Transmittance
        self.total_T=0
        #Total Diffuse Reflectance
        self.total_Rd=0
    
    def score_absorption(self,photon,weight:float):
        """Scores the absorption of the photon weight in the tissue.
           The weight drop is dw = w*mua/(mua+mus).
           
           The dropped weight is assigned to the absorption array
           A(r,z).

        Args:
            photon (Photon): Photon that has interacted with the tissue
            weight (float): Weight lost due to absorbance
        """
        #Obtain the coords in cylindrical system
        r_coord = np.sqrt(photon.x**2 + photon.y**2)
        z_coord = photon.z
        
        #Round to Int number to compute an homogeneous grid
        r_coord_round = int(r_coord/self.delta_r)
        z_coord_round = int(z_coord/self.delta_z)
        
        #If the r_coord is outside the world grid or it is the first interaction append it to the last element
        if r_coord_round>=self.Nr-1 or photon.scatterings==0:
            r_coord_round=-1

        #If the z_coord is outside the world grid append it to the last element
        if z_coord_round>=self.Nz-1:
            z_coord_round=-1
            
        #score the absorption
        self.A_rz[r_coord_round,z_coord_round] += weight
    
    def score_Rd(self,photon,weight:float):
        """Scores the reflectance of the photon when exiting the tissue back to the light source layer.
                     
           The exited photon weight is assigned to the reflectance array
           Rd(r,alpha). Where alpha is the angle between the photon exiting
           direction and the (-z) axis (when Rd).

        Args:
            photon (Photon): Photon that has interacted with the tissue
            weight (float): Weight that has exited the tissue
        """
        #Obtain the coords in cylindrical system
        r_coord = np.sqrt(photon.x**2 + photon.y**2)
        alpha_coord = np.arccos(-photon.uz)
        
        #Round to Int number to compute an homogeneous grid
        r_coord_round = int(r_coord/self.delta_r)
        alpha_coord_round = int(alpha_coord/self.delta_a)
        
        #If the r_coord is outside the world grid or it is the first interaction append it to the last element
        if r_coord_round>=self.Nr-1 or photon.scatterings==0:
            r_coord_round=-1

        #If the a_coord is outside the world grid append it to the last element
        if alpha_coord_round>=self.Na-1:
            alpha_coord_round=-1
            
        #score the Rd
        self.Rd_ra[r_coord_round,alpha_coord_round] += weight
        
    def score_Td(self,photon,weight:float):
        """Scores the Transmittance of the photon when escaping the tissue in the last layer.
                     
           The exited photon weight is assigned to the transmittance array
           Td(r,alpha). Where alpha is the angle between the photon exiting
           direction and the (z) axis (when Td).

        Args:
            photon (Photon): Photon that has interacted with the tissue
            weight (float): Weight that has exited the tissue
        """
        #Obtain the coords in cylindrical system
        r_coord = np.sqrt(photon.x**2 + photon.y**2)
        alpha_coord = np.arccos(photon.uz)
        
        #Round to Int number to compute an homogeneous grid
        r_coord_round = int(r_coord/self.delta_r)
        alpha_coord_round = int(alpha_coord/self.delta_a)
        
        #If the r_coord is outside the world grid or it is the first interaction append it to the last element
        if r_coord_round>=self.Nr-1 or photon.scatterings==0:
            r_coord_round=-1

        #If the a_coord is outside the world grid append it to the last element
        if alpha_coord_round>=self.Na-1:
            alpha_coord_round=-1
            
        #score the Rd
        self.Td_ra[r_coord_round,alpha_coord_round] += weight
        
    def score_specular(self,specular_r):
        self.total_Rsp+=specular_r
    
    def score_transmittance(self,transmittance):
        #TODO score transmittance when photon scapes to the last layer
        #transmittance + portion of w that scaped
        pass

    def score_reflectance(self,reflectance):
        #TODO score reflectance when photon scapes back to the source of light
        #reflectance + portion of w that scaped
        pass
    
    def get_Rd(self,n_photon,inFunctionOf="raw"):
        """Get Rd in function of different elements
           Also scales the result to obtain photon probability per area

        Args:
            n_photon (int): Number of photons in the simulation
            inFunctionOf (str, optional): "raw","r","alpha","sum". Defaults to "raw".

        Returns:
            Rd
        """        
        
        if inFunctionOf=="raw":
            scale1 = 4.0*np.pi*np.pi*self.delta_r*np.sin(self.delta_a/2)*self.delta_r*n_photon
            Rd_ra=self.Rd_ra
            for i_r in range(self.Nr):
                for i_a in range(self.Na):
                    scale2=1/((i_r+0.5)*np.sin(2*(i_a+0.5)*self.delta_a)*scale1)
                    Rd_ra[i_r][i_a] *= scale2
            return Rd_ra
        
        elif inFunctionOf=="r":
            scale1=2*np.pi*self.delta_r*self.delta_r*n_photon
            Rd_r=np.sum(self.Rd_ra,axis=1)
            for i_r in range(self.Nr):
                scale2=1/((i_r+0.5)*scale1)
                Rd_r[i_r] *= scale2
            return Rd_r
        
        elif inFunctionOf=="alpha":
            scale1=2*np.pi*self.delta_a*n_photon
            Rd_a=np.sum(self.Rd_ra,axis=0)
            for i_a in range(self.Na):
                scale2=1/(np.sin((i_a+0.5)*self.delta_a)*scale1)
                Rd_a[i_a] *= scale2
            return Rd_a
        
        elif inFunctionOf=="sum":
            Rd_r=np.sum(self.Rd_ra,axis=1)
            return np.sum(Rd_r,axis=0)/n_photon
        
        else:
            print("Bad argument provided to get_Rd(). Only 'raw','r','alpha' and 'sum' are supported!")
        
    def get_Td(self,n_photon,inFunctionOf="raw"):
        """Get Td in function of different elements
           Also scales the result to obtain photon probability per area
        Args:
            n_photon (int): Number of photons in the simulation
            inFunctionOf (str, optional): "raw","r","alpha","sum". Defaults to "raw".

        Returns:
            Td
        """
        if inFunctionOf=="raw":
            scale1 = 4.0*np.pi*np.pi*self.delta_r*np.sin(self.delta_a/2)*self.delta_r*n_photon
            Td_ra=self.Td_ra
            for i_r in range(self.Nr):
                for i_a in range(self.Na):
                    scale2=1/((i_r+0.5)*np.sin(2*(i_a+0.5)*self.delta_a)*scale1)
                    Td_ra[i_r][i_a] *= scale2
            return Td_ra
        
        elif inFunctionOf=="r":
            scale1=2*np.pi*self.delta_r*self.delta_r*n_photon
            Td_r=np.sum(self.Td_ra,axis=1)
            for i_r in range(self.Nr):
                scale2=1/((i_r+0.5)*scale1)
                Td_r[i_r] *= scale2
            return Td_r
        
        elif inFunctionOf=="alpha":
            scale1=2*np.pi*self.delta_a*n_photon
            Td_a=np.sum(self.Td_ra,axis=0)
            for i_a in range(self.Na):
                scale2=1/(np.sin((i_a+0.5)*self.delta_a)*scale1)
                Td_a[i_a] *= scale2
            return Td_a
        
        elif inFunctionOf=="sum":
            Td_r=np.sum(self.Td_ra,axis=1)
            return np.sum(Td_r,axis=0)/n_photon
        else:
            print("Bad argument provided to get_Td(). Only 'raw','r','alpha' and 'sum' are supported!")
    
    def get_Abs(self,n_photon,inFunctionOf="raw"):
        """Returns the absorbance values in function of the given property.
        

        Args:
            n_photon (int): Number of photons in the simulation
            inFunctionOf (str, optional): "raw","z","layer","sum". Defaults to "raw".

        Returns:
            Absorbance
        """
        if inFunctionOf=="raw":
            A_rz=self.A_rz
            scale1=2*np.pi*self.delta_r*self.delta_r*self.delta_z*n_photon
            for i_z in range(self.Nz):
                for i_r in range(self.Nr):
                    A_rz[i_r][i_z] /=(i_r+0.5)*scale1
            return A_rz
        elif inFunctionOf=="z":
            A_z=np.sum(self.A_rz,axis=0)
            scale1=1/(self.delta_z*n_photon)
            for i_z in range(self.Nz):
                A_z[i_z] *= scale1
            return A_z
        elif inFunctionOf=="layer":
            A_z=np.sum(self.A_rz,axis=0)
            A_l=[]
            for i_z in range(self.Nz):
                A_l[self._izToLayer(i_z)] +=A_z[i_z]/n_photon
            return A_l
        elif inFunctionOf=="sum":
            A_z=np.sum(self.A_rz,axis=0)
            return np.sum(A_z,axis=0)/n_photon
        
    def get_fluence(self,n_photon,inFunctionOf="raw"):
        """Returns the photon fluence(cm⁻²)

        Args:
            n_photon (_type_): _description_
            inFunctionOf (str, optional): _description_. Defaults to "raw".

        Returns:
            _type_: _description_
        """
        if inFunctionOf=="raw":
            A_rz=self.get_Abs(n_photon,"raw")
            fluence_rz=A_rz
            for i_z in range(self.Nz):
                layer_index=self._izToLayer(i_z)
                mua_local=self.layers[layer_index].mua
                for i_r in range(self.Nr):
                    fluence_rz[i_r][i_z]=A_rz[i_r][i_z]/mua_local
            return fluence_rz
        elif inFunctionOf=="z":
            A_z=self.get_Abs(n_photon,"z")
            fluence_z=A_z
            for i_z in range(self.Nz):
                layer_index=self._izToLayer(i_z)
                #print(layer_index)
                mua_local=self.layers[layer_index].mua
                #print(mua_local)
                fluence_z[i_z]=A_z[i_z]/mua_local
            return fluence_z
        
    def _izToLayer(self,index_z):
        """Return the index to the layer according to the index
        	to the grid line system in z direction (Iz).

        	Use the center of box.
        """
        layer_index=0
        while (index_z+0.5)*self.delta_z>=self.layers[layer_index].end_depth and layer_index<len(self.layers):
            layer_index+=1
        return layer_index
     
    