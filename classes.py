from dis import dis
import numpy as np
from tqdm import tqdm
from loguru import logger
import time
from datetime import timedelta
#from vprof import runner
from random import random

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
    

class Photon():
    """
    Class that represents a Photon in MC simulation
    """
    def __init__(self):
        """Will init the Photon and give it its properties
        """
        self.x=0
        self.y=0
        self.z=0
        self.ux=0
        self.uy=0
        self.uz=1
        self.w=1
        self.s=0
        self.dead=False
        self.scatterings=0
        self.layer_index=0
    

    
    def compute_movement(self,layers:list[Layer],grid)->None:
        """
        Computes the movement of the photon to the different series of layers.
        It will also score the relevant physical quantities into the world grid.
        """            
        if self.layer_index == 0 and self.isAlive():
            """
            # First movement of the photon from the source
                        
            We need to check if there is a refractive-index-mismatched interface between the "air" medium and the "tissue" medium.
            If so, specular reflectance will occur.
            
                R_sp= (n1-n2)²/(n1+n2)² (Eq: 3.2)
            
            In this way, the photon weight is decreased by R_sp when entering the new medium and changes layer
            """
            
            n1=layers[self.layer_index].n
            n2=layers[self.layer_index+1].n
            R_sp=((n1-n2)/(n1+n2))**2
            
            grid.score_specular(R_sp)
            
            self.w-=R_sp
            self.layer_index=1
        else:
            """
            # Photon moving through tissue
            
            We need to compute the distance between the current photon location and the boundary of the current layer in the direction that
            the photon propagation is happening. (Eq: 3.26)
            
                                  _
                                  | (layer.ini_depth-self.z)/self.uz  ----- if self.uz<0
                distance_boundary=| np.inf                            ----- if self.uz==0
                                  | (layer.end_depth-self.z)/self.uz  ----- if self.uz>0
                                  -           
            """
            
            if self.s==0:
                self.compute_s(layers)
            
            if self.uz<0: distance_to_boundary=(layers[self.layer_index].init_depth-self.z)/self.uz
            elif self.uz==0: distance_to_boundary=np.inf
            else: distance_to_boundary=(layers[self.layer_index].end_depth-self.z)/self.uz
            
            #print(self.s)
            #print(distance_to_boundary)
            #print()
            #
            #exit()
            #
            #print(layers[self.layer_index].mut)
            """
            Now we need to decide if self.s is greater than distance_to_boundary. In this case, the photon will hit the boundary, so we
            must move the photon to the boundary and update self.s  
            """
            if self.s<distance_to_boundary*layers[self.layer_index].mut:
                #Photon does not hit the boundary, the step fits in the current tissue layer. Update the Photon location by (s/mut) and set self.s=0 to be regenerated from random.
                self.x += self.ux * self.s/layers[self.layer_index].mut
                self.y += self.uy * self.s/layers[self.layer_index].mut
                self.z += self.uz * self.s/layers[self.layer_index].mut
                self.s = 0
                """
                Now there is a need to update the absorb and scatter properties of the tissue
                
                    -(Eq: 3.18)
                """
                
                #Compute absorption 
                self._computeAbsorption(layers,grid)
                
                #Compute scattering
                self._computeScattering(layers)
                
            else:
                #Photon hits the boundary
                self.x += self.ux * distance_to_boundary
                self.y += self.uy * distance_to_boundary
                self.z += self.uz * distance_to_boundary
                self.s -= distance_to_boundary*layers[self.layer_index].mut
                        
                """
                Now we need to compute the probability of the photon being internally reflected, depending on the angle of incidence. The value of the angle is 
                calculated with: 
                    
                    alpha_i=(cos(|self.uz|))⁻¹   (Eq: 3.28)
                
                Then, applying Snell's Law:

                    ni*sin(alphai)=nt*sin(alphat)  (Eq: 3.29)

                We can compute the critical angle like:
                
                    alpha_c=(sin(nt/ni))⁻¹

                If alpha_i > alpha_c, then the internal reflectance is 1. Otherwise, we must compute it with Fresnel's:
                
                    R(alpha_i) = 1/2[ (sin(alpha_i-alpha_t)²/sin(alpha_i+alpha_t)²)  +  (tan(alpha_i-alpha_t)²/tan(alpha_i+alpha_t)²) ]   (Eq: 3.30)
                
                """
                               
                
                #in order to find the next layer we need to know if the photon is has a Z component going upwards or downwards
                
                goingUp=False
                if self.uz < 0: #Photon is moving upwards
                    next_layer_index = self.layer_index - 1
                    goingUp = True
                else: #Photon is moving downwards
                    next_layer_index = self.layer_index + 1
                
                n_i=layers[self.layer_index].n
                n_t=layers[next_layer_index].n
                                
                alpha_i=np.arccos(np.abs(self.uz))
                alpha_c=np.arcsin(n_t/n_i)
                
                if n_i>n_t and alpha_i>alpha_c:
                    alpha_t = np.pi / 2
                    #Internal reflectance
                    R_alpha_i=1
                else:
                    alpha_t= np.arcsin(n_i * np.sin(alpha_i) / n_t)
                    if alpha_i==alpha_t and alpha_i==0: #Fresnel nromal incidence
                        R_alpha_i=np.abs((n_i-n_t)/(n_i+n_t))**2
                    elif n_i==n_t:
                        #matched boundary
                        R_alpha_i=0
                    elif -self.uz<1.0e-6:
                        R_alpha_i=1
                    else:
                        #General
                        #R_alpha_i=(1/2)*(pow((np.sin(alpha_i-alpha_t)/np.sin(alpha_i+alpha_t)),2)+pow((np.tan(alpha_i-alpha_t)/np.tan(alpha_i+alpha_t)),2))
                        R_alpha_i_1=(np.sin(alpha_i-alpha_t)/np.sin(alpha_i+alpha_t))**2
                        R_alpha_i_2=(np.tan(alpha_i-alpha_t)/np.tan(alpha_i+alpha_t))**2
                        R_alpha_i=R_alpha_i_1+R_alpha_i_2
                        R_alpha_i=R_alpha_i/2
                """
                We determine whether the photon is internally reflected or transmitted by generating a random number and comparing with R_alpha_i (Eq: 3.31)
                """
                #rnd=np.random.rand()
                rnd=random()
                if rnd <= R_alpha_i:
                    #The photon is internally reflected, so the directional cosines must be updated by reversing the z component
                    self.uz = -self.uz
                    #Then back to the first stage in the movements
                else:
                    #The photon is transmitted
                    """
                    In this case we need to know if the photon has entered another layer of tissue or exited to the ambient medium. 
                    If the photon is transmitted to the next layer of tissue, it must continue propagation with an updated direction 
                    and step size. The new directional cosines are derived from (Eq: 3.34)
                    """
                    if goingUp:
                        self.layer_index-=1
                    else:
                        self.layer_index+=1
                    if self._isPhotonExiting(layers):
                        #Photon has exited to ambient medium
                        self.dead=True
                        if goingUp:
                            #Exiting upwards, we need to score Rd (same as A_rz but with r and alpha)
                            if self.scatterings>0:
                                grid.score_Rd(self,self.w)
                            else:
                                #TODO score in the unscattered Reflectance array
                                pass
                        else:
                            #Exiting downwards, we need to score T
                            #if self.scatterings>0:
                            grid.score_Td(self,self.w)
                            
                    else:
                        #Photon is transmitted to another layer
                        self.ux=self.ux * n_i/n_t
                        self.uy=self.uy * n_i/n_t
                        self.uz=np.sign(self.uz)*np.cos(alpha_t)      
                
    def compute_s(self,layers):
        """Compute the photon step size from a uniformly distributed random variable.
            The same as MCML.
        """
        epsilon=0.000001 #avoid zero
        rnd=random()+epsilon
        self.s = -np.log(rnd)#/(layers[self.layer_index].mua/layers[self.layer_index].mus) 
          
    def _computeScattering(self,layers:list[Layer]):
        """Computes scattering of the photon.
           Scores the value in the world grid.
        """
        #First compute the cos(\theta) by applying Eq:3.22
        #rnd=np.random.rand()
        rnd=random()
        if layers[self.layer_index].g==0:
            cos_theta=2*rnd-1
        else:
            cos_theta=(1/(2*layers[self.layer_index].g))*(1+(layers[self.layer_index].g**2)-((1-layers[self.layer_index].g**2)/(1-layers[self.layer_index].g+2*layers[self.layer_index].g*rnd))**2)
        
        if cos_theta<-1:cos_theta=-1
        elif cos_theta>1:cos_theta=1
        
        #Computes the theta and azimutal angles
        #rnd=np.random.rand()
        rnd=random()
        theta=np.arccos(cos_theta)
        azimutal_phi=2*np.pi*rnd

        #Now compute the new direction of the photon (Eq:3.24) and (Eq:3.25)
        sin_theta=np.sin(theta)
        sin_phi=np.sin(azimutal_phi)
        cos_phi=np.cos(azimutal_phi)
        
        
        if np.abs(self.uz)>0.99999:
            self.ux=sin_theta*np.cos(azimutal_phi)
            self.uy=sin_theta*sin_phi
            self.uz=np.sign(self.uz)*cos_theta
        else:
            self.ux=sin_theta*(self.ux*self.uz*cos_phi-self.uy*sin_phi)/(np.sqrt(1-self.uz**2)) + self.ux*cos_theta
            self.uy=sin_theta*(self.uy*self.uz*cos_phi+self.ux*sin_phi)/(np.sqrt(1-self.uz**2)) + self.uy*cos_theta
            self.uz=-sin_theta*cos_phi*(np.sqrt(1-self.uz**2)) + self.uz*cos_theta
        
        #score an scattering event happening
        self.scatterings += 1
    
    def _computeAbsorption(self,layers:list[Layer],grid):
        """Computes absorption of the photon and reduces photon weight.
           Scores the value in the world grid.
        """
        #Computes the portion of weight to be absorbed
        absorption_weight=(layers[self.layer_index].mua/layers[self.layer_index].mut)*self.w
        #Reduce photon weight
        self.w -= absorption_weight          
        grid.score_absorption(self,absorption_weight)          
                                        
    def _isPhotonExiting(self,layers:list[Layer]):
        """Returns True when the photon is exiting (Layer 0 or Layer max)
        """
        isPhotonExiting=False
        #if self.layer_index == 0 or self.layer_index==len(layers)-1:
        #    isPhotonExiting = True
        if layers[self.layer_index].isOutLayer:
            isPhotonExiting = True
        return isPhotonExiting
                   
    def isAlive(self)->bool:
        """True if photon is alive. False otherwise
        """
        return not self.dead
      
    def terminate(self,m):
        """Plays Russian Roulette with M as assistants and 1 bullet to determine if the photon dies or not
           The probability to survive is '1' in 'm'.
        Args:
            m (int): Assistants to Russian roulette
        """
        #rnd=np.random.rand()
        rnd=random()
        if rnd<=(1/m):
            self.w=m*self.w
        else:
            self.w=0
            self.dead=True
        
        
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
    
    def score_absorption(self,photon:Photon,weight:float):
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
    
    def score_Rd(self,photon:Photon,weight:float):
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
        
    def score_Td(self,photon:Photon,weight:float):
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
            logger.error("Bad argument provided to get_Rd(). Only 'raw','r','alpha' and 'sum' are supported!")
        
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
            logger.error("Bad argument provided to get_Td(). Only 'raw','r','alpha' and 'sum' are supported!")
    
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
            scale1=2*np.pi+self.delta_r*self.delta_r*self.delta_z*n_photon
            for i_z in range(self.Nz):
                for i_r in range(self.Nr):
                    A_rz[i_r][i_z] /=(i_r+0.5)*scale1
            return A_rz
        elif inFunctionOf=="z":
            A_z=np.sum(self.A_rz,axis=0)
            scale1=1/(self.delta_z*n_photon)
            for i_z in range(self.Nz):
                A_z[i_z] *= scale1
        elif inFunctionOf=="layer":
            A_z=np.sum(self.A_rz,axis=0)
            A_l=[]
            for i_z in range(self.Nz):
                A_l[self._izToLayer(i_z)] +=A_z[i_z]
            return A_l
        elif inFunctionOf=="sum":
            A_z=np.sum(self.A_rz,axis=0)
            return np.sum(A_z,axis=0)/n_photon
        
    def _izToLayer(self,index_z):
        """Return the index to the layer according to the index
        	to the grid line system in z direction (Iz).

        	Use the center of box.
        """
        layer_index=0
        while (index_z+0.5)*self.delta_z>=self.layers[layer_index].end_depth and layer_index<len(self.layers):
            layer_index+=1
        return layer_index
        
class MC():
    def __init__(self,n_photons,layers,grid:WorldGrid):
        self.n_photons=int(n_photons)
        self.layers=layers
        self.grid=grid
        #logger.debug("Creating Monte Carlo simulation with ({}) photons and ({}) layers...",n_photons,len(layers))

        
    def run_experiment(self,w_threshold,m,debug=True):
        """Will create all the photons and run the experiment under the given conditions

        Args:
            w_threshold (float): threshold for the weight of the photon to be considered as not relevant for the simulation
            m (int): Roussian Roulette variable to determine the probability of survival of the photon
        
        Returns:
            grid (WorldGrid): Object with all the recorded information of the experiment
            time (str): Elapsed time
        """
        start_time = time.time()
        #for photon_i in tqdm(range(self.n_photons)):
        for photon_i in range(self.n_photons):    
            photon=Photon()
            while photon.isAlive():
                photon.compute_movement(self.layers,self.grid)
                #runner.run(photon.compute_movement, 'cmhp', args=(self.layers,self.grid), host='localhost', port=8000)
       
                #Accounts for the low weight photons to be rouletted or die
                if photon.isAlive():
                    if photon.w < w_threshold:
                        photon.terminate(m)    
        elapsed_time=time.time()-start_time
        formated_time=str(timedelta(seconds=elapsed_time))
        
        return self.grid,elapsed_time
        