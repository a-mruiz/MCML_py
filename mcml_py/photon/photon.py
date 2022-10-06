# import numpy as np
# from random import gauss, random
# from ..layer import Layer
# from ..grid import WorldGrid


# class Photon():
#     """
#     Class that represents a Photon in MC simulation
#     """
#     def __init__(self, typeOfLight="laser",ledRadius=0.5,ledHeight=10):
#         """Will init the Photon and give it its properties
#         Args:
#             typeOfLight (str): Determine the launching of the photon dependin on the different type of illumination, "laser" or "led". Defaults to "laser".
#                 This implies that the "led" case will launch photons with different initial positions (according to the led emitting surface) and
#                 that the initial directions of this photons will also change, as led bulbs generate photons in different directions.
#             ledRadius (double): This parameter only applies in the case of using a LED light source. It defines the led emitting radius (cm). Defaults to 0.5.
#         """
#         self.s=0
#         self.dead=False
#         self.scatterings=0
#         self.layer_index=0
#         self.typeOfLight=typeOfLight
#         self.ledRadius=ledRadius
#         self.ledHeight=ledHeight
#         self.w=1
        
#         if typeOfLight=="laser":
#             self.x=0
#             self.y=0
#             self.z=-1
#             self.ux=0
#             self.uy=0
#             self.uz=1
#         elif typeOfLight=="led":
#             """
#             Some considerations are needed here: 
#                 - This method will follow the dissertation in "Monte-carlo Simulation of Light Propagation considering Characteristic of Near-infrared LED and Evaluation on Tissue Phantom"
#                 - The led light will launch photons uniformly in all the led emitting surface
#                 - Photons will be launched with different angles.
#                 - #TODO
#             """
#             #First, generate random point inside the led emitting surface (https://stackoverflow.com/questions/5837572/generate-a-random-point-within-a-circle-uniformly)
#             r = ledRadius * np.sqrt(random())
#             theta = random() *2*np.pi
            
#             self.x=r*np.cos(theta)
#             self.y=r*np.sin(theta)        
        
#             #The 'z' value is going to be fixed for now to 0, no curvature in the led surface is contempled
#             self.z=0
            
#             #Now we need to randomly sample the incident angle of the photon (or the launching angle of the photon)
            
#             #phi is referencing the angle from Z on the XZ plane
#             phi=np.mean(gauss(0,90))
#             #Theta is referencing the angle from X in the XY plane
#             theta=360*random()
            
            
#             #From cylindrical to cartesian coordinates and unitary vectors->
#             self.ux=np.sin(phi)*np.cos(theta)+np.cos(phi)*np.cos(theta)-np.sin(phi)
#             self.uy=np.sin(phi)*np.sin(theta)+np.cos(phi)*np.sin(theta)+np.cos(phi)
#             self.uz=np.cos(phi)-np.sin(phi)
            
#             print(self.ux)
#             print(self.uy)
#             print(self.uz)

            
#             x=self.w*np.sin(phi)*np.cos(theta)
#             y=self.w*np.sin(phi)*np.sin(theta)
#             z=self.w*np.cos(phi)
            
                        
            
            
            

    
#     def compute_movement(self,layers:list[Layer],grid:WorldGrid)->None:
#         """
#         Computes the movement of the photon to the different series of layers.
#         It will also score the relevant physical quantities into the world grid.
#         """            
#         if self.layer_index == 0 and self.isAlive():
#             """
#             # First movement of the photon from the source
                        
#             We need to check if there is a refractive-index-mismatched interface between the "air" medium and the "tissue" medium.
#             If so, specular reflectance will occur.
            
#                 R_sp= (n1-n2)²/(n1+n2)² (Eq: 3.2)
            
#             In this way, the photon weight is decreased by R_sp when entering the new medium and changes layer
#             """
            
#             n1=layers[self.layer_index].n
#             n2=layers[self.layer_index+1].n
#             R_sp=((n1-n2)/(n1+n2))**2
            
#             grid.score_specular(R_sp)
            
#             self.w-=R_sp
#             self.layer_index=1
#             if self.typeOfLight=="laser":
#                 #IMPORTANT THIS IS ONLY VALID FOR THE USE CASE OF LASER PHOTONS
#                 #move the photon to the init of the tissue
#                 self.z=layers[self.layer_index].init_depth
#             elif self.typeOfLight=="led":
#                 #TODO compute the position of the photon when arriving to the tissue and the directional cosines
#                 pass
#         else:
#             """
#             # Photon moving through tissue
            
#             We need to compute the distance between the current photon location and the boundary of the current layer in the direction that
#             the photon propagation is happening. (Eq: 3.26)
            
#                                   _
#                                   | (layer.ini_depth-self.z)/self.uz  ----- if self.uz<0
#                 distance_boundary=| np.inf                            ----- if self.uz==0
#                                   | (layer.end_depth-self.z)/self.uz  ----- if self.uz>0
#                                   -           
#             """
            
#             if self.s==0:
#                 self.compute_s(layers)
            
#             if self.uz<0: distance_to_boundary=(layers[self.layer_index].init_depth-self.z)/self.uz
#             elif self.uz==0: distance_to_boundary=np.inf
#             else: distance_to_boundary=(layers[self.layer_index].end_depth-self.z)/self.uz
            
#             #print(self.s)
#             #print(distance_to_boundary)
#             #print()
#             #
#             #exit()
#             #
#             #print(layers[self.layer_index].mut)
#             """
#             Now we need to decide if self.s is greater than distance_to_boundary. In this case, the photon will hit the boundary, so we
#             must move the photon to the boundary and update self.s  
#             """
#             if self.s<distance_to_boundary*layers[self.layer_index].mut:
#                 #Photon does not hit the boundary, the step fits in the current tissue layer. Update the Photon location by (s/mut) and set self.s=0 to be regenerated from random.
#                 self.x += self.ux * self.s/layers[self.layer_index].mut
#                 self.y += self.uy * self.s/layers[self.layer_index].mut
#                 self.z += self.uz * self.s/layers[self.layer_index].mut
#                 self.s = 0
#                 """
#                 Now there is a need to update the absorb and scatter properties of the tissue
                
#                     -(Eq: 3.18)
#                 """
                
#                 #Compute absorption 
#                 self._computeAbsorption(layers,grid)
                
#                 #Compute scattering
#                 self._computeScattering(layers)
                
#             else:
#                 #Photon hits the boundary
#                 self.x += self.ux * distance_to_boundary
#                 self.y += self.uy * distance_to_boundary
#                 self.z += self.uz * distance_to_boundary
#                 self.s -= distance_to_boundary*layers[self.layer_index].mut
                        
#                 """
#                 Now we need to compute the probability of the photon being internally reflected, depending on the angle of incidence. The value of the angle is 
#                 calculated with: 
                    
#                     alpha_i=(cos(|self.uz|))⁻¹   (Eq: 3.28)
                
#                 Then, applying Snell's Law:

#                     ni*sin(alphai)=nt*sin(alphat)  (Eq: 3.29)

#                 We can compute the critical angle like:
                
#                     alpha_c=(sin(nt/ni))⁻¹

#                 If alpha_i > alpha_c, then the internal reflectance is 1. Otherwise, we must compute it with Fresnel's:
                
#                     R(alpha_i) = 1/2[ (sin(alpha_i-alpha_t)²/sin(alpha_i+alpha_t)²)  +  (tan(alpha_i-alpha_t)²/tan(alpha_i+alpha_t)²) ]   (Eq: 3.30)
                
#                 """
                               
                
#                 #in order to find the next layer we need to know if the photon is has a Z component going upwards or downwards
                
#                 goingUp=False
#                 if self.uz < 0: #Photon is moving upwards
#                     next_layer_index = self.layer_index - 1
#                     goingUp = True
#                 else: #Photon is moving downwards
#                     next_layer_index = self.layer_index + 1
                
#                 n_i=layers[self.layer_index].n
#                 n_t=layers[next_layer_index].n
                                
#                 alpha_i=np.arccos(np.abs(self.uz))
#                 alpha_c=np.arcsin(n_t/n_i)
                
#                 if n_i>n_t and alpha_i>alpha_c:
#                     alpha_t = np.pi / 2
#                     #Internal reflectance
#                     R_alpha_i=1
#                 else:
#                     alpha_t= np.arcsin(n_i * np.sin(alpha_i) / n_t)
#                     if alpha_i==alpha_t and alpha_i==0: #Fresnel nromal incidence
#                         R_alpha_i=np.abs((n_i-n_t)/(n_i+n_t))**2
#                     elif n_i==n_t:
#                         #matched boundary
#                         R_alpha_i=0
#                     elif -self.uz<1.0e-6:
#                         R_alpha_i=1
#                     else:
#                         #General
#                         #R_alpha_i=(1/2)*(pow((np.sin(alpha_i-alpha_t)/np.sin(alpha_i+alpha_t)),2)+pow((np.tan(alpha_i-alpha_t)/np.tan(alpha_i+alpha_t)),2))
#                         R_alpha_i_1=(np.sin(alpha_i-alpha_t)/np.sin(alpha_i+alpha_t))**2
#                         R_alpha_i_2=(np.tan(alpha_i-alpha_t)/np.tan(alpha_i+alpha_t))**2
#                         R_alpha_i=R_alpha_i_1+R_alpha_i_2
#                         R_alpha_i=R_alpha_i/2
#                 """
#                 We determine whether the photon is internally reflected or transmitted by generating a random number and comparing with R_alpha_i (Eq: 3.31)
#                 """
#                 #rnd=np.random.rand()
#                 rnd=random()
#                 if rnd <= R_alpha_i:
#                     #The photon is internally reflected, so the directional cosines must be updated by reversing the z component
#                     self.uz = -self.uz
#                     #Then back to the first stage in the movements
#                 else:
#                     #The photon is transmitted
#                     """
#                     In this case we need to know if the photon has entered another layer of tissue or exited to the ambient medium. 
#                     If the photon is transmitted to the next layer of tissue, it must continue propagation with an updated direction 
#                     and step size. The new directional cosines are derived from (Eq: 3.34)
#                     """
#                     if goingUp:
#                         self.layer_index-=1
#                     else:
#                         self.layer_index+=1
#                     if self._isPhotonExiting(layers):
#                         #Photon has exited to ambient medium
#                         self.dead=True
#                         if goingUp:
#                             #Exiting upwards, we need to score Rd (same as A_rz but with r and alpha)
#                             if self.scatterings>0:
#                                 grid.score_Rd(self,self.w)
#                             else:
#                                 #TODO score in the unscattered Reflectance array
#                                 pass
#                         else:
#                             #Exiting downwards, we need to score T
#                             if self.scatterings>0:
#                                 grid.score_Td(self,self.w)
                            
#                     else:
#                         #Photon is transmitted to another layer
#                         self.ux=self.ux * n_i/n_t
#                         self.uy=self.uy * n_i/n_t
#                         self.uz=np.sign(self.uz)*np.cos(alpha_t)      
                
#     def compute_s(self,layers):
#         """Compute the photon step size from a uniformly distributed random variable.
#             The same as MCML.
#         """
#         epsilon=0.000001 #avoid zero
#         rnd=random()+epsilon
#         self.s = -np.log(rnd)#/(layers[self.layer_index].mua/layers[self.layer_index].mus) 
          
#     def _computeScattering(self,layers:list[Layer]):
#         """Computes scattering of the photon.
#            Scores the value in the world grid.
#         """
#         #First compute the cos(\theta) by applying Eq:3.22
#         #rnd=np.random.rand()
#         rnd=random()
#         if layers[self.layer_index].g==0:
#             cos_theta=2*rnd-1
#         else:
#             cos_theta=(1/(2*layers[self.layer_index].g))*(1+(layers[self.layer_index].g**2)-((1-layers[self.layer_index].g**2)/(1-layers[self.layer_index].g+2*layers[self.layer_index].g*rnd))**2)
        
#         if cos_theta<-1:cos_theta=-1
#         elif cos_theta>1:cos_theta=1
        
#         #Computes the theta and azimutal angles
#         #rnd=np.random.rand()
#         rnd=random()
#         theta=np.arccos(cos_theta)
#         azimutal_phi=2*np.pi*rnd

#         #Now compute the new direction of the photon (Eq:3.24) and (Eq:3.25)
#         sin_theta=np.sin(theta)
#         sin_phi=np.sin(azimutal_phi)
#         cos_phi=np.cos(azimutal_phi)
        
        
#         if np.abs(self.uz)>0.99999:
#             self.ux=sin_theta*np.cos(azimutal_phi)
#             self.uy=sin_theta*sin_phi
#             self.uz=np.sign(self.uz)*cos_theta
#         else:
#             self.ux=sin_theta*(self.ux*self.uz*cos_phi-self.uy*sin_phi)/(np.sqrt(1-self.uz**2)) + self.ux*cos_theta
#             self.uy=sin_theta*(self.uy*self.uz*cos_phi+self.ux*sin_phi)/(np.sqrt(1-self.uz**2)) + self.uy*cos_theta
#             self.uz=-sin_theta*cos_phi*(np.sqrt(1-self.uz**2)) + self.uz*cos_theta
        
#         #score an scattering event happening
#         self.scatterings += 1
    
#     def _computeAbsorption(self,layers:list[Layer],grid):
#         """Computes absorption of the photon and reduces photon weight.
#            Scores the value in the world grid.
#         """
#         #Computes the portion of weight to be absorbed
#         absorption_weight=(layers[self.layer_index].mua/layers[self.layer_index].mut)*self.w
#         #Reduce photon weight
#         self.w -= absorption_weight          
#         grid.score_absorption(self,absorption_weight)          
                                        
#     def _isPhotonExiting(self,layers:list[Layer]):
#         """Returns True when the photon is exiting (Layer 0 or Layer max)
#         """
#         isPhotonExiting=False
#         #if self.layer_index == 0 or self.layer_index==len(layers)-1:
#         #    isPhotonExiting = True
#         if layers[self.layer_index].isOutLayer:
#             isPhotonExiting = True
#         return isPhotonExiting
                   
#     def isAlive(self)->bool:
#         """True if photon is alive. False otherwise
#         """
#         return not self.dead
      
#     def terminate(self,m):
#         """Plays Russian Roulette with M as assistants and 1 bullet to determine if the photon dies or not
#            The probability to survive is '1' in 'm'.
#         Args:
#             m (int): Assistants to Russian roulette
#         """
#         #rnd=np.random.rand()
#         rnd=random()
#         if rnd<=(1/m):
#             self.w=m*self.w
#         else:
#             self.w=0
#             self.dead=True
        
  
  
  
  
  
  
  
  
import numpy as np
from random import gauss, random
from ..layer import Layer
from ..grid import WorldGrid
from ..lights import LightSource

class Photon():
    """
    Class that represents a Photon in MC simulation
    """
    def __init__(self,lightSource:LightSource):
        """Will init the Photon and give it its properties
        Args:
            typeOfLight (str): Determine the launching of the photon dependin on the different type of illumination, "laser" or "led". Defaults to "laser".
                This implies that the "led" case will launch photons with different initial positions (according to the led emitting surface) and
                that the initial directions of this photons will also change, as led bulbs generate photons in different directions.
            ledRadius (double): This parameter only applies in the case of using a LED light source. It defines the led emitting radius (cm). Defaults to 0.5.
        """
        self.s=0
        self.dead=False
        self.scatterings=0
        self.layer_index=0
        self.typeOfLight=lightSource.getType()
        self.ledRadius=ledRadius
        self.ledHeight=ledHeight
        self.w=1
        
        if self.typeOfLight=="laser":
            self.x=0
            self.y=0
            self.z=-1
            self.ux=0
            self.uy=0
            self.uz=1
        elif self.typeOfLight=="led":
            """
            Some considerations are needed here: 
                - This method will follow the dissertation in "Monte-carlo Simulation of Light Propagation considering Characteristic of Near-infrared LED and Evaluation on Tissue Phantom"
                - The led light will launch photons uniformly in all the led emitting surface
                - Photons will be launched with different angles.
                - #TODO
            """
            #First, generate random point inside the led emitting surface (https://stackoverflow.com/questions/5837572/generate-a-random-point-within-a-circle-uniformly)
            r = ledRadius * np.sqrt(random())
            theta = random() *2*np.pi
            
            self.x=r*np.cos(theta)
            self.y=r*np.sin(theta)        
        
            #The 'z' value is going to be fixed for now to 0, no curvature in the led surface is contempled
            self.z=0
            
            #Now we need to randomly sample the incident angle of the photon (or the launching angle of the photon)
            
            #Phi is referencing the angle from Z on the XZ plane
            phi=np.mean(gauss(0,90))
            #Theta is referencing the angle from X in the XY plane
            theta=360*random()
            
            
            #From cylindrical to cartesian coordinates and unitary vectors->
            self.ux=np.sin(phi)*np.cos(theta)+np.cos(phi)*np.cos(theta)-np.sin(phi)
            self.uy=np.sin(phi)*np.sin(theta)+np.cos(phi)*np.sin(theta)+np.cos(phi)
            self.uz=np.cos(phi)-np.sin(phi)
            
            print(self.ux)
            print(self.uy)
            print(self.uz)

            
            x=self.w*np.sin(phi)*np.cos(theta)
            y=self.w*np.sin(phi)*np.sin(theta)
            z=self.w*np.cos(phi)
            
                        
            
            
            

    
    def compute_movement(self,layers:list[Layer],grid:WorldGrid)->None:
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
            if self.typeOfLight=="laser":
                #IMPORTANT THIS IS ONLY VALID FOR THE USE CASE OF LASER PHOTONS
                #move the photon to the init of the tissue
                self.z=layers[self.layer_index].init_depth
            elif self.typeOfLight=="led":
                #TODO compute the position of the photon when arriving to the tissue and the directional cosines
                pass
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
                            if self.scatterings>0:
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
        
  