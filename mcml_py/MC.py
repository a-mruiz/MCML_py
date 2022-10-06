from .photon import Photon
from .layer import Layer
from .grid import WorldGrid
from .lights import LaserSource,LedSource

from datetime import timedelta
import time

class MC():
    """Monte Carlo simulation
    """
    def __init__(self,n_photons:int,layers:list[Layer],grid:WorldGrid):
        """Inits the Monte Carlo simulation with the specified number of photons, layers and the given recording grid

        Args:
            n_photons (int): Number of photons to simulate
            layers ([Layer]): List of layers (from top to bottom) that represent the tissue 
            grid (WorldGrid): Grid that will record the photon interaction with the tissue
        """
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
        
        return self.grid,formated_time