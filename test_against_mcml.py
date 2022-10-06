from classes import *
import matplotlib.pyplot as plt
#First Experiment

layer_1=Layer(id=0,init_depth=-1,end_depth=0,n=1,mua=0,mus=0,g=1,isOutLayer=True)
layer_2=Layer(id=1,init_depth=0,end_depth=0.02,n=1,mua=10,mus=90,g=0.75,isOutLayer=False)
layer_3=Layer(id=2,init_depth=np.inf,end_depth=np.inf,n=1,mua=0,mus=0,g=1,isOutLayer=True)

layers=[layer_1,layer_2,layer_3]

n_photons=50000
w_threshold=0.0001
m=10
delta_z=0.01
delta_r=0.01
Nz=40
Nr=50
Na=1


logger.debug("Starting experiments to check results against MCML...")

rds=[]
tds=[]
for i in tqdm(range(10)):
    grid=WorldGrid(delta_r,delta_r,Nr,Nz,Na,layers)
    simulation=MC(n_photons,layers,grid)
    result,_=simulation.run_experiment(w_threshold,m)
    rds.append(result.get_Rd(n_photons,inFunctionOf="sum"))
    tds.append(result.get_Td(n_photons,inFunctionOf="sum"))
    
Rd_simulated_1=np.mean(rds)
Td_simulated_1=np.mean(tds)
Rd_expected_1=0.09734
Td_expected_1=0.66096




#Second Experiment


layer_1=Layer(id=0,init_depth=-1,end_depth=0,n=1,mua=0,mus=0,g=1,isOutLayer=True)
layer_2=Layer(id=1,init_depth=0,end_depth=0.02,n=1.5,mua=10,mus=90,g=0,isOutLayer=False)
layer_3=Layer(id=2,init_depth=np.inf,end_depth=np.inf,n=1,mua=0,mus=0,g=1,isOutLayer=True)

layers=[layer_1,layer_2,layer_3]

n_photons=50000
w_threshold=0.0001
m=10
delta_z=0.01
delta_r=0.01
Nz=40
Nr=50
Na=1

rds=[]
for i in tqdm(range(10)):
    grid=WorldGrid(delta_r,delta_r,Nr,Nz,Na,layers)
    simulation=MC(n_photons,layers,grid)
    result,_=simulation.run_experiment(w_threshold,m)
    Rd_simulated=result.get_Rd(n_photons,inFunctionOf="sum")
    rds.append(Rd_simulated)
    
Rd_simulated_2=np.mean(rds)
Rd_expected_2=0.25907








logger.info("Total Rd by first simulation: {}",str(Rd_simulated_1))
logger.info("Expected Rd by MCML: {}",str(Rd_expected_1))
logger.info("Error-> {}",np.abs(Rd_simulated_1-Rd_expected_1))
logger.info("Total Td by first simulation: {}",str(Td_simulated_1))
logger.info("Expected Td by MCML: {}",str(Td_expected_1))
logger.info("Error-> {}",np.abs(Td_simulated_1-Td_expected_1))
logger.info("")
logger.info("Total Rd by second simulation: {}",str(Rd_simulated_2))
logger.info("Expected Rd by MCML: {}",str(Rd_expected_2))
logger.info("Error-> {}",np.abs(Rd_simulated_2-Rd_expected_2))



#Third Experiment


layer_1=Layer(id=0,init_depth=-1,end_depth=0,n=1,mua=0,mus=0,g=1,isOutLayer=True)
layer_2=Layer(id=1,init_depth=0,end_depth=0.02,n=1,mua=10,mus=90,g=0.75,isOutLayer=False)
layer_3=Layer(id=2,init_depth=np.inf,end_depth=np.inf,n=1,mua=0,mus=0,g=1,isOutLayer=True)

layers=[layer_1,layer_2,layer_3]

n_photons=50000
w_threshold=0.0001
m=10
delta_z=0.01
delta_r=0.01
Nz=40
Nr=50
Na=30

rds=[]
tds=[]
for i in tqdm(range(10)):
    grid=WorldGrid(delta_r,delta_r,Nr,Nz,Na,layers)
    simulation=MC(n_photons,layers,grid)
    result,_=simulation.run_experiment(w_threshold,m)
    rds.append(result.get_Rd(n_photons,inFunctionOf="alpha"))
    tds.append(result.get_Td(n_photons,inFunctionOf="alpha"))
    

Rd_a=np.mean(rds,axis=0)
Td_a=np.mean(tds,axis=0)

#Extracted by hand from the paper
Rd_a_mcml=[0.0207,0.02,0.0198,0.0205,0.0185,0.019,0.021,0.0195,0.0208,0.021,0.021,0.021,0.0215,0.02,
           0.021,0.0205,0.0202,0.0205,0.0195,0.019,0.0185,0.017,0.016,0.014,0.012,0.01,0.008,0.0055,0.003,0.001]


plt.figure("Comparison with literature")

plt.subplot(211)
x=np.linspace(0,np.pi/2,Na)
plt.plot(x,Rd_a, label='Simulated Rd_alpha')
plt.plot(x,Rd_a_mcml, label='MCML Rd_alpha')
plt.legend()

plt.subplot(212)
plt.plot(x, Td_a,label='Simulated Td_alpha')

plt.legend()
plt.show()










