from classes import *
import matplotlib.pyplot as plt

#Third Experiment


n_photons=100000
w_threshold=0.0001
m=10
delta_r=0.01
Nz=200
Nr=50
Na=30

end_depth_2=1
delta_z=end_depth_2/Nz



layer_1=Layer(id=0,init_depth=-1,end_depth=0,n=1,mua=0,mus=0,g=1,isOutLayer=True)
layer_2=Layer(id=1,init_depth=0,end_depth=end_depth_2,n=1,mua=0.1,mus=100,g=0.9,isOutLayer=False)
layer_3=Layer(id=2,init_depth=np.inf,end_depth=np.inf,n=1,mua=0.0000001,mus=0,g=1,isOutLayer=True)

layers=[layer_1,layer_2,layer_3]



rds=[]
tds=[]
fluence=[]
fluence_2=[]
for i in tqdm(range(1)):
    grid=WorldGrid(delta_r,delta_r,Nr,Nz,Na,layers)
    simulation=MC(n_photons,layers,grid)
    result,_=simulation.run_experiment(w_threshold,m)
    rds.append(result.get_Rd(n_photons,inFunctionOf="alpha"))
    tds.append(result.get_Td(n_photons,inFunctionOf="alpha"))
    fluence.append(result.get_fluence(n_photons,"z"))
    fluence_2.append(result.get_Abs(n_photons, "z"))
Rd_a=np.mean(rds,axis=0)
Td_a=np.mean(tds,axis=0)
fluence=np.mean(fluence,axis=0)
fluence_2=np.mean(fluence_2,axis=0)

#Extracted by hand from the paper
Rd_a_mcml=[0.0207,0.02,0.0198,0.0205,0.0185,0.019,0.021,0.0195,0.0208,0.021,0.021,0.021,0.0215,0.02,
           0.021,0.0205,0.0202,0.0205,0.0195,0.019,0.0185,0.017,0.016,0.014,0.012,0.01,0.008,0.0055,0.003,0.001]


plt.figure("Comparison with literature")

plt.subplot(311)
x=np.linspace(0,np.pi/2,Na)
plt.plot(x,Rd_a, label='Simulated Rd_alpha')
plt.plot(x,Rd_a_mcml, label='MCML Rd_alpha')
plt.legend()

plt.subplot(312)
plt.plot(x, Td_a,label='Simulated Td_alpha')
plt.legend()

x=np.linspace(0,1,Nz)
plt.subplot(313)
plt.plot(x, fluence,label='Simulated Fluence')
plt.plot(x, fluence_2/0.1,label='Simulated Fluence_2')
plt.yscale('log')
plt.legend()


plt.show()

print(fluence)
print(Td_a)