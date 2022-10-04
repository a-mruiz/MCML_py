from classes import *

layer_1=Layer(id=0,init_depth=-1,end_depth=0,n=1,mua=0,mus=0,g=1,isOutLayer=True)
layer_2=Layer(id=1,init_depth=0,end_depth=0.02,n=1,mua=10,mus=90,g=0.75,isOutLayer=False)
layer_3=Layer(id=2,init_depth=np.inf,end_depth=np.inf,n=1,mua=0,mus=0,g=1,isOutLayer=True)

layers=[layer_1,layer_2,layer_3]

n_photons=1e6
w_threshold=0.0001
m=10
delta_z=0.01
delta_r=0.01
Nz=40
Nr=50
Na=1

grid=WorldGrid(delta_r,delta_r,Nr,Nz,Na,layers)



simulation=MC(n_photons,layers,grid)

result=simulation.run_experiment(w_threshold,m)

#runner.run(simulation.run_experiment, 'cmhp', args=(0.001, 10), host='localhost', port=8000)

#print(result.get_absorption(1e5))

#print(result.get_Rd(inFunctionOf="raw"))
#print(result.get_Rd(inFunctionOf="r"))
#print(result.get_Rd(inFunctionOf="alpha"))
print("Total Rd-> "+str(result.get_Rd(n_photons,inFunctionOf="sum")))





