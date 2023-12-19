import msprime as msp
import numpy as np
import os
#msp.simulate(seed,samplesize,ne,mutation,recob)
#ts=msp.simulate(random_seed=42,sample_size=2,Ne=7000,length=10,mutation_rate=10e-2,recombination_rate=10e-2)



sim={'train':1}



for i in sim:
	os.chdir('./'+i)
	file_rate =str(i)+'_rate.npy'
	lis=[]
	lia=[]
	rate_array = np.random.uniform(low=10e-13, high=10e-7, size=sim[i])
	print(rate_array)
	print(len(rate_array))
	print(max(rate_array))
	print(list(rate_array).index(max(rate_array)))
	print(min(rate_array))
	print(list(rate_array).index(min(rate_array)))
	for simNum in range(sim[i]):
		print(i, simNum)
		#rate_array = np.random.uniform(low=10e-8, high=10e-7, size=1)
		#rate_array = [10e-6]
		#rate_array = [10e-7]
		#rate_array = [10e-8]
		'''	
		if i=='train':
			rate_array=[10e-6]
		if i=='test':
			rate_array = [10e-7]
		if i=='validation':
			rate_array = [10e-8]
		

		print(len(rate_array))
		print(len(list(range(haps_length))))
		rate_map = msp.RateMap(
		    position=list(range(haps_length)),
		    rate=rate_array
		)
		

		rate_map = msp.RateMap(
		    position=[0, 10,20, 30],rate=[0, 0.1, 0.6])
		
		'''

		file_h = str(simNum)+'_haps.npy'
		file_P = str(simNum)+'_pos.npy'
	



		ts=msp.simulate(sample_size=20,Ne=70000,length=10000,mutation_rate=10e-5,recombination_rate=rate_array[simNum])
		#print(ts)
		#print(rate_array[0])
		
		H= ts.genotype_matrix()
	
		'''for s in ts.sites():
			print(s.position)
			print(s)
		'''


		np.save(file_h,H)
		lia.append(ts.num_sites)
		print(ts.num_sites)
		P = np.array([s.position for s in ts.sites()],dtype='float32')

		lis.append(rate_array[simNum])
		np.save(file_h,H)
		#np.save(file_P,P)
	
	import matplotlib.pyplot as plt
	plt.scatter(lis,lia)
	plt.savefig(str(i)+'.png')
	#plt.show()
	arr=np.array(lis)
	np.save(file_rate,arr)
	os.chdir('..')
