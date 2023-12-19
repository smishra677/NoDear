import msprime as msp
import numpy as np
import os
#msp.simulate(seed,samplesize,ne,mutation,recob)
#ts=msp.simulate(random_seed=42,sample_size=2,Ne=7000,length=10,mutation_rate=10e-2,recombination_rate=10e-2)



#sim={'train':1000,'test':50,'validation':500}

sim = {'train':1}

for i in sim:
	os.chdir('./'+i)



	file_rate =str(i)+'_rate.npy'
	lis=[]
	Ne= 7000
	theta =4*Ne*(10e-8)*1000
	for simNum in range(sim[i]):

		file_h = str(simNum)+'_haps.npy'
		file_P = str(simNum)+'_pos.npy'
	

		rate_array = np.random.uniform(low=10e-8, high=10e-7, size=1)

		rho = 4*Ne*rate_array*1000

		os.system('.././ms 20 1 -t '+str(theta)+' -r '+str(rho)+' 1000 >file')

		f = open('file')
		o= f.read()
		f.close()
		d = o.split('//')
		H= []
		P=[]
		for i in d[1:]:
			sp =(i.strip().split('\n'))

			if len(sp[1])>1:
				sp_pos =sp[1].split(':')[1].split()
				val = sp[3:]
				val = [list(ii) for ii in val]
				H.append(val)

				P += sp_pos


		H =np.array(H[0])
		H =H.astype('float32')
		P =np.array(P)
		P = P.astype('float32')


		np.save(file_h,H)


		lis.append(rate_array[0])
		np.save(file_h,H)
		np.save(file_P,P)

		lis.append(rate_array[0])
	arr=np.array(lis)
	np.save(file_rate,arr)
	os.chdir('..')