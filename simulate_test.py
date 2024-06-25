import msprime as msp
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
#msp.simulate(seed,samplesize,ne,mutation,recob)
#ts=msp.simulate(random_seed=42,sample_size=2,Ne=7000,length=10,mutation_rate=10e-2,recombination_rate=10e-2)


def np_har(n):
    su_=[]
    for i in range(1,n,1):
        su_.append(1/i)
    return sum(su_)

#sim={'train_50k':6000,'train_1k':6000,'train_5k':6000,'train_10k':6000,'train_100k':6000,'train_150k':6000,'train_200k':6000}
#sim={'train_10k':6000,'train_100k':6000,'train_150k':6000,'train_200k':6000}

sim={'test_50k':120}



for i in sim:
	os.chdir('./'+i)
	file_rate =str(i)+'_rate.npy'
	file_rate_recom =str(i)+'_rec_rate.npy'
	lis=[]
	lia=[]
	lis1=[]
	rate_array_1 = np.random.uniform(low=10e-13, high=10e-12, size=20)
	rate_array_2 = np.random.uniform(low=10e-12, high=10e-11, size=20)
	rate_array_3=np.random.uniform(low=10e-11, high=10e-10, size=20)
	rate_array_4= np.random.uniform(low=10e-10, high=10e-09, size=20)
	rate_array_5= np.random.uniform(low=10e-09, high=10e-08, size=20)
	rate_array_6= np.random.uniform(low=10e-08, high=10e-07, size=20)


	rate_array= np.concatenate([rate_array_1, rate_array_2, rate_array_3,rate_array_4,rate_array_5,rate_array_6], axis=0)


	#mean= np.mean(rate_array)
	#std =np.std(rate_array)

	#bins = np.linspace(10e-13, 10e-7, 100000)
	#print(bins)


	#bin_means,bins = np.histogram(rate_array, bins)
	#bins= bins[1:]

	#plt.clf()
	#bin_means= bin_means/bin_means.sum()
	#plt.bar(range(100000-1), bin_means)
	#plt.savefig(str(i)+'_density_.png')
	#print(rate_array)
	#print(len(rate_array))
	#print(max(rate_array))
	#print(list(rate_array).index(max(rate_array)))
	#print(min(rate_array))
	#print(list(rate_array).index(min(rate_array)))
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
	


		if i.split('_')[1]=='10k':
			print('10k')
			tsa=msp.sim_ancestry(samples=20,ploidy=1,population_size=70000,sequence_length=10000,recombination_rate=rate_array[simNum],model=msp.StandardCoalescent())
			ts= msp.sim_mutations(tsa,rate=1e-8,model="jc69")
		elif i.split('_')[1]=='50k':
			print('50k')
			tsa=msp.sim_ancestry(samples=20,ploidy=1, population_size=70000,sequence_length=50000,recombination_rate=rate_array[simNum],model=msp.StandardCoalescent())
			ts= msp.sim_mutations(tsa,rate=1e-8,model="jc69")
		elif i.split('_')[1]=='5k':
			print('5k')
			tsa=msp.sim_ancestry(samples=20,ploidy=1,population_size=70000,sequence_length=5000,recombination_rate=rate_array[simNum],model=msp.StandardCoalescent())
			ts= msp.sim_mutations(tsa,rate=1e-8,model="jc69")

		elif i.split('_')[1]=='1k':
			print('1k')
			tsa=msp.sim_ancestry(samples=20,ploidy=1,population_size=70000,sequence_length=1000,recombination_rate=rate_array[simNum],model=msp.StandardCoalescent())
			ts= msp.sim_mutations(tsa,rate=1e-8,model="jc69")
		elif i.split('_')[1]=='150k':
			print('150k')
			tsa=msp.sim_ancestry(samples=20,ploidy=1,population_size=70000,sequence_length=150000,recombination_rate=rate_array[simNum],model=msp.StandardCoalescent())
			ts= msp.sim_mutations(tsa,rate=1e-8,model="jc69")

		elif i.split('_')[1]=='100k':
			print('100k')
			tsa=msp.sim_ancestry(samples=20,ploidy=1,population_size=70000,sequence_length=100000,recombination_rate=rate_array[simNum],model=msp.StandardCoalescent())
			ts= msp.sim_mutations(tsa,rate=1e-8,model="jc69")
		else:
			print('200k')
			tsa=msp.sim_ancestry(samples=20,ploidy=1,population_size=70000,sequence_length=200000,recombination_rate=rate_array[simNum],model=msp.StandardCoalescent())
			ts= msp.sim_mutations(tsa,rate=1e-8,model="jc69")
		#print(ts)
		#print(rate_array[0])
		H= ts.genotype_matrix()


		'''for s in ts.sites():
			print(s.position)
			print(s)
		'''



		np.save(file_h,H)
		lia.append(ts.num_sites)
		
		P = np.array([s.position for s in ts.sites()],dtype='float32')
		print(ts.num_sites)

		rho= ((ts.num_trees)/np_har(20))

		with open(str(simNum)+".vcf", "w") as vcf_file:

			vcf_file.write("##fileformat=VCFv4.1\n")
			vcf_file.write("##source=msprime\n")

			vcf_file.write("##contig=<ID=1,length=50000>\n")  
			vcf_file.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
			vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
			for ii in range(ts.num_samples // 2):
				vcf_file.write(f"\tdiploid{ii}")
			vcf_file.write("\n")

			
			haplotypes = list(ts.haplotypes())
			for variant in ts.variants():
				
				vcf_file.write(f"1\t{int(variant.site.position)}\t.\t{variant.alleles[0]}\t{variant.alleles[1]}\t.\t.\t.\tGT")
				for ij in range(0, ts.num_samples, 2):
					
					diplotype = f"{variant.genotypes[ij]}|{variant.genotypes[ij + 1]}"
					vcf_file.write(f"\t{diplotype}")
				vcf_file.write("\n")

		#with open(str(simNum)+'.vcf', 'w') as file:
			#for var in ts.variants():
				#file.write(f"{var.site.position}\t{var.alleles}\t{var.genotypes}\n")


		lis.append(rho)
		lis1.append(rate_array[simNum])
		#lis.append(rate_array[simNum])
		np.save(file_h,H)
		#np.save(file_P,P)
	
	import matplotlib.pyplot as plt
	plt.clf()
	plt.scatter(lis,lia)
	plt.savefig(str(i)+'.png')
	#plt.show()
	arr=np.array(lis)
	np.save(file_rate,arr)

	arr1=np.array(lis1)
	np.save(file_rate_recom,arr1)
	os.chdir('..')
