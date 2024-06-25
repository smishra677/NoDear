import msprime as msp
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import SVG, display

N0 = 70000
sim = {'non_eq_train': 6000}

def np_har(n):
    su_ = []
    for i in range(1, n, 1):
        su_.append(1/i)
    return sum(su_)

for i in sim:
    print(i)
    os.chdir('./' + i)
    file_rate = str(i) + '_rate.npy'
    file_rate_recom = str(i) + '_rec_rate.npy'
    file_rate_R = str(i) + '_R.npy'
    lis = []
    lis1 = []
    lia = []
    lis_R=[]

    rate_array_1 = np.random.uniform(low=10e-13, high=10e-12, size=1000)
    rate_array_2 = np.random.uniform(low=10e-12, high=10e-11, size=1000)
    rate_array_3 = np.random.uniform(low=10e-11, high=10e-10, size=1000)
    rate_array_4 = np.random.uniform(low=10e-10, high=10e-09, size=1000)
    rate_array_5 = np.random.uniform(low=10e-09, high=10e-08, size=1000)
    rate_array_6 = np.random.uniform(low=10e-08, high=10e-07, size=1000)
    rate_array = np.concatenate([rate_array_1,rate_array_2,rate_array_3, rate_array_4, rate_array_5,rate_array_6], axis=0)
    
    #rate_array = np.concatenate([rate_array_1], axis=0)
    for simNum in range(sim[i]):
        changes = []
        changes_rate = str(simNum) + '.csv'
        file_h = str(simNum) + '_haps.npy'
        file_P = str(simNum) + '_pos.npy'
        if True:
            N0 = 70000
            demography = msp.Demography()
            demography.add_population(name="A", initial_size=70000.0)

            print(i, simNum)



            poi = np.random.poisson(lam=3.0, size=1)

            for ij in poi:
                exp = np.random.exponential(scale=1.0, size=poi)

                change = [('ABC', 0, N0, 'path', 0)]
                running_time = 0
                for t in range(len(exp)):
                    running_time = exp[t] + running_time
                    samp = N0 + np.random.randn() * 50000
                    while samp < 0:
                        samp = N0 + np.random.randn() * 50000

                    N0 = samp
                    change.append(('ABC', running_time, N0, 'path', 0))

                    demography.add_population_parameters_change(time=running_time, initial_size=N0, population='A')
            changes.append(change)

        df = pd.DataFrame(change, columns=['label', 'x', 'y', 'plot_type', 'plot_num'])
        df.to_csv(changes_rate, index=False)

        tsa = msp.sim_ancestry(samples=20, ploidy=1, demography=demography, sequence_length=50000, recombination_rate=rate_array[simNum], model=msp.StandardCoalescent())
        ts = msp.sim_mutations(tsa, rate=1e-8, model="jc69")
        H = ts.genotype_matrix()
        lia.append(ts.num_sites)
        P = np.array([s.position for s in ts.sites()], dtype='float32')

        rho = ((ts.num_trees) / np_har(20))

        with open(str(simNum) + ".vcf", "w") as vcf_file:
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


        lis_R.append(ts.num_trees)
        
        lis1.append(rate_array[simNum])
        lis.append(rho)
        np.save(file_h, H)

    plt.clf()
    plt.scatter(lis1, lia)
    plt.savefig(str(i) + '.png')
    arr = np.array(lis)
    np.save(file_rate, arr)

    arr1 = np.array(lis1)
    np.save(file_rate_recom, arr1)


    arr2 = np.array(lis_R)
    np.save(file_rate_R, arr2)
    os.chdir('..')
