import msprime as msp
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

def np_har(n):
    su_=[]
    for i in range(1,n,1):
        su_.append(1/i)
    return sum(su_)

sim = {'train_50k': 120}
os.chdir('./test_sim')

for i in sim:
    print(i)
    lis = []
    lis1 = []
    lia = []
    file_rate = 'test_rate.npy'

    rate_array_1 = np.random.uniform(low=10e-13, high=10e-12, size=20)
    rate_array_2 = np.random.uniform(low=10e-12, high=10e-11, size=20)
    rate_array_3 = np.random.uniform(low=10e-11, high=10e-10, size=20)
    rate_array_4 = np.random.uniform(low=10e-10, high=10e-09, size=20)
    rate_array_5 = np.random.uniform(low=10e-09, high=10e-08, size=20)
    rate_array_6 = np.random.uniform(low=10e-08, high=10e-07, size=20)
    rate_array= np.concatenate([rate_array_1, rate_array_2, rate_array_3,rate_array_4,rate_array_5,rate_array_6], axis=0)

    for simNum in range(sim[i]):
        print(i, simNum)

        file_h = str(simNum) + '_haps.npy'
        file_P = str(simNum) + '_pos.npy'

        tsa = msp.sim_ancestry(samples=20, ploidy=1, population_size=70000, sequence_length=50000, recombination_rate=rate_array[simNum], model=msp.StandardCoalescent())
        ts = msp.sim_mutations(tsa, rate=1e-8, model="jc69")

        H = ts.genotype_matrix()

        haps = []
        for i in ts.haplotypes():
            haps.append(i)

        sequence_IDs = []
        for i in range(len(haps)):
            sequence_IDs.append(f'sample_{ts.samples()[i]}_pop_{ts.node(i).population}')

        with open(str(simNum) + '.fasta', 'w') as f:
            for i in range(len(haps)):
                f.write(f'>{sequence_IDs[i]}\n{haps[i]}\n')
                #asads=1

        np.save(file_h, H)
        lia.append(ts.num_sites)

        P = np.array([s.position for s in ts.sites()], dtype='float32')
        print(ts.num_sites)
        rho= ((ts.num_trees)/np_har(20))
        lis.append(rho)
        lis1.append(rate_array[simNum])
        np.save(file_h, H)

        os.system('./LDhelmet-main/ldhelmet find_confs --num_threads 24 -w 50 -o ./output_' + str(simNum) + '.conf ./' + str(simNum) + '.fasta')
        os.system('time ./LDhelmet-main/ldhelmet table_gen --num_threads 24 -t 5.6e-8 -r 0.0 0.1 10.0 1.0 100.0 -c ./output_' + str(simNum) + '.conf -o ./output_' + str(simNum) + '.lk')
        os.system('./LDhelmet-main/ldhelmet pade --num_threads 24 -t 5.6e-8 -x 12 --defect_threshold 40 -c ./output_' + str(simNum) + '.conf -o ./output_' + str(simNum) + '.pade')
        os.system('time ./LDhelmet-main/ldhelmet rjmcmc --num_threads 24 -l ./output_' + str(simNum) + '.lk -p ./output_' + str(simNum) + '.pade -s ./' + str(simNum) + '.fasta -b 50 -w 50  --burn_in 10000 -n 100000 -o ./output_' + str(simNum) + '.post')
        os.system('time ./LDhelmet-main/ldhelmet post_to_text -m -p 0.025 -p 0.50 -p 0.975 -o ./output_' + str(simNum) + '.txt ./output_' + str(simNum) + '.post')
        os.system('./LDhelmet-main/ldhelmet max_lk --num_threads 24 -l ./output_' + str(simNum) + '.lk -p ./output_' + str(simNum) + '.pade -s ./' + str(simNum) + '.fasta')

    plt.clf()
    plt.scatter(lis, lia)
    plt.savefig(str(simNum) + '.png')
    arr = np.array(lis)
    np.save(file_rate, arr)
    arr1=np.array(lis1)
    np.save(file_rate_recom,arr1)
    os.chdir('..')
