import msprime as msp
import numpy as np
import os
from IPython.display import SVG, display

N0 = 70000
sim = {'non_eq': 1002}


def np_har(n):
    su_=[]
    for i in range(1,n,1):
        su_.append(1/i)
    return sum(su_)



for i in sim:
    print(i)
    os.chdir('./' + i)
    file_rate = str(i) + '_rate.npy'
    changes_rate = str(i) + '_change.npy'
    changes=[]
    lis = []
    lia = []
    rate_array_1 = np.random.uniform(low=10e-13, high=10e-12, size=167)
    rate_array_2 = np.random.uniform(low=10e-12, high=10e-11, size=167)
    rate_array_3 = np.random.uniform(low=10e-11, high=10e-10, size=167)
    rate_array_4 = np.random.uniform(low=10e-10, high=10e-09, size=167)
    rate_array_5 = np.random.uniform(low=10e-09, high=10e-08, size=167)
    rate_array_6 = np.random.uniform(low=10e-08, high=10e-07, size=167)

    rate_array = np.concatenate([rate_array_1, rate_array_2, rate_array_3, rate_array_4, rate_array_5, rate_array_6], axis=0)

    for simNum in range(sim[i]):

        N0 = 70000
        demography = msp.Demography()
        demography.add_population(name="A", initial_size=70000.0)

        print(i, simNum)

        file_h = str(simNum) + '_haps.npy'
        file_P = str(simNum) + '_pos.npy'

        poi = np.random.poisson(lam=3.0, size=1)

        for ij in poi:
            exp = np.random.exponential(scale=1.0, size=poi)

            change = [(0,N0)]
            running_time=0
            for t in range(len(exp)):
                running_time=exp[t]+running_time
                samp = N0 + np.random.randn() * 50000
                while samp < 0:
                    samp = N0 + np.random.randn() * 50000

                N0 = samp
                change.append((running_time, N0))

                demography.add_population_parameters_change(time=running_time, initial_size=N0, population='A')
        changes.append(change)
        tsa = msp.sim_ancestry(samples=20, ploidy=1, demography=demography, sequence_length=50000, recombination_rate=rate_array[simNum], model=msp.StandardCoalescent())
        ts = msp.sim_mutations(tsa, rate=1e-8, model="jc69")
        #display(SVG(ts.draw_svg(y_axis=True, size=(300, 200))))
        H = ts.genotype_matrix()

        np.save(file_h, H)
        lia.append(ts.num_sites)

        P = np.array([s.position for s in ts.sites()], dtype='float32')
        #print(ts.num_sites)

        rho = ((ts.num_trees) / np_har(20))

        print(rho)

        lis.append(rho)
        np.save(file_h, H)

    import matplotlib.pyplot as plt
    plt.clf()
    plt.scatter(lis, lia)
    plt.savefig(str(i) + '.png')
    arr = np.array(lis)
    np.save(file_rate, arr)
    print(changes)
    np.save(changes_rate,changes)
    os.chdir('..')
