
import os
import random
import numpy as np
import pickle
import absl.logging
absl.logging.set_verbosity(absl.logging.ERROR)
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter

def normalizeTargets(i,file):
    rate_file = str(i)+'_rate.npy'
    infoFilename= os.path.join(file,rate_file)
    nTargets_1 = pickle.load(open(infoFilename))

    norm = 'zscore'
    if(norm == 'zscore'):
        tar_mean = np.mean(nTargets_1,axis=0)
        tar_sd = np.std(nTargets_1,axis=0)
        nTargets =nTargets_1- tar_mean
        nTargets = np.divide(nTargets,tar_sd,out=np.zeros_like(nTargets),where=tar_sd!=0)

    return nTargets[i]

def z_normalizeTargets(error,file):
    infoFilename= os.path.join('/N/u/samishr/Carbonate/Desktop/Recombination_landscape/ReLERNN-master/examples/example_output_pool/'+file+'/',"info.p")
    nTargets_1 = pickle.load(open(infoFilename,"rb"))['rho']

    tar_mean = np.mean(nTargets_1,axis=0)
    tar_sd = np.std(nTargets_1,axis=0)
    erro = np.dot(error,tar_sd)+tar_mean


    return erro


def SNP_WINDOW(w,r,i,file):
 
    respectiveNormalizedTargets = [normalizeTargets(i,file)]
    targets = np.array(respectiveNormalizedTargets)
    lis=[]

    for i in w:
        co = dict(Counter(i))
        if len(co.keys())==2:
            lis.append(co[min(co,key=co.get)]/len(i))
    bins = np.linspace(0, 1, 20)

    bin_means,bins = np.histogram(lis, bins, normed=True, density=True)
    bins= bins[1:]
    bin_means= bin_means/bin_means.sum()
    

def __data_generation1(batchTreeIndices,val):

        haps = []
        pos = []
        for treeIndex in batchTreeIndices:
            Hfilepath = os.path.join('/N/u/samishr/Carbonate/Desktop/Recombination_landscape/ReLERNN-master/examples/example_output_pool/'+val+'/val/',str(treeIndex) + "_haps.npy")
            Pfilepath = os.path.join('/N/u/samishr/Carbonate/Desktop/Recombination_landscape/ReLERNN-master/examples/example_output_pool/'+val+'/val/',str(treeIndex) + "_pos.npy")
            H = np.load(Hfilepath)
            P = np.load(Pfilepath)
            haps.append(H)
            pos.append(P)

       
        return haps,pos


            
def __getitem__1(idx,val):
    X,x1= __data_generation1([idx],val)
    return X,x1