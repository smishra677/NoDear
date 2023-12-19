import msprime as msp
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter




fp =pd.read_csv('/N/u/samishr/Quartz/Desktop/LDhelmet-main/recombination_maps/dmel_rec_maps/RAL_chr2L.txt', sep='\t', lineterminator='\r')
fp.columns
#print(fp.columns)
#print(msp.RateMap.read_hapmap(fp))