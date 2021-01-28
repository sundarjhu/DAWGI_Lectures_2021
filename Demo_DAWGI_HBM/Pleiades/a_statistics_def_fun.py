from scipy.optimize import minimize
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt

#############################################
################# HISTOGRAM #################
#############################################

def a_stat(a):

    Na=len(a)

    perc = []
    for i in range(Na):

        p = np.percentile(a[i],[0,50,100])

        perc.append([p[0], p[1], p[2]])

    perc = np.array(perc)
#    Nzeros=len(perc[0])
#    perc[0] = np.zeros(Nzeros)
#    perc[-1] = np.zeros(Nzeros)

    return perc



