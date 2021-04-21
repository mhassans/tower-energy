import pandas as pd
import numpy as np
from funcs import findBestFit, sortAndNormalize

N_div=8

TCs = pd.read_csv('input/cellPositions/TCPositions_sctintillator.csv', sep=' ')

thresh ={
        37 : 1.57,
        38 : 1.58,
        39 : 1.59,
        40 : 1.60,
        41 : 1.79,
        42 :1.81,
        43 :1.83,
        44 : 1.85,
        45 : 1.86,
        46 : 1.88,
        47 : 1.90,
        48 : 1.92,
        49 : 1.93,
        50 :1.95
        }

eta_low = []
eta_high = []
eta_mid = []

for j in range(37,51):
    m = (TCs['triggercelleta']*-1)[(TCs['layer']==j) & (TCs['triggercelliphi']==1)]
    eta_low.append(pd.array(m)[-1] - (pd.array(m)[-2] - pd.array(m)[-1])/2)
    eta_high.append(pd.array(m)[0] + (pd.array(m)[0] - pd.array(m)[1])/2)
    for i in range(len(pd.array(m))):
        if pd.array(m)[i]<thresh[j]:
            eta_mid.append(0.5*(pd.array(m)[i]+pd.array(m)[i-1]))
            break

borders = pd.DataFrame({
    "layer": list(range(37,51)),
    "eta_low": eta_low,
    "eta_mid": eta_mid,
    "eta_high": eta_high,    
})

test = np.array([4, 2.3, 5.1, 0.1, 5, 7, 12, 1, 8, 3])

sort_index = np.argsort(test)#save indices before sorting

towerSortedNormed = sortAndNormalize(towerSmoothed, N_div) #returns 1D np array of float type with sum=N_div
                    
bestFit, isDegenerate = findBestFit(towerSortedNormed, N_div) #returns 1D np array of integer type with sum=N_div
if (isDegenerate):
    print('check: degenerate fit result!')
    print('')
    print(20*'-')
towerFit_array = np.zeros(len(test))
for fit_index in range(len(bestFit)):#undo sort (retrieve original index)
    towerFit_array[ sort_index[-1 - fit_index] ] = bestFit[-1 - fit_index]

