import pandas as pd
import numpy as np
from funcs import findBestFit, sortAndNormalize

N_div = 16
half_N_div = N_div//2 #Splitting over two phi slices (of 5 deg) assumed to be the same. So optmiziation performed on one only.

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



modules = {}
modulesFit = {}
towerEtaLines = []

for i in range(15, 28):
    towerEtaLines.append(i*0.0870)

for layer in borders['layer']:
    tower_u0 = np.zeros(len(towerEtaLines)-1)
    tower_u1 = np.zeros(len(towerEtaLines)-1)
    low = borders['eta_low'][borders['layer']==layer].iloc[0]
    mid = borders['eta_mid'][borders['layer']==layer].iloc[0]
    high = borders['eta_high'][borders['layer']==layer].iloc[0]
    
    for index in range(len(towerEtaLines)-1):
        if ( (towerEtaLines[index] < mid) and (towerEtaLines[index+1] > low )):
            tower_u1[index] = min(towerEtaLines[index+1], mid) - max(towerEtaLines[index], low)
        if ( (towerEtaLines[index] < high) and (towerEtaLines[index+1] > mid )):
            tower_u0[index] = min(towerEtaLines[index+1], high) - max(towerEtaLines[index], mid)
    
    modules['l'+str(layer)+'-u0'] = tower_u0
    modules['l'+str(layer)+'-u1'] = tower_u1


for module in modules:
    sort_index = np.argsort(modules[module])#save indices before sorting
    
    towerSortedNormed = sortAndNormalize(modules[module], half_N_div) #returns 1D np array of float type with sum=half_N_div
                        
    bestFit, isDegenerate = findBestFit(towerSortedNormed, half_N_div) #returns 1D np array of integer type with sum=half_N_div
    if (isDegenerate):
        print('check: degenerate fit result!')
        print(module)
        print(20*'-')
    towerFit = np.zeros(len(modules[module]))
    for fit_index in range(len(bestFit)):#undo sort (retrieve original index)
        towerFit[ sort_index[-1 - fit_index] ] = bestFit[-1 - fit_index]
    
    modulesFit[module] = towerFit

