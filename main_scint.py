import pandas as pd
import numpy as np
from funcs import findBestFit, sortAndNormalize, weight, getModulesPerBundle, find_v

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



luSlices = {} #layer and u slices
luSlicesFit = {}
towerEtaLines = []

for i in range(15, 28):
    towerEtaLines.append(i*0.0870)

for layer in borders['layer']:
    slice_u0 = np.zeros(len(towerEtaLines)-1)
    slice_u1 = np.zeros(len(towerEtaLines)-1)
    low = borders['eta_low'][borders['layer']==layer].iloc[0]
    mid = borders['eta_mid'][borders['layer']==layer].iloc[0]
    high = borders['eta_high'][borders['layer']==layer].iloc[0]
    
    for index in range(len(towerEtaLines)-1):
        
        if ( (towerEtaLines[index] < mid) and (towerEtaLines[index+1] > low )):
            highEtaEdge=min(towerEtaLines[index+1], mid)
            lowEtaEdge=max(towerEtaLines[index], low)
            slice_u1[index] = (highEtaEdge - lowEtaEdge) * weight(highEtaEdge, lowEtaEdge, noWeight=False) * 1000
        
        if ( (towerEtaLines[index] < high) and (towerEtaLines[index+1] > mid )):
            highEtaEdge=min(towerEtaLines[index+1], high)
            lowEtaEdge=max(towerEtaLines[index], mid)
            slice_u0[index] = (highEtaEdge - lowEtaEdge) * weight(highEtaEdge, lowEtaEdge, noWeight=False) * 1000
    
    luSlices['l'+str(layer)+'-u0'] = slice_u0
    luSlices['l'+str(layer)+'-u1'] = slice_u1


for luSlice in luSlices:
    sort_index = np.argsort(luSlices[luSlice])#save indices before sorting
    
    luSliceSortedNormed = sortAndNormalize(luSlices[luSlice], half_N_div) #returns 1D np array of float type with sum=half_N_div
                        
    bestFit, isDegenerate = findBestFit(luSliceSortedNormed, half_N_div) #returns 1D np array of integer type with sum=half_N_div
    if (isDegenerate):
        print('check: degenerate fit result!')
        print(luSlice)
        print(20*'-')
    sliceFit = np.zeros(len(luSlices[luSlice]))
    for fit_index in range(len(bestFit)):#undo sort (retrieve original index)
        sliceFit[ sort_index[-1 - fit_index] ] = bestFit[-1 - fit_index]
    
    luSlicesFit[luSlice] = sliceFit


with open('input/allocation/allocation_20210421.txt') as f:
    lines = [line.rstrip('\n') for line in f]
    f.close()

bundlesScint = getModulesPerBundle(lines, isScintil=True)

modules = []
for bundle in bundlesScint:
    modules += bundlesScint[bundle]

towers = []
for towerPhi in range(24):
    for towerEta in range(-2,10):
        towers.append('had-eta'+str(towerEta)+'-phi'+str(towerPhi))

parMtx = pd.DataFrame(0, index=towers, columns=modules)

for module in modules:
    shares = luSlicesFit[module[:6]]
    for index in range(len(shares)):
        parMtx.at['had-eta'+str(index-2)+'-phi'+str(2*find_v(module)), module] = shares[index]
        parMtx.at['had-eta'+str(index-2)+'-phi'+str(1+2*find_v(module)), module] = shares[index]

parMtx_perBundle = {}

for i in range(len(bundlesScint)):
    parMtx_perBundle[i] = parMtx[parMtx.columns.intersection(bundlesScint[i])]




