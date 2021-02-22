import ROOT
import numpy as np
import pandas as pd
import sys
import math
import yaml
from  root_numpy import hist2array, array2hist
from funcs import partitions, chisquare, calc_kernel, getModulesPerBundle,\
                    getParMtxPerBundle, writeParMtxPerBundleToFile, writeTowerPerModuleToFile,\
                    getModulesWithTC, applyKernel, sortAndNormalize, findBestFit, isDegenerate

def param_mtx(inputdir, SC_position_file, outputdir, param_mtx_em_name, param_mtx_had_name,\
                inputdir_bundlefile, bundles_file_path, do2DHists):

    cells = pd.read_csv(inputdir + SC_position_file , sep=' ') 
    cells.columns= ["layer","waferu","waferv","triggercellu","triggercellv","SC_eta","SC_phi"]
    
    cells = cells[(cells.layer % 2 == 1) | (cells.layer >28)].reset_index(drop=True)#Only use trigger layers. 
    cells["SC_phi"] = cells["SC_phi"].replace(0, 1e-5) #Force SCs on border phi=0 to fill positive-phi bins.
    
    N_div = 8 # Divide module sum to (1/N_div)'s
    
    etaBinStep = 0.0870
    minBinEta = 16
    maxBinEta = 38
    minEta = minBinEta * etaBinStep
    maxEta = maxBinEta * etaBinStep
    nBinsEta = maxBinEta - minBinEta
    
    phiBinStep = 2*math.pi/72
    minBinPhi = -7
    maxBinPhi = 30
    minPhi = minBinPhi * phiBinStep
    maxPhi = maxBinPhi * phiBinStep
    nBinsPhi = maxBinPhi - minBinPhi
        
    
    tower = {}#key: tuple showing module ID. value: ROOT 2D hist showing how many SCs per tower
    towerFit = {} #key: tuple showing module ID. value: ROOT 2D hist showing how module sum split to (1/N_div)
    inclusive_towerFit = ROOT.TH2D("inclusive_towerFit","",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
    inclusive_numOfModulesPerTower = ROOT.TH2D("inclusive_numOfModulesPerTower","",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
    
    kernel = calc_kernel(5) #input=n; then kernel is nxn matrix for smoothening. n should be odd.
    
    param_mtx_em = pd.DataFrame() #paramter matrix (module vs tower) for CE-E
    param_mtx_had = pd.DataFrame() #paramter matrix (module vs tower) for CE-H
    param_mtx = {0:param_mtx_em, 1:param_mtx_had}
    
    modulesWithTC = getModulesWithTC(inputdir_bundlefile + bundles_file_path)
                    #Some partial modules have SC but not TC ('c' shaped). The line below finds modules with TC

    #for l in range(1, 1+int(np.max(cells['layer'])) ): #layer number
    for l in [1, 33]: #layer number
        if (l <= 28 and l%2 == 0): #only using trigger layers 
            continue
        print('layer= ', l)
        for u in range(1+np.max(cells['waferu'])): #wafer u
            for v in range(1+np.max(cells['waferv'])): #wafer v
                tower[u,v,l] = ROOT.TH2D("tower_u"+str(u)+"_v"+str(v)+"_layer"+str(l),"",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
                
                wafer_data = cells[(cells["waferu"]==u) & (cells["waferv"]==v) & (cells["layer"]==l)] 
                if (len(wafer_data)!=0) and ('l'+str(l)+'-u'+str(u)+'-v'+str(v) in modulesWithTC):
    
                    for index, row in wafer_data.iterrows():
                        tower[u,v,l].Fill(-1.0*row["SC_eta"], row["SC_phi"])#2D hist of the number of SC
                    
                    towerSmoothed = applyKernel(tower[u,v,l], kernel)#tower smoothed with kernel. returns np 2D array 
                    sort_index = np.argsort(towerSmoothed.flatten())#save indices before sorting
                    towerSortedNormed = sortAndNormalize(towerSmoothed, N_div) #returns 1D np array of float type with sum=N_div
                    bestFit, chi2_min = findBestFit(towerSortedNormed, N_div) #returns 1D np array of integer type with sum=N_div
                    if (isDegenerate(towerSortedNormed, N_div, chi2_min)):
                        print('check: degenerate fit result!')
                        print('l=', l, 'u=', u, "v=", v)
                        print(20*'-')
                    
                    towerFit_array = np.zeros(len(towerSmoothed.flatten()))
                    for fit_index in range(len(bestFit)):#undo sort (retrieve original index)
                        towerFit_array[ sort_index[-1 - fit_index] ] = bestFit[-1 - fit_index]
                    towerFit_array = towerFit_array.reshape(nBinsEta, nBinsPhi)
                    OverlapTowerCoord = [[m[0]-1, m[1]-7] for m in np.transpose(np.nonzero(towerFit_array)).tolist()] 
                                        #eta & phi coordinates of towers overlapping the module. 1 and 7 are just offset in eta and phi
                    OverlapTowerShare = towerFit_array[np.nonzero(towerFit_array)].astype(int) 
                                        #array of integers, ranging 1 to N_div. Shows the share each tower gets from module sum.
    
                    ####################Adding to dataframe#######################
                    isHad = l>28 # True if in CE-H
                    colName = 'l' + str(l) + '-u' + str(u) + '-v' +str(v) # name of new column
                    param_mtx[isHad].insert(len(param_mtx[isHad].columns), colName, np.zeros(len(param_mtx[isHad])))
                    prefixRowName = 'had' if isHad else 'em'
                    NumTowersOverlapModule = len(OverlapTowerCoord)
                    for idx in range(NumTowersOverlapModule):
                        RowName = prefixRowName + '-eta' + str(OverlapTowerCoord[idx][0]) + '-phi' + str(OverlapTowerCoord[idx][1])
                        if (not RowName in param_mtx[isHad].index):
                            param_mtx[isHad].loc[RowName] = np.zeros(len(param_mtx[isHad].columns))
                        param_mtx[isHad].at[RowName, colName] = OverlapTowerShare[idx]
                    ##############################################################

                    if (do2DHists):
                        #hist: how many (1/N_div)'s each tower gets from a module
                        towerFit[u,v,l] = ROOT.TH2D("towerFit_u"+str(u)+"_v"+str(v)+"_layer"+str(l),"",\
                                                nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
                        _ = array2hist (towerFit_array, towerFit[u,v,l])
                        inclusive_towerFit.Add(towerFit[u,v,l]) #inclusive all layers
                        
                        #hist: which towers a module overlaps with (i.e. fills with 0 or 1). "Inclusive" shows how many sums needed
                        numOfModulesPerTower = ROOT.TH2D("numOfModulesPerTower_u"\
                                                        +str(u)+"_v"+str(v)+"_layer"+str(l),"",\
                                                        nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
                        _ = array2hist((towerFit_array!=0).astype(int), numOfModulesPerTower)
                                            #(array!=0).astype(int) includes 0 & 1 only
                        inclusive_numOfModulesPerTower.Add(numOfModulesPerTower) #inclusive all layers
    
    param_mtx[0].to_pickle(outputdir + param_mtx_em_name)
    param_mtx[1].to_pickle(outputdir + param_mtx_had_name)
    
    return tower, towerFit, inclusive_towerFit, inclusive_numOfModulesPerTower

def module_per_tower(inputdir, outputdir, bundles_file_path, inputdir_paramMtx, param_mtx_em_name, param_mtx_had_name):
    with open(inputdir + bundles_file_path) as f:
        lines = [line.rstrip('\n') for line in f]
    f.close()
    bundles = getModulesPerBundle(lines)
    parMtxEM_PerBundle, parMtxHad_PerBundle = getParMtxPerBundle(bundles, inputdir_paramMtx, param_mtx_em_name, param_mtx_had_name)
    writeParMtxPerBundleToFile(outputdir, parMtxEM_PerBundle, name='CE-E')
    writeParMtxPerBundleToFile(outputdir, parMtxHad_PerBundle, name='CE-H')

def tower_per_module(outputdir, inputdir_paramMtx, param_mtx_em_name, param_mtx_had_name):
    parMtxEM = pd.read_pickle(inputdir_paramMtx + param_mtx_em_name).astype('int')
    parMtxHad = pd.read_pickle(inputdir_paramMtx + param_mtx_had_name).astype('int')
    writeTowerPerModuleToFile(outputdir, parMtxEM, parMtxHad)

def main():

    try:
        config_file = sys.argv[1]
    except IndexError:
        print("Please give a valid config file")
        exit()
    try:
        with open(config_file, 'r') as file:
            config = yaml.load(file, Loader=yaml.FullLoader)
    except EnvironmentError:
        print("Please give a valid config file")
        exit()
    
    if (config['function']['param_mtx']):
        tower, towerFit, inclusive_towerFit, inclusive_numOfModulesPerTower = param_mtx(inputdir=config['param_mtx']['inputdir'], \
                  SC_position_file=config['param_mtx']['SC_position_file'],\
                  outputdir=config['param_mtx']['outputdir'], \
                  param_mtx_em_name=config['param_mtx']['param_mtx_em_name'],\
                  param_mtx_had_name=config['param_mtx']['param_mtx_had_name'], \
                  inputdir_bundlefile=config['module_per_tower']['inputdir'],\
                  bundles_file_path=config['module_per_tower']['bundles_file'],\
                  do2DHists=config['param_mtx']['do2DHists']\
                  )

    if (config['function']['tower_per_module']):
        tower_per_module(outputdir=config['tower_per_module']['outputdir'],\
                         inputdir_paramMtx=config['param_mtx']['outputdir'],\
                         param_mtx_em_name=config['param_mtx']['param_mtx_em_name'],\
                         param_mtx_had_name=config['param_mtx']['param_mtx_had_name']\
                         )

    if (config['function']['module_per_tower']):
        module_per_tower(inputdir=config['module_per_tower']['inputdir'],\
                         outputdir=config['module_per_tower']['outputdir'],\
                         bundles_file_path=config['module_per_tower']['bundles_file'],\
                         inputdir_paramMtx=config['param_mtx']['outputdir'],\
                         param_mtx_em_name=config['param_mtx']['param_mtx_em_name'],\
                         param_mtx_had_name=config['param_mtx']['param_mtx_had_name']\
                         )
    
    return tower, towerFit, inclusive_towerFit, inclusive_numOfModulesPerTower



if __name__ == "__main__":
    tower, towerFit, inclusive_towerFit, inclusive_numOfModulesPerTower = main()
