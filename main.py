import ROOT
import numpy as np
import pandas as pd
import sys
import time
import math
import yaml
from  root_numpy import hist2array, array2hist
from funcs import partitions, chisquare, calc_kernel, getModulesPerBundle,\
                    getParMtxPerBundle, writeParMtxPerBundleToFile, writeTowerPerModuleToFile,\
                    getModulesWithTC, applyKernel, sortAndNormalize, findBestFit, SaveHist

def param_mtx(inputdir, SC_position_file, outputdir, param_mtx_em_name, param_mtx_had_name,\
                inputdir_bundlefile, bundles_file_path, do2DHists):

    last_CE_E_layer = 28

    cells = pd.read_csv(inputdir + SC_position_file , sep=' ') 
    cells.columns= ["layer","waferu","waferv","triggercellu","triggercellv","SC_eta","SC_phi"]
    
    cells = cells[(cells.layer % 2 == 1) | (cells.layer > last_CE_E_layer)].reset_index(drop=True)#Only use trigger layers. 
    cells["SC_phi"] = cells["SC_phi"].replace(0, 1e-5) #Force SCs on border phi=0 to fill positive-phi bins.
    
    N_div = 1 # Divide module sum to (1/N_div)'s
    
    etaBinStep = 0.0870
    minBinEta = 16 #chosen conservatively for visualization
    maxBinEta = 38 #chosen conservatively for visualization 
    minEta = minBinEta * etaBinStep
    maxEta = maxBinEta * etaBinStep
    nBinsEta = maxBinEta - minBinEta
    
    phiBinStep = 2*math.pi/72
    minBinPhi = -7 #chosen conservatively for visualization 
    maxBinPhi = 30 #chosen conservatively for visualization 
    minPhi = minBinPhi * phiBinStep
    maxPhi = maxBinPhi * phiBinStep
    nBinsPhi = maxBinPhi - minBinPhi
    
    
    tower = ROOT.TH2D("tower","",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)#hist showing how much a module overlap with each tower
    towerFit = ROOT.TH2D("towerFit","",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)#hist showing how a module sum split to (1/N_div)
    numOfModulesPerTower = ROOT.TH2D("numOfModulesPerTower","",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
    
    #inclusives
    inclusive_tower = ROOT.TH2D("inclusive_tower","",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
    inclusive_towerFit = ROOT.TH2D("inclusive_towerFit","",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
    inclusive_numOfModulesPerTower = ROOT.TH2D("inclusive_numOfModulesPerTower"\
                                     ,"",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
    inclusive_numOfModulesPerTower_OnlyEM = ROOT.TH2D("inclusive_numOfModulesPerTower_OnlyEM"\
                                            ,"",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
    inclusive_numOfModulesPerTower_OnlyHad = ROOT.TH2D("inclusive_numOfModulesPerTower_OnlyHad"\
                                            ,"",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
    
    if(do2DHists): #save all! good for debugging, but slows down
        tower_saved = {}#key: tuple showing module ID. value: ROOT 2D hist showing how much the module overlap with each tower
        towerFit_saved = {} #key: tuple showing module ID. value: ROOT 2D hist showing how module sum split to (1/N_div)
    
    kernel = calc_kernel(5) #input=n; then kernel is nxn matrix for smoothening. n should be odd.
    
    param_mtx_em = pd.DataFrame() #paramter matrix (module vs tower) for CE-E
    param_mtx_had = pd.DataFrame() #paramter matrix (module vs tower) for CE-H
    param_mtx = {0:param_mtx_em, 1:param_mtx_had}
    
    modulesWithTC = getModulesWithTC(inputdir_bundlefile + bundles_file_path)
                    #Some partial modules have SC but not TC ('c' shaped). The line below finds modules with TC

    for l in range(1, 1+int(np.max(cells['layer'])) ): #layer number
        if (l <= last_CE_E_layer and l%2 == 0): #only using trigger layers 
            continue
        print('layer= ', l)
        for u in range(1+np.max(cells['waferu'])): #wafer u
            for v in range(1+np.max(cells['waferv'])): #wafer v
                wafer_data = cells[(cells["waferu"]==u) & (cells["waferv"]==v) & (cells["layer"]==l)] 
                
                if (len(wafer_data)!=0) and ('l'+str(l)+'-u'+str(u)+'-v'+str(v) in modulesWithTC):
                    tower.Reset()
                    for index, row in wafer_data.iterrows():
                        tower.Fill(-1.0*row["SC_eta"], row["SC_phi"])#2D hist of the number of SC
                    towerSmoothed = applyKernel(tower, kernel)#tower smoothed with kernel. returns np 2D array 
                    sort_index = np.argsort(towerSmoothed.flatten())#save indices before sorting
                    towerSortedNormed = sortAndNormalize(towerSmoothed, N_div) #returns 1D np array of float type with sum=N_div
                    bestFit, isDegenerate = findBestFit(towerSortedNormed, N_div) #returns 1D np array of integer type with sum=N_div
                    if (isDegenerate):
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
                    isHad = l > last_CE_E_layer # True if in CE-H
                    colName = 'l' + str(l) + '-u' + str(u) + '-v' +str(v) # name of new column
                    param_mtx[isHad].insert(len(param_mtx[isHad].columns), colName, np.zeros(len(param_mtx[isHad])))
                    prefixRowName = 'had' if isHad else 'em'
                    NumTowersOverlapModule = len(OverlapTowerCoord)
                    for idx in range(NumTowersOverlapModule):
                        RowName = prefixRowName + '-eta' + str(OverlapTowerCoord[idx][0]) + '-phi' + str(OverlapTowerCoord[idx][1])
                        if (not RowName in param_mtx[isHad].index):
                            param_mtx[isHad].loc[RowName] = np.zeros(len(param_mtx[isHad].columns))
                        param_mtx[isHad].at[RowName, colName] = OverlapTowerShare[idx]
                    
                    #########################Visualsing###########################
                    #hist: sum all SCs in tower coordinates
                    inclusive_tower.Add(tower)
                    
                    #hist: how many (1/N_div)'s each tower gets from a module
                    towerFit.Reset()
                    _ = array2hist (towerFit_array, towerFit)
                    inclusive_towerFit.Add(towerFit) #inclusive all layers
                    
                    #hist: which towers a module overlaps with (i.e. fills with 0 or 1). "Inclusive" shows how many sums needed
                    numOfModulesPerTower.Reset()
                    _ = array2hist((towerFit_array!=0).astype(int), numOfModulesPerTower)
                                        #(array!=0).astype(int) includes 0 & 1 only
                    inclusive_numOfModulesPerTower.Add(numOfModulesPerTower) #inclusive all layers
                    if (isHad):
                        inclusive_numOfModulesPerTower_OnlyHad.Add(numOfModulesPerTower) #inclusive all CE-H layers
                    else:
                        inclusive_numOfModulesPerTower_OnlyEM.Add(numOfModulesPerTower) #inclusive all CE-E layers
                    
                    if (do2DHists): #Save hists per module
                        #copy "tower" 2D hist
                        tower_saved[u,v,l] = ROOT.TH2D("tower_saved_u"+str(u)+"_v"+str(v)+"_layer"+str(l),\
                                            "",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
                        for index, row in wafer_data.iterrows():
                            tower_saved[u,v,l].Fill(-1.0*row["SC_eta"], row["SC_phi"])#2D hist of the number of SC

                        #copy "towerFit" 2D hist
                        towerFit_saved[u,v,l] = ROOT.TH2D("towerFit_saved_u"+str(u)+"_v"+str(v)+"_layer"+str(l),"",\
                                                nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
                        _ = array2hist (towerFit_array, towerFit_saved[u,v,l])

    param_mtx[0] = param_mtx[0].drop([x for x in param_mtx[0].index if (x[:8]=='em-eta18')])#remove eta>3.045 (Not removed from hists)
    param_mtx[1] = param_mtx[1].drop([x for x in param_mtx[1].index if (x[:9]=='had-eta18')])#remove eta>3.045 (Not removed from hists)
    param_mtx[0].to_pickle(outputdir + param_mtx_em_name)
    param_mtx[1].to_pickle(outputdir + param_mtx_had_name)
    
    inclusive_towerFit.Scale(1./N_div) #normalize
    SaveHist(inclusive_towerFit, outputdir+'/plots/', 'inclusive_towerFit_1Over'+str(N_div)+'s', 'root')
                                                                            #how the module sum (energy) is distributed
    SaveHist(inclusive_tower, outputdir+'/plots/', 'inclusive_tower_1Over'+str(N_div)+'s', 'root') #just to show SC occupation
    SaveHist(inclusive_numOfModulesPerTower, outputdir+'/plots/', \
                    'inclusive_numOfModulesPerTower_1Over'+str(N_div)+'s', 'root') #How many sums per tower
    SaveHist(inclusive_numOfModulesPerTower_OnlyHad, outputdir+'/plots/', \
                    'inclusive_numOfModulesPerTower_OnlyHad_1Over'+str(N_div)+'s', 'root') #How many sums per tower in CE-H
    SaveHist(inclusive_numOfModulesPerTower_OnlyEM, outputdir+'/plots/', \
                    'inclusive_numOfModulesPerTower_OnlyEM_1Over'+str(N_div)+'s', 'root') #How many sums per tower in CE-E

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
        param_mtx(inputdir=config['param_mtx']['inputdir'], \
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

if __name__ == "__main__":
    start = time.time()
    main()
    print('The program ran in', time.time() - start, 'seconds!')
