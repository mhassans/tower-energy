import ROOT
import numpy as np
import pandas as pd
import sys
import time
import math
import yaml
from  root_numpy import hist2array, array2hist
from funcs import partitions, chisquare, calc_kernel, getModulesPerBundle,\
                    getParMtxPerBundle_Silic, getParMtxPerBundle_Scint, writeParMtxPerBundleToFile, writeTowerPerModuleToFile,\
                    getModulesWithTC, applyKernel, sortAndNormalize, findBestFit, SaveHist,\
                    findBestFit, sortAndNormalize, weight, getModulesPerBundle, find_v

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()

def silicons(N_div, inputdir, SCsPosition_file, outputdir, param_mtx_em_name, param_mtx_had_name_silic,\
                inputdir_bundlefile, bundles_file_path, do2DHists):
    """
    Module sums will be split to (1/N_div)'s.
    """

    last_CE_E_layer = 28

    cells = pd.read_csv(inputdir + SCsPosition_file , sep=' ') 
    cells.columns= ["layer","waferu","waferv","triggercellu","triggercellv","SC_eta","SC_phi"]
    
    cells = cells[(cells.layer % 2 == 1) | (cells.layer > last_CE_E_layer)].reset_index(drop=True)#Only use trigger layers. 
    cells["SC_phi"] = cells["SC_phi"].replace(0, 1e-5) #Force SCs on border phi=0 to fill positive-phi bins.
    
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
    param_mtx[1].to_pickle(outputdir + param_mtx_had_name_silic)
    
    inclusive_towerFit.Scale(1./N_div) #normalize
    SaveHist(inclusive_towerFit, outputdir+'/plots/', 'inclusive_towerFit_1Over'+str(N_div)+'s', 'root', AddGrid=False)
                                                                            #how the module sum (energy) is distributed
    SaveHist(inclusive_tower, outputdir+'/plots/', 'inclusive_tower_1Over'+str(N_div)+'s', 'root', AddGrid=False) #just to show SC occupation
    SaveHist(inclusive_numOfModulesPerTower, outputdir+'/plots/', \
                    'inclusive_numOfModulesPerTower_1Over'+str(N_div)+'s', 'root', AddGrid=False) #How many sums per tower
    SaveHist(inclusive_numOfModulesPerTower_OnlyHad, outputdir+'/plots/', \
                    'inclusive_numOfModulesPerTower_OnlyHad_1Over'+str(N_div)+'s', 'root', AddGrid=False) #How many sums per tower in CE-H
    SaveHist(inclusive_numOfModulesPerTower_OnlyEM, outputdir+'/plots/', \
                    'inclusive_numOfModulesPerTower_OnlyEM_1Over'+str(N_div)+'s', 'root', AddGrid=False) #How many sums per tower in CE-E

def scintillators(N_div, TCsPosition_path, bundles_path, output_path):
    """
    Module sums will be split to (1/N_div)'s.
    """
    half_N_div = N_div//2 #Splitting over two phi slices (of 5 deg) assumed to be the same. So optmiziation performed on one only.
    etaBinStep = 0.0870
    
    TCs = pd.read_csv(TCsPosition_path, sep=' ')
    
    thresh ={
            37 : 1.57,
            38 : 1.58,
            39 : 1.59,
            40 : 1.60,
            41 : 1.79,
            42 : 1.81,
            43 : 1.83,
            44 : 1.85,
            45 : 1.86,
            46 : 1.88,
            47 : 1.90,
            48 : 1.92,
            49 : 1.93,
            50 : 1.95
            }
    
    eta_low = [] #(estimate for) the lowest-eta edge of scintillators (for all layers)
    eta_high = [] #(estimate for) highest-eta edge of scintillators (for all layers)
    eta_mid = [] #(estimate for) the border of u=0 and u=1 (for all layers)
    
    for j in range(37,51): #layers
        m = (TCs['triggercelleta']*-1)[(TCs['layer']==j) & (TCs['triggercelliphi']==1)]
        eta_low.append(pd.array(m)[-1] - (pd.array(m)[-2] - pd.array(m)[-1])/2) #2nd term added because last TC is not exactly on the border!
        eta_high.append(pd.array(m)[0] + (pd.array(m)[0] - pd.array(m)[1])/2) #2nd term added because last TC is not exactly on the border!
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
    
    
    
    luSlices = {} #Scints are divided into layer(l) and eta/u slices. No need to do it for phi/v as it is symmetrical. luSLices will be like: {'l50-u1': [0., 0.03, 3.50, 2.85, 2.33, 1.91, 1.58, 0.26, 0., 0. , 0. , 0.], ...}
    luSlicesFit = {} #shows the final allocation (Fit with integers summed to 8). Members will be like: {'l42-u1': [0, 2, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0],...}
    towerEtaLines = []
    
    for i in range(15, 28): #the eta range of scintillators
        towerEtaLines.append(i*etaBinStep)
    
    NormConst = 1000 # to avoid calculation on small numbers
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
                slice_u1[index] = (highEtaEdge - lowEtaEdge) * weight(highEtaEdge, lowEtaEdge, noWeight=False) * NormConst
            
            if ( (towerEtaLines[index] < high) and (towerEtaLines[index+1] > mid )):
                highEtaEdge=min(towerEtaLines[index+1], high)
                lowEtaEdge=max(towerEtaLines[index], mid)
                slice_u0[index] = (highEtaEdge - lowEtaEdge) * weight(highEtaEdge, lowEtaEdge, noWeight=False) * NormConst
        
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
    
    #using bundles' mapping-file ONLY for retrieving module names, e.g. 'scint-l46-u1-v6', 'scint-l49-u0-v3', ...
    with open(bundles_path) as f:
        lines = [line.rstrip('\n') for line in f]
        f.close()
    
    bundlesScint = getModulesPerBundle(lines, isScintil=True)
    
    modules = []
    for bundle in bundlesScint:
        modules += bundlesScint[bundle]
    towers = []
    for towerPhi in range(24): #tower phi bins in one sector
        for towerEta in range(-2,10): #tower eta bins. defined such that silicon layers start at towerEta=0.
            towers.append('had-eta'+str(towerEta)+'-phi'+str(towerPhi))
    
    parMtx = pd.DataFrame(0, index=towers, columns=modules)
    
    for module in modules:
        shares = luSlicesFit[module[6:12]]
        for index in range(len(shares)):
            parMtx.at['had-eta'+str(index-2)+'-phi'+str(2*find_v(module)), module] = shares[index]
            parMtx.at['had-eta'+str(index-2)+'-phi'+str(1+2*find_v(module)), module] = shares[index]
    
    parMtx.to_pickle(output_path)

def tower_per_module(N_div_silic, N_div_scint, outputdir, inputdir_paramMtx, param_mtx_em_name, param_mtx_had_name_silic, param_mtx_had_name_scint):
    parMtxEM = pd.read_pickle(inputdir_paramMtx + param_mtx_em_name).astype('int')
    parMtxHadSilic = pd.read_pickle(inputdir_paramMtx + param_mtx_had_name_silic).astype('int')
    parMtxHadScint = pd.read_pickle(inputdir_paramMtx + param_mtx_had_name_scint).astype('int')
    writeTowerPerModuleToFile(N_div_silic, N_div_scint, outputdir, parMtxEM, parMtxHadSilic, parMtxHadScint)

def module_per_tower(inputdir, outputdir, bundles_file_path, inputdir_paramMtx, param_mtx_em_name,\
                        param_mtx_had_name_silic, param_mtx_had_name_scint):
    with open(inputdir + bundles_file_path) as f:
        lines = [line.rstrip('\n') for line in f]
    f.close()
    
    bundlesSilic = getModulesPerBundle(lines, isScintil=False)
    bundlesScint = getModulesPerBundle(lines, isScintil=True)
    
    parMtxEM_PerBundle, parMtxHadSilic_PerBundle = getParMtxPerBundle_Silic(bundlesSilic, inputdir_paramMtx,\
                                                                        param_mtx_em_name, param_mtx_had_name_silic)
    parMtxHadScint_PerBundle = getParMtxPerBundle_Scint(bundlesScint, inputdir_paramMtx, param_mtx_had_name_scint)

    #merging silicon and scintillators:
    parMtxHad_PerBundle = {}
    for i in range(len(bundlesScint)):
        parMtxHadSilic_PerBundle[i] = parMtxHadSilic_PerBundle[i].multiply(2) #WARNING!!!\
                                                #To make all CE-H have the same denominator (16).Should\
                                                #be changed if the N_div_scint or N_div_silic changes in the config.
        parMtxHad_PerBundle[i] = pd.concat([parMtxHadSilic_PerBundle[i],parMtxHadScint_PerBundle[i]], axis=1).fillna(0).astype('int')
    print("==================================================")
    print("WARNING: In the current version, scintillators and silicons are divided by 1/16's and 1/8's,"\
           + " respectively. Energy shares of silicons in CE-H are multiplied by 2, so that all CE-E are"\
           + " divided by 8 and all CE-H are divided by 16 regardless of silicon/scint type."\
           + " IF THE DIVISORS ARE DIFFERENT, THE CODE NEEDS TO BE CHANGED. This change applied only in"\
           + " the module_per_tower function, as requested by Ante. No change is needed in other functions."\
           + " This change is currently hardcoded, should later be implemented in a clever way.")
    print("==================================================")
    
    writeParMtxPerBundleToFile(outputdir, parMtxEM_PerBundle, name='CE-E')
    writeParMtxPerBundleToFile(outputdir, parMtxHad_PerBundle, name='CE-H')

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
    
    if (config['mainFuncs']['silicons']):
        silicons(N_div=config['silicons']['N_div_silic'],\
                 inputdir=config['silicons']['inputdir'], \
                 SCsPosition_file=config['silicons']['SCsPosition_file'],\
                 outputdir=config['silicons']['outputdir'], \
                 param_mtx_em_name=config['silicons']['param_mtx_em_name'],\
                 param_mtx_had_name_silic=config['silicons']['param_mtx_had_name_silic'], \
                 inputdir_bundlefile=config['module_per_tower']['inputdir'],\
                 bundles_file_path=config['module_per_tower']['bundles_file'],\
                 do2DHists=config['silicons']['do2DHists']\
                 )

    if(config['mainFuncs']['scintillators']):
        scintillators(N_div=config['scintillators']['N_div_scint'],\
                      TCsPosition_path=config['silicons']['inputdir']+config['scintillators']['TCsPosition_file'],\
                      bundles_path=config['module_per_tower']['inputdir']+config['module_per_tower']['bundles_file'],\
                      output_path=config['silicons']['outputdir']+config['silicons']['param_mtx_had_name_scint']
                      )

    if (config['mainFuncs']['tower_per_module']):
        tower_per_module(N_div_silic=config['silicons']['N_div_silic'],\
                         N_div_scint=config['scintillators']['N_div_scint'],\
                         outputdir=config['tower_per_module']['outputdir'],\
                         inputdir_paramMtx=config['silicons']['outputdir'],\
                         param_mtx_em_name=config['silicons']['param_mtx_em_name'],\
                         param_mtx_had_name_silic=config['silicons']['param_mtx_had_name_silic'],\
                         param_mtx_had_name_scint=config['silicons']['param_mtx_had_name_scint']\
                         )

    if (config['mainFuncs']['module_per_tower']):
        module_per_tower(inputdir=config['module_per_tower']['inputdir'],\
                         outputdir=config['module_per_tower']['outputdir'],\
                         bundles_file_path=config['module_per_tower']['bundles_file'],\
                         inputdir_paramMtx=config['silicons']['outputdir'],\
                         param_mtx_em_name=config['silicons']['param_mtx_em_name'],\
                         param_mtx_had_name_silic=config['silicons']['param_mtx_had_name_silic'],\
                         param_mtx_had_name_scint=config['silicons']['param_mtx_had_name_scint']\
                         )

if __name__ == "__main__":
    start = time.time()
    main()
    print('The program ran in', time.time() - start, 'seconds!')
