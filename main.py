print("start importing main")
import ROOT
import sys
import math
import yaml
from  root_numpy import hist2array, array2hist
from funcs import partitions, chisquare, calc_kernel, getModulesPerBundle,\
                    getParMtxPerBundle, writeParMtxPerBundleToFile, writeTowerPerModuleToFile,\
                    getModulesWithTC, applyKernel, sortAndNormalize, findBestFit, isDegenerate
print("start importing np in main")
import numpy as np
print("start importing pd in main")
import pandas as pd
print("finish importing main")

def param_mtx(inputdir, SC_position_file, outputdir, param_mtx_em_name, param_mtx_had_name,\
                inputdir_bundlefile, bundles_file_path):

    cells = pd.read_csv(inputdir + SC_position_file , sep=' ') 
    cells.columns= ["layer","waferu","waferv","triggercellu","triggercellv","SC_eta","SC_phi"]
    
    cells = cells[(cells.layer % 2 == 1) | (cells.layer >28)].reset_index(drop=True)#Only use trigger layers. 
    cells["SC_phi"] = cells["SC_phi"].replace(0, 1e-5) #Force SCs on border phi=0 to fill positive phi bins.
    
    N_div = 1
    
    etaBinStep = 0.0870
    minEta = 16 * etaBinStep
    maxEta = 38 * etaBinStep
    nBinsEta = 38 - 16
    
    phiBinStep = 2*math.pi/72
    minPhi = -7 * phiBinStep
    maxPhi = 30 * phiBinStep
    nBinsPhi = 30 - (-7)
        
    inclusive_fit_TC = ROOT.TH2D("inclusive_fit_TC","",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
    inclusive_numOfModulesPerTower = ROOT.TH2D("inclusive_numOfModulesPerTower","",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
    
    tower = {}
    towerFit = {}
    tower_bestfit_flat_array = {}
    numOfModulesPerTower = {}
    
    kernel = calc_kernel(5) #nxn matrix for smoothening. The input(i.e. n) should be odd.
    
    param_mtx_em = pd.DataFrame() #paramter matrix (module vs tower) for CE-E
    param_mtx_had = pd.DataFrame() #paramter matrix (module vs tower) for CE-H
    param_mtx = {0:param_mtx_em, 1:param_mtx_had}
    
    modulesWithTC = getModulesWithTC(inputdir_bundlefile + bundles_file_path)
                    #Some partial modules have SC but not TC ('c' shaped). The line below finds modules with TC

    for l in range(1, 1+int(np.max(cells['layer'])) ): #layer number
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
                    
                    ##########
                    ###tower_flat_array[u,v,l] = towerSmoothed.flatten()
                    
                    ###fit_TC[u,v,l] = np.zeros(len(tower_flat_array[u,v,l]))
                    ###fit_TC_hist[u,v,l] = ROOT.TH2D("fitTC_u"+str(u)+"_v"+str(v)+"_layer"+str(l),"",\
                    ###                        nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
                    
                    ###sort_index = np.argsort(tower_flat_array[u,v,l]) #save positions before sorting
                    ###tower_flat_array[u,v,l].sort()
                    
                    ###if (tower_flat_array[u,v,l][-1 * N_div] == tower_flat_array[u,v,l][-1 * N_div -1]):
                    ###    print('check: degenerate!')
                    ###tower_flat_array_highest[u,v,l] = tower_flat_array[u,v,l][-1 * N_div:] #keep N_div highest values
                    ###tower_flat_array_highest[u,v,l] = tower_flat_array_highest[u,v,l][tower_flat_array_highest[u,v,l]!=0] 
                                                    #remove zeros prevents possiblity giving energy to towers with zero overlap.
                    ###factor = tower_flat_array_highest[u,v,l].sum()
                    ###tower_flat_array_highest[u,v,l] = (tower_flat_array_highest[u,v,l]/factor)*N_div
                            
                    ###if len(tower_flat_array_highest[u,v,l])==0:
                    ###    print(20*'*' + 'ERROR: Should not be the case?: ' + 20*'*')
                    ###    sys.exit(1)
                   ########### 
                    ###chi2_min = 1000000.0 #a very large number
                    ###for p in partitions(N_div, len(tower_flat_array_highest[u,v,l])):
                    ###    chi2_temp = chisquare(p,tower_flat_array_highest[u,v,l])
                    ###    if (chi2_temp < chi2_min):
                    ###        chi2_min = chi2_temp
                    ###        tower_bestfit_flat_array[u,v,l] = p 
                    
                    ###FIXME CHECK DEGENERATE SOLUTIONS        
                    
                    towerFit_array = np.zeros(len(towerSmoothed.flatten()))
                    towerFit[u,v,l] = ROOT.TH2D("fitTC_u"+str(u)+"_v"+str(v)+"_layer"+str(l),"",\
                                            nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
               
                    for fit_index in range(len(tower_bestfit_flat_array[u,v,l])):#undo sort (retrieve original index)
                        towerFit_array[ sort_index[-1 - fit_index] ] = tower_bestfit_flat_array[u,v,l][-1 - fit_index]
                    
                    towerFit_array = towerFit_array.reshape(nBinsEta, nBinsPhi)
                    _ = array2hist (towerFit_array, towerFit[u,v,l])
    
                    numOfModulesPerTower[u,v,l] = ROOT.TH2D("numOfModulesPerTower_u"\
                                                    +str(u)+"_v"+str(v)+"_layer"+str(l),"",\
                                                    nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
                    _ = array2hist((towerFit_array!=0).astype(int), numOfModulesPerTower[u,v,l])#(array!=0).astype(int) includes 0 & 1 only
    
                    eta_phi_fit_TC = [[m[0]-1, m[1]-7] for m in np.transpose(np.nonzero(towerFit_array)).tolist()]
                    values_fit_TC = towerFit_array[np.nonzero(towerFit_array)].astype(int)
                    NumTowersOverlapModule = len(eta_phi_fit_TC)
    
                    ####################DataFrame#######################
                    isHad = l>28 # True if in CE-H
                    colName = 'l' + str(l) + '-u' + str(u) + '-v' +str(v) # name of new column
                    param_mtx[isHad].insert(len(param_mtx[isHad].columns), colName, np.zeros(len(param_mtx[isHad])))
                    prefixRowName = 'had' if isHad else 'em'
                    for idx in range(NumTowersOverlapModule):
                        RowName = prefixRowName + '-eta' + str(eta_phi_fit_TC[idx][0]) + '-phi' + str(eta_phi_fit_TC[idx][1])
                        if (not RowName in param_mtx[isHad].index):
                            param_mtx[isHad].loc[RowName] = np.zeros(len(param_mtx[isHad].columns))
                        param_mtx[isHad].at[RowName, colName] = values_fit_TC[idx]
                    ######################DataFrame#####################
                    
                    inclusive_fit_TC.Add(towerFit[u,v,l])
                    inclusive_numOfModulesPerTower.Add(numOfModulesPerTower[u,v,l])
    
    param_mtx[0].to_pickle(outputdir + param_mtx_em_name)
    param_mtx[1].to_pickle(outputdir + param_mtx_had_name)
    
    return inclusive_numOfModulesPerTower

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
        inclusive_numOfModulesPerTower = param_mtx(inputdir=config['param_mtx']['inputdir'], \
                  SC_position_file=config['param_mtx']['SC_position_file'],\
                  outputdir=config['param_mtx']['outputdir'], \
                  param_mtx_em_name=config['param_mtx']['param_mtx_em_name'],\
                  param_mtx_had_name=config['param_mtx']['param_mtx_had_name'], \
                  inputdir_bundlefile=config['module_per_tower']['inputdir'],\
                  bundles_file_path=config['module_per_tower']['bundles_file']\
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
    
    return inclusive_numOfModulesPerTower



if __name__ == "__main__":
    print('main started!')
    inclusive_numOfModulesPerTower = main()
