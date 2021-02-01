import pandas as pd
import numpy as np
import ROOT
import sys
import math
import yaml
from  root_numpy import hist2array, array2hist
from funcs import partitions, chisquare, calc_kernel, getModulesPerBundle
from scipy import ndimage

def param_mtx(inputdir, SC_position_file, outputdir, debugging):

    cells = pd.read_csv(inputdir + SC_position_file , sep=' ') 
    cells.columns= ["layer","waferu","waferv","triggercellu","triggercellv","SC_eta","SC_phi"]
     
    N_div = 8
    
    bin_step = 2*math.pi/72
    
    minEta = 16 * bin_step
    maxEta = 38 * bin_step
    nBinsEta = 38 - 16
    
    minPhi = -7 * bin_step 
    maxPhi = 30 * bin_step
    nBinsPhi = 30 - (-7)
        
    inclusive = ROOT.TH2D("inclusive","",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
    inclusive_fit_TC = ROOT.TH2D("inclusive_fit_TC","",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
    inclusive_numOfModulesPerTower = ROOT.TH2D("inclusive_numOfModulesPerTower","",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
    
    if (debugging):
        nonZero_Overlap = ROOT.TH1D("nonZero_Overlap","",45,0,45)
        nonZero_OverlapVsEta = ROOT.TH2D("nonZero_OverlapVsEta","",nBinsEta,minEta,maxEta,45,0,45)
    
    nonZero_Overlap_afterFit = ROOT.TH1D("nonZero_Overlap_afterFit","",10,0,10)
    
    tower = {}
    fit_TC_hist = {}
    fit_TC = {}
    #tower_bestfit = {}
    tower_array = {}
    tower_array_Kernel = {}
    tower_array_Kernel_hist = {}
    tower_flat_array = {}
    tower_flat_array_highest = {}
    tower_bestfit_flat_array = {}
    numOfModulesPerTower = {}
    numOfModulesPerTowerPerLayer = ROOT.TH1D("numOfModulesPerTowerPerLayer","",6,0,6)
    inclusive_numOfModulesPerTowerPerLayer = ROOT.TH2D("inclusive_numOfModulesPerTowerPerLayer","",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
    
    kernel = calc_kernel(5)
    
    param_mtx_em = pd.DataFrame() #paramter matrix (module vs tower) for CE-E
    param_mtx_had = pd.DataFrame() #paramter matrix (module vs tower) for CE-H
    param_mtx = {0:param_mtx_em, 1:param_mtx_had}
    
    #with open("./output/towers_per_module/splitModuleSumsOverTowers.txt", "w") as f:
        #f.write("layer waferu waferv numOfTowers ListOf:eta-phi-fraction\n")
    for l in range(1, 1+int(np.max(cells['layer'])) ): #layer number
        if (l <= 28 and l%2 == 0): #only using trigger layers 
            continue
        print('layer= ', l)
        inclusive_numOfModulesPerTowerPerLayer.Reset()
        for u in range(1+np.max(cells['waferu'])): #wafer u
            for v in range(1+np.max(cells['waferv'])): #wafer v
                tower[u,v,l] = ROOT.TH2D("tower_u"+str(u)+"_v"+str(v)+"_layer"+str(l),"",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
                
                wafer_data = cells[(cells["waferu"]==u) & (cells["waferv"]==v) & (cells["layer"]==l)] 
                if len(wafer_data)!=0:
    
                    for index, row in wafer_data.iterrows():
                        tower[u,v,l].Fill(-1.0*row["SC_eta"], row["SC_phi"])
                    
                    tower_array[u,v,l] = hist2array(tower[u,v,l])
                    
                    tower_array_Kernel_hist[u,v,l] = ROOT.TH2D("fitTCKernel_u"+str(u)+"_v"+str(v)+"_layer"+str(l),"",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
                    tower_array_Kernel[u,v,l] = ndimage.correlate(tower_array[u,v,l], kernel, mode='constant', cval = 0.0)  
                    for i in range(tower_array[u,v,l].shape[0]):
                        for j in range(tower_array[u,v,l].shape[1]):
                            if tower_array[u,v,l][i][j] == 0:
                                tower_array_Kernel[u,v,l][i][j] = 0
    
                    _ = array2hist (tower_array_Kernel[u,v,l], tower_array_Kernel_hist[u,v,l])
                    tower_flat_array[u,v,l] = tower_array_Kernel[u,v,l].flatten()
                    #tower_flat_array[u,v,l] = tower_array[u,v,l].flatten()
                    
                    fit_TC[u,v,l] = np.zeros(len(tower_flat_array[u,v,l]))
                    fit_TC_hist[u,v,l] = ROOT.TH2D("fitTC_u"+str(u)+"_v"+str(v)+"_layer"+str(l),"",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
                    sort_index = np.argsort(tower_flat_array[u,v,l])
                    
                    tower_flat_array[u,v,l].sort()
                    
                    if(debugging):
                        nonZero_Overlap_numOfTowers = len(tower_flat_array[u,v,l][tower_flat_array[u,v,l]!=0])
                        nonZero_Overlap.Fill(nonZero_Overlap_numOfTowers)
                        nonZero_OverlapVsEta.Fill(wafer_data["SC_eta"].mean()*-1.0, nonZero_Overlap_numOfTowers)
    
                        if tower_flat_array[u,v,l][-1 * N_div] == tower_flat_array[u,v,l][-1 * N_div - 1]\
                           and tower_flat_array[u,v,l][-1 * N_div]!=0:
                            print("--------------------------------------------------")
                            print("degeneracy in choosing the", N_div, "'th highest value.")
                            print("u=", u, ", v=", v, ", l=", l)
                            print(tower_flat_array[u,v,l][-3*N_div:])
                            print("--------------------------------------------------")
                        
                        if tower_flat_array[u,v,l][-1 * N_div] == tower_flat_array[u,v,l][-1 * N_div + 1]\
                           and tower_flat_array[u,v,l][-1 * N_div]!=0:
                            print("++++++++++++++++++++++++++++++++++++++++++++++++++")
                            print("degeneracy in choosing the", N_div - 1, "'th highest value.")
                            print("u=", u, ", v=", v, ", l=", l)
                            print(tower_flat_array[u,v,l][-3*N_div:])
                            print("++++++++++++++++++++++++++++++++++++++++++++++++++")
    
                    tower_flat_array_highest[u,v,l] = tower_flat_array[u,v,l][-1 * N_div:] #keep N_div highest values
                    tower_flat_array_highest[u,v,l] = tower_flat_array_highest[u,v,l][tower_flat_array_highest[u,v,l]!=0] #remove zeros prevents possiblity giving energy to towers with zero overlap.
                    factor = tower_flat_array_highest[u,v,l].sum()
                    tower_flat_array_highest[u,v,l] = (tower_flat_array_highest[u,v,l]/factor)*N_div
                            
                    if len(tower_flat_array_highest[u,v,l])==0:
                        print("Should not be the case?")
                    
                    chi2_min = 1000000.0
                    for p in partitions(N_div, len(tower_flat_array_highest[u,v,l])):
                        chi2_temp = chisquare(p,tower_flat_array_highest[u,v,l])
                        if (debugging):
                            if (abs(chi2_temp - chi2_min)<1.e-8):
                                print("DEGENERATE!!!!!!", " u=", u, ", v=", v, ", l=", l)
                                print(tower_flat_array[u,v,l][-3*N_div:])
                                print('-'*20)
                        if (chi2_temp < chi2_min):
                            chi2_min = chi2_temp
                            tower_bestfit_flat_array[u,v,l] = p 
               
                    for fit_index in range(len(tower_bestfit_flat_array[u,v,l])):
                        fit_TC[u,v,l][ sort_index[-1 - fit_index] ] = tower_bestfit_flat_array[u,v,l][-1 - fit_index]
                    
                    fit_TC[u,v,l] = fit_TC[u,v,l].reshape(nBinsEta, nBinsPhi)
                    _ = array2hist (fit_TC[u,v,l], fit_TC_hist[u,v,l])
    
                    numOfModulesPerTower[u,v,l] = ROOT.TH2D("numOfModulesPerTower_u"+str(u)+"_v"+str(v)+"_layer"+str(l),"",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
                    _ = array2hist((fit_TC[u,v,l]!=0).astype(int), numOfModulesPerTower[u,v,l])    
    
    
                    eta_phi_fit_TC = [[l[0]-1, l[1]-7] for l in np.transpose(np.nonzero(fit_TC[u,v,l])).tolist()]
                    values_fit_TC = fit_TC[u,v,l][np.nonzero(fit_TC[u,v,l])].astype(int)
                    NumTowersOverlapModule = len(eta_phi_fit_TC)
    
                    nonZero_Overlap_afterFit.Fill(NumTowersOverlapModule)
    
                    #f.write("{} {} {} ".format(l, u, v))#, tuple(zip(eta_phi_fit_TC, values_fit_TC))))
                    #f.write("{} ".format(NumTowersOverlapModule))
                    #for idx in range(NumTowersOverlapModule):
                    #    f.write("{} {} {} ".format(eta_phi_fit_TC[idx][0], eta_phi_fit_TC[idx][1], values_fit_TC[idx]))
                    #f.write("\n")
                    
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
                    ######################DataFrame Ends#####################
                    
                    inclusive_fit_TC.Add(fit_TC_hist[u,v,l])
                    inclusive.Add(tower[u,v,l])
                    inclusive_numOfModulesPerTowerPerLayer.Add(numOfModulesPerTower[u,v,l])
                    #if(numOfModulesPerTower[u,v,l].GetBinContent(2, 17)!=0):
                    #    print("4 towers", " etaphi: ", eta_phi_fit_TC, " u=", u, ", v=", v, ", l=", l)
                    inclusive_numOfModulesPerTower.Add(numOfModulesPerTower[u,v,l])
    
        for row in range(1, 1 + nBinsEta):
            for col in range(1, 1 + nBinsPhi):
                #if(inclusive_numOfModulesPerTowerPerLayer.GetBinContent(row, col)!=0):
                if(inclusive_numOfModulesPerTowerPerLayer.GetBinContent(row, col)==4):
                    numOfModulesPerTowerPerLayer.Fill(inclusive_numOfModulesPerTowerPerLayer.GetBinContent(row, col)) 
                    if(inclusive_numOfModulesPerTowerPerLayer.GetBinContent(row, col)==4):
                        print(row, col, l)
    
    
    
    param_mtx[0].to_pickle(outputdir + '/param_mtx_em.pkl')
    param_mtx[1].to_pickle(outputdir + '/param_mtx_had.pkl')

def module_per_tower(inputdir, outputdir, bundles_file_path):
    with open(inputdir + bundles_file_path) as f:
        lines = [line.rstrip('\n') for line in f]
    f.close()
    bundles = getModulesPerBundle(lines)
    #writeFileParamMtxPerBundle()


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
        param_mtx(inputdir=config['param_mtx']['inputdir'], SC_position_file=config['param_mtx']['SC_position_file'],\
                        outputdir=config['param_mtx']['outputdir'], debugging=config['debugging'])

    #if (config['function']['tower_per_module']):
    #    tower_per_module()

    if (config['function']['module_per_tower']):
        module_per_tower(inputdir=config['module_per_tower']['inputdir'], outputdir=config['module_per_tower']['outputdir'],\
                              bundles_file_path=config['module_per_tower']['bundles_file'])

if __name__ == "__main__":
    main()
