import pandas as pd
import numpy as np
import ROOT
import sys
import math
import yaml
from  root_numpy import hist2array, array2hist
from funcs import partitions, chisquare, calc_kernel

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

    cells = pd.read_csv('cellPositions/sensorCell_positions.txt', sep=' ') 
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

    if (config['debugging']):
        first_max_Overlap = ROOT.TH1D("max_Overlap","",45,0,45)
        second_max_Overlap = ROOT.TH1D("second_max_Overlap","",45,0,45)
        third_max_Overlap = ROOT.TH1D("third_max_Overlap","",45,0,45)
        nonZero_Overlap = ROOT.TH1D("nonZero_Overlap","",45,0,45)
        nonZero_OverlapVsEta = ROOT.TH2D("nonZero_OverlapVsEta","",nBinsEta,minEta,maxEta,45,0,45)
        counter_degen = 0 
        counter_total = 0 

    tower = {}
    fit_TC_hist = {}
    fit_TC = {}
    #tower_bestfit = {}
    tower_array = {}
    tower_array_gausKernel = {}
    tower_array_gausKernel_hist = {}
    tower_flat_array = {}
    tower_flat_array_highest = {}
    tower_bestfit_flat_array = {}


    kernel = calc_kernel(5)




    with open("splitModuleSumsOverTowers.txt", "w") as f:
        f.write("layer waferu waferv eta-phi-fraction\n")
    
        for l in range(1, 1+int(np.max(SCs['layer'])) ): #layer number
    #    for l in range(1, 2): #layer number
            for u in range(1+np.max(SCs['waferu'])): #wafer u
                for v in range(1+np.max(SCs['waferv'])): #wafer v
                    tower[u,v,l] = ROOT.TH2D("tower_u"+str(u)+"_v"+str(v)+"_layer"+str(l),"",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
                    
                    wafer_data = SCs[(SCs["waferu"]==u) & (SCs["waferv"]==v) & (SCs["layer"]==l)] 
                    if len(wafer_data)!=0:
                        counter_total += 1
        
                        for index, row in wafer_data.iterrows():
                            tower[u,v,l].Fill(-1.0*row["SC_eta"], row["SC_phi"])
                        
                        tower_array[u,v,l] = hist2array(tower[u,v,l])
                        
                        
                        #NEW
                        tower_array_gausKernel_hist[u,v,l] = ROOT.TH2D("fitTCKernel_u"+str(u)+"_v"+str(v)+"_layer"+str(l),"",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
                        tower_array_gausKernel[u,v,l] = ndimage.correlate(tower_array[u,v,l], kernel, mode='constant', cval = 0.0)  
                        _ = array2hist (tower_array_gausKernel[u,v,l], tower_array_gausKernel_hist[u,v,l])
                        tower_flat_array[u,v,l] = tower_array_gausKernel[u,v,l].flatten()
                        #New finish
    
                        
                        #tower_flat_array[u,v,l] = tower_array[u,v,l].flatten()
                        
                        fit_TC[u,v,l] = np.zeros(len(tower_flat_array[u,v,l]))
                        fit_TC_hist[u,v,l] = ROOT.TH2D("fitTC_u"+str(u)+"_v"+str(v)+"_layer"+str(l),"",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
                        sort_index = np.argsort(tower_flat_array[u,v,l])
                        
                        tower_flat_array[u,v,l].sort()
                        
                        #first_max_Overlap.Fill(tower_flat_array[u,v,l][-1])
                        #second_max_Overlap.Fill(tower_flat_array[u,v,l][-2])
                        #third_max_Overlap.Fill(tower_flat_array[u,v,l][-3])
                        
                        nonZero_Overlap_numOfTowers = len(tower_flat_array[u,v,l][tower_flat_array[u,v,l]!=0])
                        if (nonZero_Overlap_numOfTowers>125):
                            print("numTowers=", nonZero_Overlap_numOfTowers, ", u=", u, ", v=", v, ", l=", l)
                        #nonZero_Overlap.Fill(nonZero_Overlap_numOfTowers)
                        #nonZero_OverlapVsEta.Fill(wafer_data["SC_eta"].mean()*-1.0, nonZero_Overlap_numOfTowers)
                        #if len(wafer_data)==9*48:
                        #    nonZero_OverlapVsEta_fullsHighDensity.Fill(wafer_data["SC_eta"].mean()*-1.0, nonZero_Overlap_numOfTowers)
        
                        if ((tower_flat_array[u,v,l][-1 * N_div] == tower_flat_array[u,v,l][-1 * N_div - 1])\
                                    and (tower_flat_array[u,v,l][-1 * N_div])>tower_flat_array[u,v,l][-1]/100.0):
                            print("--------------------------------------------------")
                            print("degeneracy in choosing the", N_div, "'th highest value.")
                            print("u=", u, ", v=", v, ", l=", l)
                            print(tower_flat_array[u,v,l][-3*N_div:])
                            print("--------------------------------------------------")
                            counter_degen += 1
                        
                        if ((tower_flat_array[u,v,l][-1 * N_div] == tower_flat_array[u,v,l][-1 * N_div + 1])\
                                    and (tower_flat_array[u,v,l][-1 * N_div])>tower_flat_array[u,v,l][-1]/100.0):
                            print("++++++++++++++++++++++++++++++++++++++++++++++++++")
                            print("degeneracy in choosing the", N_div - 1, "'th highest value.")
                            print("u=", u, ", v=", v, ", l=", l)
                            print(tower_flat_array[u,v,l][-3*N_div:])
                            print("++++++++++++++++++++++++++++++++++++++++++++++++++")
                            counter_degen += 1
        
                        tower_flat_array_highest[u,v,l] = tower_flat_array[u,v,l][-1 * N_div:] #keep N_div highest values
                        tower_flat_array_highest[u,v,l] = tower_flat_array_highest[u,v,l][tower_flat_array_highest[u,v,l]>tower_flat_array_highest[u,v,l][-1]/100.] #remove small numbers
                        factor = tower_flat_array_highest[u,v,l].sum()
                        tower_flat_array_highest[u,v,l] = (tower_flat_array_highest[u,v,l]/factor)*N_div
                                
                        
                        chi2_min = 1000000.0
                        if len(tower_flat_array_highest[u,v,l])!=0:
                            for p in partitions(N_div, len(tower_flat_array_highest[u,v,l])):
                                chi2_temp = chisquare(p,tower_flat_array_highest[u,v,l])
                                if (chi2_temp < chi2_min):
                                    chi2_min = chi2_temp
                                    tower_bestfit_flat_array[u,v,l] = p 
                   
                            for fit_index in range(len(tower_bestfit_flat_array[u,v,l])):
                                fit_TC[u,v,l][ sort_index[-1 - fit_index] ] = tower_bestfit_flat_array[u,v,l][-1 - fit_index]
                            
                            fit_TC[u,v,l] = fit_TC[u,v,l].reshape(nBinsEta, nBinsPhi)
                            _ = array2hist (fit_TC[u,v,l], fit_TC_hist[u,v,l])
    
                        eta_phi_fit_TC = [[l[0]-1, l[1]-7] for l in np.transpose(np.nonzero(fit_TC[u,v,l])).tolist()]
                        values_fit_TC = fit_TC[u,v,l][np.nonzero(fit_TC[u,v,l])].astype(int)
        
                        f.write("{} {} {} {}\n".format(l, u, v, tuple(zip(eta_phi_fit_TC, values_fit_TC))))
        
                        #tower_bestfit[u,v,l] = ROOT.TH2D("tower_u"+str(u)+"_v"+str(v)+"_layer"+str(l),"",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
                        
                        inclusive_fit_TC.Add(fit_TC_hist[u,v,l])
                        inclusive.Add(tower[u,v,l])







main()
