import pandas as pd
import numpy as np
import ROOT
import math
from  root_numpy import hist2array, array2hist
from funcs import partitions, chisquare
#nBinsPhi = 32 #ROverZ_silicon_6_6_1->GetNbinsY() = 32
#minPhi = -0.34906585 #ROverZ_silicon_6_6_1->ProjectionY()->GetBinLowEdge(1) = -0.34906585
#maxPhi = 2.4434610 #ROverZ_silicon_6_6_1->ProjectionY()->GetBinLowEdge(1+ROverZ_silicon_6_6_1->GetNbinsY()) = 2.4434610
#
#nBinsROverZ = 42
#minROverZ = 0.076
#maxROverZ = 0.58

N_div = 8

bin_step = 2*math.pi/72
minEta = 16 * bin_step
maxEta = 38 * bin_step
nBinsEta = 38 - 16

minPhi = -4 * bin_step 
maxPhi = 26 * bin_step
nBinsPhi = 26 - (-4)

TCs = pd.read_csv('TCPositions/TCPositions_Zminus_siliconOnly.csv', sep=' ')
tower = {}
fit_TC_hist = {}
#tower_bestfit = {}
tower_array = {}
tower_flat_array = {}
tower_flat_array_highest = {}
tower_bestfit_flat_array = {}
inclusive = ROOT.TH2D("inclusive","",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)

first_max_Overlap = ROOT.TH1D("max_Overlap","",50,0,50)
second_max_Overlap = ROOT.TH1D("second_max_Overlap","",25,0,25)
third_max_Overlap = ROOT.TH1D("third_max_Overlap","",20,0,20)

nonZero_Overlap = ROOT.TH1D("nonZero_Overlap","",30,0,30)
nonZero_OverlapVsEta = ROOT.TH2D("nonZero_OverlapVsEta","",30,0,30,nBinsEta,minEta,maxEta)

for u in range(np.max(TCs['waferu'])): #wafer u
    for v in range(np.max(TCs['waferv'])): #wafer v
        #for l in range(1, 2): #layer number
        for l in range(1, int(np.max(TCs['layer'])) ): #layer number
            tower[u,v,l] = ROOT.TH2D("tower_u"+str(u)+"_v"+str(v)+"_layer"+str(l),"",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
            
            for index, row in TCs[(TCs["waferu"]==u) & (TCs["waferv"]==v) & (TCs["layer"]==l)].iterrows():
                tower[u,v,l].Fill(-1.0*row["triggercelleta"], row["triggercellphi"])
            
            tower_array[u,v,l] = hist2array(tower[u,v,l])
            tower_flat_array[u,v,l] = tower_array[u,v,l].flatten()
            
            fit_TC = np.zeros(len(tower_flat_array[u,v,l]))
            fit_TC_hist[u,v,l] = ROOT.TH2D("fitTC_u"+str(u)+"_v"+str(v)+"_layer"+str(l),"",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
            sort_index = np.argsort(tower_flat_array[u,v,l])
            
            tower_flat_array[u,v,l].sort()
            
            first_max_Overlap.Fill(tower_flat_array[u,v,l][-1])
            second_max_Overlap.Fill(tower_flat_array[u,v,l][-2])
            third_max_Overlap.Fill(tower_flat_array[u,v,l][-3])
            
            nonZero_Overlap_numOfTowers = len(tower_flat_array[u,v,l][tower_flat_array[u,v,l]!=0])
            if (nonZero_Overlap_numOfTowers>20):
                print("numTowers=", nonZero_Overlap_numOfTowers, ", u=", u, ", v=", v, ", l=", l)
            nonZero_Overlap.Fill(nonZero_Overlap_numOfTowers)
            nonZero_OverlapVsEta.Fill(nonZero_Overlap_numOfTowers,\
                                 TCs[(TCs["waferu"]==u) & (TCs["waferv"]==v) & (TCs["layer"]==l)]["triggercelleta"].mean()*-1.0)

            tower_flat_array_highest[u,v,l] = tower_flat_array[u,v,l][-1 * N_div:] #keep N_div highest values
            tower_flat_array_highest[u,v,l] = tower_flat_array_highest[u,v,l][tower_flat_array_highest[u,v,l]!=0] #remove zeros
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
                    fit_TC[ sort_index[-1 - fit_index] ] = tower_bestfit_flat_array[u,v,l][-1 - fit_index]
                
                fit_TC = fit_TC.reshape(nBinsEta, nBinsPhi)
                _ = array2hist (fit_TC, fit_TC_hist[u,v,l])



            #tower_bestfit[u,v,l] = ROOT.TH2D("tower_u"+str(u)+"_v"+str(v)+"_layer"+str(l),"",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
            


            inclusive.Add(tower[u,v,l])

