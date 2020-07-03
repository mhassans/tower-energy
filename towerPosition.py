import pandas as pd
import numpy as np
import ROOT

nBinsPhi = 32 #ROverZ_silicon_6_6_1->GetNbinsY() = 32
minPhi = -0.34906585 #ROverZ_silicon_6_6_1->ProjectionY()->GetBinLowEdge(1) = -0.34906585
maxPhi = 2.4434610 #ROverZ_silicon_6_6_1->ProjectionY()->GetBinLowEdge(1+ROverZ_silicon_6_6_1->GetNbinsY()) = 2.4434610

nBinsROverZ = 42
minROverZ = 0.076
maxROverZ = 0.58

TCs = pd.read_csv('TCPositions/TCPositions_Zminus_siliconOnly.csv', sep=' ')
TCs["ROverZ"] = 1.0/np.sinh(abs(TCs["triggercelleta"]))
tower = {}
inclusive = ROOT.TH2D("inclusive","",nBinsROverZ,minROverZ,maxROverZ, nBinsPhi,minPhi,maxPhi)

for u in range(np.max(TCs['waferu'])): #wafer u
    for v in range(np.max(TCs['waferv'])): #wafer v
        for l in range(1, int(np.max(TCs['layer'])) ): #layer number
            tower[u,v,l] = ROOT.TH2D("tower_u"+str(u)+"_v"+str(v)+"_layer"+str(l),"",nBinsROverZ,minROverZ,maxROverZ, nBinsPhi,minPhi,maxPhi)
            for index, row in TCs[(TCs["waferu"]==u) & (TCs["waferv"]==v) & (TCs["layer"]==l)].iterrows():
                tower[u,v,l].Fill(row["ROverZ"], row["triggercellphi"])
            inclusive.Add(tower[u,v,l])

