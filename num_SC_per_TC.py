import pandas as pd
import numpy as np
import ROOT
hist = ROOT.TH2D("hist","",300,1.48,3.23, 12,0,12)
SCs = pd.read_csv('TCPositions/sensorCell_positions.txt', sep=' ')
SCs.columns= ["layer","waferu","waferv","triggercellu","triggercellv","SC_eta","SC_phi"]
#SCs = SCs[abs(SCs['eta'])>3]
#sample_size = 12
#hist.Fill(len(SCs[(SCs['layer']==1)&(SCs['waferu']==3)&(SCs['waferv']==3)&(SCs['triggercellu']==3)&(SCs['triggercellv']==3)]))

#for m in np.random.choice(1+int(np.max(SCs['layer'])), 25, replace = False):
#    print(m)
#    for i in np.random.choice(1+np.max(SCs['waferu']), 5 ,  replace = False):
#        for j in np.random.choice(1+np.max(SCs['waferv']), 5, replace = False):
#            for k in np.random.choice( 1+np.max(SCs['triggercellu']), 3, replace = False):
#                for l in np.random.choice( 1+np.max(SCs['triggercellv']), 3, replace = False):

for m in range(1, (1+int(np.max(SCs['layer'])))):
    print(m)
    for i in range(0, (1+np.max(SCs['waferu']))):
        for j in range(0, (1+np.max(SCs['waferv']))):
            for k in range(0, (1+np.max(SCs['triggercellu']))):
                for l in range(0, (1+np.max(SCs['triggercellv']))):
                    num_sCell = len(SCs[(SCs['layer']==m)&(SCs['waferu']==i)&(SCs['waferv']==j)&(SCs['triggercellu']==k)&(SCs['triggercellv']==l)])
                    eta_TC = -1.0*np.mean(SCs['SC_eta'][(SCs['layer']==m)&(SCs['waferu']==i)&(SCs['waferv']==j)&(SCs['triggercellu']==k)&(SCs['triggercellv']==l)])
                    #print(num_sCell, m,i,j,k,l)
                    if num_sCell!=0:
                        hist.Fill(eta_TC, num_sCell)
                    if num_sCell!=9 and num_sCell!=4 and eta_TC<3.0 and eta_TC>1.6:
                        print("num_sCell=",num_sCell,"eta_TC=",eta_TC,"layer=",m, " waferu=",i, " waferv=",j, " triggercellu=",k, " triggercellv=",l)
                        if num_sCell==8:
                            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

f1 = ROOT.TFile("sensor_cells.root","RECREATE")
hist.Write()       
f1.Close()
