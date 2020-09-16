import pandas as pd
import numpy as np
import ROOT
hist = ROOT.TH2D("hist","",300,1.48,3.23, 12,0,12)
TCs = pd.read_csv('TCPositions/sensorCell_positions.txt', sep=' ')
TCs.columns= ["layer","b","c","d","e","SC_eta","SC_phi"]
#TCs = TCs[abs(TCs['eta'])>3]
#sample_size = 12
#hist.Fill(len(TCs[(TCs['layer']==1)&(TCs['b']==3)&(TCs['c']==3)&(TCs['d']==3)&(TCs['e']==3)]))

#for m in np.random.choice(1+int(np.max(TCs['layer'])), 25, replace = False):
#    print(m)
#    for i in np.random.choice(1+np.max(TCs['b']), 5 ,  replace = False):
#        for j in np.random.choice(1+np.max(TCs['c']), 5, replace = False):
#            for k in np.random.choice( 1+np.max(TCs['d']), 3, replace = False):
#                for l in np.random.choice( 1+np.max(TCs['e']), 3, replace = False):

for m in range(1, (1+int(np.max(TCs['layer'])))):
    print(m)
    for i in range(0, (1+np.max(TCs['b']))):
        for j in range(0, (1+np.max(TCs['c']))):
            for k in range(0, (1+np.max(TCs['d']))):
                for l in range(0, (1+np.max(TCs['e']))):
                    num_sCell = len(TCs[(TCs['layer']==m)&(TCs['b']==i)&(TCs['c']==j)&(TCs['d']==k)&(TCs['e']==l)])
                    eta_TC = -1.0*np.mean(TCs['SC_eta'][(TCs['layer']==m)&(TCs['b']==i)&(TCs['c']==j)&(TCs['d']==k)&(TCs['e']==l)])
                    #print(num_sCell, m,i,j,k,l)
                    if num_sCell!=0:
                        hist.Fill(eta_TC, num_sCell)
                    if num_sCell!=9 and num_sCell!=4 and eta_TC<3.0 and eta_TC>1.6:
                        print("num_sCell=",num_sCell,"eta_TC=",eta_TC,"layer=",m, " b=",i, " c=",j, " d=",k, " e=",l)
                        if num_sCell==8:
                            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

f1 = ROOT.TFile("sensor_cells.root","RECREATE")
hist.Write()       
f1.Close()
