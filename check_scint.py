import numpy as np
import pandas as pd
import ROOT

TCs = pd.read_csv('TCPositions/TCPositions_sctintillator.csv', sep=' ')
ietaVslayer = ROOT.TH2D("ietaVslayer","",14,37,51,22,0,22)

for index, row in TCs.iterrows():
    if(row['triggercelliphi']==1):
        ietaVslayer.Fill(row['layer'], row['triggercellieta']*-1)

c1 = ROOT.TCanvas()
ietaVslayer.GetXaxis().SetTitle('layer number')
ietaVslayer.GetYaxis().SetTitle('ieta')
ietaVslayer.SetTitle('ieta vs layer number for HGCal scintillator layers')
ietaVslayer.Draw("colz")
c1.SetGrid()
ietaVslayer.GetXaxis().SetTickLength(0)
ietaVslayer.GetYaxis().SetTickLength(0)
ietaVslayer.GetXaxis().SetNdivisions(14)
ietaVslayer.GetYaxis().SetNdivisions(22)
ROOT.gStyle.SetOptStat(0)

