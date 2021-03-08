import sys
import matplotlib.pyplot as plt
import pandas as pd
import math
import yaml
from funcs_plotter import find_eta, find_phi
#from funcs import getParMtxPerBundle, getModulesPerBundle, SaveHist #FIXME:The commented lines cause crash when imported with plt
#import ROOT

#ROOT.gStyle.SetOptStat(0)
#ROOT.gROOT.SetBatch()

def plotSCsOverTower(SC_file, outputdir):
    etaBinStep = 0.0870
    minBinEta = 24
    maxBinEta = 37
    
    phiBinStep = 2*math.pi/72
    minBinPhi = -4
    maxBinPhi = 8

    eta_ticks = [round(a * etaBinStep,3) for a in range(minBinEta, maxBinEta)] 
    phi_ticks = [round(b * phiBinStep,3) for b in range(minBinPhi, maxBinPhi)] 
    
    cells = pd.read_csv(SC_file, sep=' ')
    cells.columns= ["layer","waferu","waferv","triggercellu","triggercellv","SC_eta","SC_phi"]

    cells["SC_phi"] = cells["SC_phi"].replace(0, 1e-5)

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    
    layer = 1
    for u, v in [(2,0),(3,0),(3,1)]: # u and v: module coordinates
        cell_u_v_layer = cells[(cells['waferu'] == u) & (cells['waferv'] == v) & (cells['layer'] == layer)]
        plt.scatter(cell_u_v_layer['SC_eta']*-1.0, cell_u_v_layer['SC_phi'])
    
    ax.set_xticks(eta_ticks)
    ax.set_xticklabels(eta_ticks, rotation = 45)
    ax.set_yticks(phi_ticks)
    ax.grid(which='both')
    ax.set_xlabel('eta')
    ax.set_ylabel('phi')
    plt.show()
    fig.savefig(outputdir+'SCs_layer'+str(layer)+'.png', dpi=300)
    

def towersPerStage1(lpgbtMappingsFile, inputdir_paramMtx, param_mtx_em_name, param_mtx_had_name, outputdir):
    
    etaBinStep = 0.0870
    minBinEta = 16
    maxBinEta = 38
    minEta = minBinEta * etaBinStep
    maxEta = maxBinEta * etaBinStep
    nBinsEta = maxBinEta - minBinEta
    
    phiBinStep = 2*math.pi/72
    minBinPhi = -7
    maxBinPhi = 30
    minPhi = minBinPhi * phiBinStep
    maxPhi = maxBinPhi * phiBinStep
    nBinsPhi = maxBinPhi - minBinPhi
        
    inclusive = ROOT.TH2D("towerFit_inclusive","",nBinsEta,minEta,maxEta,nBinsPhi,minPhi,maxPhi)

    perStage1_towerFit_EM_below30deg = {} #would be like {0:2D Hist of towerFit of 1st Stage1 FPGA towers, 1:...2nd...,..., 23:...24th..}
    perStage1_towerFit_EM_above30deg = {} #would be like {0:2D Hist of towerFit of 1st Stage1 FPGA towers, 1:...2nd...,..., 23:...24th..}
    perStage1_towerFit_Had_below30deg = {}
    perStage1_towerFit_Had_above30deg = {}
    
    with open(lpgbtMappingsFile) as f:
        lines = [line.rstrip('\n') for line in f]
        f.close()
    bundles = getModulesPerBundle(lines)
    parMtxEM_PerBundle, parMtxHad_PerBundle = getParMtxPerBundle(bundles, inputdir_paramMtx, param_mtx_em_name, param_mtx_had_name)
    for bundle in parMtxEM_PerBundle: #bundle = 0,1,2,...,23
        sumPerTower = parMtxEM_PerBundle[bundle][(parMtxEM_PerBundle[bundle]!=0).any(axis=1)]
        sumPerTower = sumPerTower.astype(bool).astype(int).sum(axis=1)
        perStage1_towerFit_EM_below30deg[bundle] = ROOT.TH2D("towerFit_EM_below30deg_Stage1FPGA_num"+str(bundle),""\
                                                  ,nBinsEta,minEta,maxEta,nBinsPhi,minPhi,maxPhi)
        perStage1_towerFit_EM_above30deg[bundle] = ROOT.TH2D("towerFit_EM_above30deg_Stage1FPGA_num"+str(bundle),""\
                                                  ,nBinsEta,minEta,maxEta,nBinsPhi,minPhi,maxPhi)
        for tower in sumPerTower.index:
            if(find_phi(tower)<=5):
                perStage1_towerFit_EM_below30deg[bundle].SetBinContent(find_eta(tower)+2, find_phi(tower)+8, sumPerTower[tower])
                                                                                       #2 and 8 are just offset
            else:
                perStage1_towerFit_EM_above30deg[bundle].SetBinContent(find_eta(tower)+2, find_phi(tower)+8, sumPerTower[tower])
                                                                                       #2 and 8 are just offset
    
    for bundle in parMtxHad_PerBundle: #bundle = 0,1,2,...,23
        sumPerTower = parMtxHad_PerBundle[bundle][(parMtxHad_PerBundle[bundle]!=0).any(axis=1)]
        sumPerTower = sumPerTower.astype(bool).astype(int).sum(axis=1)
        perStage1_towerFit_Had_below30deg[bundle] = ROOT.TH2D("towerFit_Had_below30deg_Stage1FPGA_num"+str(bundle),""\
                                                   ,nBinsEta,minEta,maxEta,nBinsPhi,minPhi,maxPhi)
        perStage1_towerFit_Had_above30deg[bundle] = ROOT.TH2D("towerFit_Had_above30deg_Stage1FPGA_num"+str(bundle),""\
                                                   ,nBinsEta,minEta,maxEta,nBinsPhi,minPhi,maxPhi)
        for tower in sumPerTower.index:
            if(find_phi(tower)<=5):
                perStage1_towerFit_Had_below30deg[bundle].SetBinContent(find_eta(tower)+2, find_phi(tower)+8, sumPerTower[tower])
                                                                                       #2 and 8 are just offset
            else:
                perStage1_towerFit_Had_above30deg[bundle].SetBinContent(find_eta(tower)+2, find_phi(tower)+8, sumPerTower[tower])
                                                                                       #2 and 8 are just offset
    
    ####Printing
    fileType = 'pdf'
    for bundle in parMtxEM_PerBundle:
        perStage1_towerFit_EM_below30deg[bundle].GetXaxis().SetTitle('eta')
        perStage1_towerFit_EM_below30deg[bundle].GetYaxis().SetTitle('phi')
        SaveHist(perStage1_towerFit_EM_below30deg[bundle], outputdir, 'towerFit_EM_below30deg_Stage1FPGA_num'+str(bundle), fileType)
        
        perStage1_towerFit_EM_above30deg[bundle].GetXaxis().SetTitle('eta')
        perStage1_towerFit_EM_above30deg[bundle].GetYaxis().SetTitle('phi')
        SaveHist(perStage1_towerFit_EM_above30deg[bundle], outputdir, 'towerFit_EM_above30deg_Stage1FPGA_num'+str(bundle), fileType)
        
        perStage1_towerFit_Had_below30deg[bundle].GetXaxis().SetTitle('eta')
        perStage1_towerFit_Had_below30deg[bundle].GetYaxis().SetTitle('phi')
        SaveHist(perStage1_towerFit_Had_below30deg[bundle], outputdir, 'towerFit_Had_below30deg_Stage1FPGA_num'+str(bundle), fileType)
        
        perStage1_towerFit_Had_above30deg[bundle].GetXaxis().SetTitle('eta')
        perStage1_towerFit_Had_above30deg[bundle].GetYaxis().SetTitle('phi')
        SaveHist(perStage1_towerFit_Had_above30deg[bundle], outputdir, 'towerFit_Had_above30deg_Stage1FPGA_num'+str(bundle), fileType)

        inclusive.Add(perStage1_towerFit_EM_below30deg[bundle])
        inclusive.Add(perStage1_towerFit_EM_above30deg[bundle])
        inclusive.Add(perStage1_towerFit_Had_below30deg[bundle])
        inclusive.Add(perStage1_towerFit_Had_above30deg[bundle])
        inclusive.GetXaxis().SetTitle('eta')
        inclusive.GetYaxis().SetTitle('phi')
        
    SaveHist(inclusive, outputdir, 'inclusive', fileType)
    

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
    
    if (config['plotting']['towersPerStage1']):
        towersPerStage1(lpgbtMappingsFile=config['module_per_tower']['inputdir']+config['module_per_tower']['bundles_file'],\
                        inputdir_paramMtx=config['param_mtx']['outputdir'],\
                        param_mtx_em_name=config['param_mtx']['param_mtx_em_name'],\
                        param_mtx_had_name=config['param_mtx']['param_mtx_had_name'],\
                        outputdir=config['towersPerStage1']['outputdir']\
                        )
    
    if (config['plotting']['plotSCsOverTower']):
        plotSCsOverTower(SC_file=config['param_mtx']['inputdir']+config['param_mtx']['SC_position_file'],\
                        outputdir=config['plotSCsOverTower']['outputdir']\
                        )

if __name__ == "__main__":
    main()
