import sys
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import math
import yaml
import matplotlib.patches as patches
from funcs import getParMtxPerBundle, getModulesPerBundle, SaveHist, getPerStage1TowerHists
import ROOT

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()

def plotSCsOverTower_singleModule(SC_file, outputdir): #for quick plotting single module position wrt tower coordinate
    layer = 11
    u = 1
    v = 3
    
    etaBinStep = 0.0870
    minBinEta = 16
    maxBinEta = 38
    
    phiBinStep = 2*math.pi/72
    minBinPhi = -2
    maxBinPhi = 24

    eta_ticks = [round(a * etaBinStep,3) for a in range(minBinEta, maxBinEta)] 
    phi_ticks = [round(b * phiBinStep,3) for b in range(minBinPhi, maxBinPhi)] 
    
    cells = pd.read_csv(SC_file, sep=' ')
    cells.columns= ["layer","waferu","waferv","triggercellu","triggercellv","SC_eta","SC_phi"]

    cells["SC_phi"] = cells["SC_phi"].replace(0, 1e-5)

    if not Path(outputdir).exists():
        Path(outputdir).mkdir(parents=True)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    
    cell_u_v_layer = cells[(cells['waferu'] == u) & (cells['waferv'] == v) & (cells['layer'] == layer)]
    plt.scatter(cell_u_v_layer['SC_eta']*-1.0, cell_u_v_layer['SC_phi'],s=10)

    ax.set_xticks(eta_ticks)
    ax.set_xticklabels(eta_ticks, rotation = 45)
    ax.set_yticks(phi_ticks)
    ax.grid(which='both')
    ax.set_xlabel('eta')
    ax.set_ylabel('phi')
    plt.tight_layout()
    fig.savefig(outputdir+'l'+str(layer)+'-u'+str(u)+'-v'+str(v)+'.png', dpi=300)

def plotSCsOverTower(SC_file, outputdir):
    etaBinStep = 0.0870
    minBinEta = 27
    maxBinEta = 39
    
    phiBinStep = 2*math.pi/72
    minBinPhi = -4
    maxBinPhi = 8

    etaBinToExamine = 31 #tower coordinate
    phiBinToExamine = 1 #tower coordinate

    eta_ticks = [round(a * etaBinStep,3) for a in range(minBinEta, maxBinEta)] 
    phi_ticks = [round(b * phiBinStep,3) for b in range(minBinPhi, maxBinPhi)] 
    
    cells = pd.read_csv(SC_file, sep=' ')
    cells.columns= ["layer","waferu","waferv","triggercellu","triggercellv","SC_eta","SC_phi"]

    cells["SC_phi"] = cells["SC_phi"].replace(0, 1e-5)

    outputdir +='/tower_EtaStep'+str(etaBinToExamine)+'_phiStep'+str(phiBinToExamine)+'/'
    if not Path(outputdir).exists():
        Path(outputdir).mkdir(parents=True)
    for layer in range(1,29):
        if ((layer%2==0) and (layer<=28)):
            continue
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        
        for u, v in [(2,0),(3,0),(3,1)]: # u and v: module coordinates
            cell_u_v_layer = cells[(cells['waferu'] == u) & (cells['waferv'] == v) & (cells['layer'] == layer)]
            plt.scatter(cell_u_v_layer['SC_eta']*-1.0, cell_u_v_layer['SC_phi'],s=10)

        rect = patches.Rectangle((etaBinStep*31, phiBinStep*1), etaBinStep, phiBinStep, linewidth=2, edgecolor='black', facecolor='none')
        ax.add_patch(rect)
        
        ax.set_xticks(eta_ticks)
        ax.set_xticklabels(eta_ticks, rotation = 45)
        ax.set_yticks(phi_ticks)
        ax.grid(which='both')
        ax.set_xlabel('eta')
        ax.set_ylabel('phi')
        ax.set_axisbelow(True)
        plt.tight_layout()
        fig.savefig(outputdir+'SCs_layer'+str(layer)+'.png', dpi=300)
        plt.clf()

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
    inclusive.SetTitle(";eta;phi;required number of sums")
    
    with open(lpgbtMappingsFile) as f:
        lines = [line.rstrip('\n') for line in f]
        f.close()
    bundles = getModulesPerBundle(lines, isScintil=False)
    parMtxEM_PerBundle, parMtxHad_PerBundle = getParMtxPerBundle(bundles, inputdir_paramMtx, param_mtx_em_name, param_mtx_had_name)
   
    coord = [nBinsEta,minEta,maxEta,nBinsPhi,minPhi,maxPhi]
    #below are dictionaries in this format: {0:2D Hist of towerFit of 1st Stage1 FPGA towers, 1:...2nd...,..., 23rd:...}
    perStage1_towerFit_EM_below30deg, perStage1_towerFit_EM_above30deg = getPerStage1TowerHists(parMtxEM_PerBundle, coord, name='EM')
    perStage1_towerFit_Had_below30deg, perStage1_towerFit_Had_above30deg = getPerStage1TowerHists(parMtxHad_PerBundle, coord, name='Had')
    
    ####Printing
    fileType = 'pdf'
    if not Path(outputdir).exists():
        Path(outputdir).mkdir(parents=True)
    for bundle in parMtxEM_PerBundle:
        SaveHist(perStage1_towerFit_EM_below30deg[bundle], outputdir, 'towerFit_EM_below30deg_Stage1FPGA_num'+str(bundle), fileType)
        SaveHist(perStage1_towerFit_EM_above30deg[bundle], outputdir, 'towerFit_EM_above30deg_Stage1FPGA_num'+str(bundle), fileType)
        SaveHist(perStage1_towerFit_Had_below30deg[bundle], outputdir, 'towerFit_Had_below30deg_Stage1FPGA_num'+str(bundle), fileType)
        SaveHist(perStage1_towerFit_Had_above30deg[bundle], outputdir, 'towerFit_Had_above30deg_Stage1FPGA_num'+str(bundle), fileType)

        inclusive.Add(perStage1_towerFit_EM_below30deg[bundle])
        inclusive.Add(perStage1_towerFit_EM_above30deg[bundle])
        inclusive.Add(perStage1_towerFit_Had_below30deg[bundle])
        inclusive.Add(perStage1_towerFit_Had_above30deg[bundle])
        
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
    
    if (config['plotterFuncs']['towersPerStage1']):
        towersPerStage1(lpgbtMappingsFile=config['module_per_tower']['inputdir']+config['module_per_tower']['bundles_file'],\
                        inputdir_paramMtx=config['param_mtx']['outputdir'],\
                        param_mtx_em_name=config['param_mtx']['param_mtx_em_name'],\
                        param_mtx_had_name=config['param_mtx']['param_mtx_had_name'],\
                        outputdir=config['towersPerStage1']['outputdir']\
                        )
    
    if (config['plotterFuncs']['plotSCsOverTower']):
        plotSCsOverTower(SC_file=config['param_mtx']['inputdir']+config['param_mtx']['SC_position_file'],\
                        outputdir=config['plotSCsOverTower']['outputdir']\
                        )

    if (config['plotterFuncs']['plotSCsOverTower_singleModule']):
        plotSCsOverTower_singleModule(SC_file=config['param_mtx']['inputdir']+config['param_mtx']['SC_position_file'],\
                        outputdir=config['plotSCsOverTower_singleModule']['outputdir']\
                        )

if __name__ == "__main__":
    main()
