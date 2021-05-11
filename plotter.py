import sys
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import math
import yaml
import matplotlib.patches as patches
from funcs import getParMtxPerBundle_Silic, getParMtxPerBundle_Scint, getModulesPerBundle, SaveHist, getPerStage1TowerHists, fillTowersIncl
import ROOT

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()

def plotSCsOverTower_singleModule(SC_file, outputdir): #for quick plotting single module position wrt tower coordinate
    layer = 11
    u = 1
    v = 3
    
    etaBinStep = 0.0870
    minBinEta = 16 #to be tuned manually
    maxBinEta = 38 #to be tuned manually
    
    phiBinStep = 2*math.pi/72
    minBinPhi = -2 #to be tuned manually
    maxBinPhi = 24 #to be tuned manually

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
    minBinEta = 27 #to be tuned manually
    maxBinEta = 39 #to be tuned manually
    
    phiBinStep = 2*math.pi/72
    minBinPhi = -4 #to be tuned manually
    maxBinPhi = 8 #to be tuned manually

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

def towersPerStage1(lpgbtMappingsFile, inputdir_paramMtx, param_mtx_em_name, param_mtx_had_name_silic, param_mtx_had_name_scint, outputdir, histRange):
    
    etaBinStep = histRange['etaBinStep']
    minBinEta = histRange['minBinEta']
    maxBinEta = histRange['maxBinEta']
    minEta = minBinEta * etaBinStep
    maxEta = maxBinEta * etaBinStep
    nBinsEta = maxBinEta - minBinEta
    
    phiBinStep = eval(histRange['phiBinStep'])
    minBinPhi = histRange['minBinPhi']
    maxBinPhi = histRange['maxBinPhi']
    minPhi = minBinPhi * phiBinStep
    maxPhi = maxBinPhi * phiBinStep
    nBinsPhi = maxBinPhi - minBinPhi
        
    inclusive = ROOT.TH2D("towerFit_inclusive","",nBinsEta,minEta,maxEta,nBinsPhi,minPhi,maxPhi)
    inclusive.SetTitle(";eta;phi;required number of sums")
    
    with open(lpgbtMappingsFile) as f:
        lines = [line.rstrip('\n') for line in f]
        f.close()
    bundlesSilic = getModulesPerBundle(lines, isScintil=False)
    bundlesScint = getModulesPerBundle(lines, isScintil=True)
    
    parMtxEM_PerBundle, parMtxHadSilic_PerBundle = getParMtxPerBundle_Silic(bundlesSilic, inputdir_paramMtx,\
                                                                        param_mtx_em_name, param_mtx_had_name_silic)
    parMtxHadScint_PerBundle = getParMtxPerBundle_Scint(bundlesScint, inputdir_paramMtx, param_mtx_had_name_scint)

    #merging silicon and scintillators:
    parMtxHad_PerBundle = {}
    for i in range(len(bundlesScint)):
        parMtxHad_PerBundle[i] = pd.concat([parMtxHadSilic_PerBundle[i],parMtxHadScint_PerBundle[i]], axis=1).fillna(0)
    
    coord = [nBinsEta,minEta,maxEta,nBinsPhi,minPhi,maxPhi]
    #below are dictionaries in this format: {0:2D Hist of towerFit of 1st Stage1 FPGA towers, 1:...2nd...,..., 23rd:...}
    perStage1_towerFit_EM_below30deg, perStage1_towerFit_EM_above30deg = getPerStage1TowerHists(parMtxEM_PerBundle, coord, name='EM')
    perStage1_towerFit_Had_below30deg, perStage1_towerFit_Had_above30deg = getPerStage1TowerHists(parMtxHad_PerBundle, coord, name='Had')
    
    ####Printing
    fileType = 'pdf'
    if not Path(outputdir).exists():
        Path(outputdir).mkdir(parents=True)
    for bundle in parMtxEM_PerBundle:
        SaveHist(perStage1_towerFit_EM_below30deg[bundle], outputdir, 'towerFit_EM_below30deg_Stage1FPGA_num'+str(bundle),\
                        fileType, AddGrid=True)
        SaveHist(perStage1_towerFit_EM_above30deg[bundle], outputdir, 'towerFit_EM_above30deg_Stage1FPGA_num'+str(bundle),\
                        fileType, AddGrid=True)
        SaveHist(perStage1_towerFit_Had_below30deg[bundle], outputdir, 'towerFit_Had_below30deg_Stage1FPGA_num'+str(bundle),\
                        fileType, AddGrid=True)
        SaveHist(perStage1_towerFit_Had_above30deg[bundle], outputdir, 'towerFit_Had_above30deg_Stage1FPGA_num'+str(bundle),\
                        fileType, AddGrid=True)

        inclusive.Add(perStage1_towerFit_EM_below30deg[bundle])
        inclusive.Add(perStage1_towerFit_EM_above30deg[bundle])
        inclusive.Add(perStage1_towerFit_Had_below30deg[bundle])
        inclusive.Add(perStage1_towerFit_Had_above30deg[bundle])
        
    SaveHist(inclusive, outputdir, 'inclusive', fileType, AddGrid=True)

def TotNumOfSumsPerTower(outputdir, histRange):

    etaBinStep = histRange['etaBinStep']
    minBinEta = histRange['minBinEta']
    maxBinEta = histRange['maxBinEta']
    minEta = minBinEta * etaBinStep
    maxEta = maxBinEta * etaBinStep
    nBinsEta = maxBinEta - minBinEta
    
    phiBinStep = eval(histRange['phiBinStep'])
    minBinPhi = histRange['minBinPhi']
    maxBinPhi = histRange['maxBinPhi']
    minPhi = minBinPhi * phiBinStep
    maxPhi = maxBinPhi * phiBinStep
    nBinsPhi = maxBinPhi - minBinPhi

    parMtxHadScint = pd.read_pickle('./output/tower_module_mapping_array/param_mtx_had_scint.pkl')
    parMtxHadSilic = pd.read_pickle('./output/tower_module_mapping_array/param_mtx_had_silic.pkl')
    parMtxEM = pd.read_pickle('./output/tower_module_mapping_array/param_mtx_em.pkl')
    
    HadSilic = ROOT.TH2D("HadSilic","",nBinsEta,minEta,maxEta,nBinsPhi,minPhi,maxPhi)
    HadScint = ROOT.TH2D("HadScint","",nBinsEta,minEta,maxEta,nBinsPhi,minPhi,maxPhi)
    Had = ROOT.TH2D("Had","",nBinsEta,minEta,maxEta,nBinsPhi,minPhi,maxPhi) #=HadScint+HadSilic
    Had.SetTitle("Distribution of total number of sums per tower;eta;phi;required number of sums")
    EM = ROOT.TH2D("EM","",nBinsEta,minEta,maxEta,nBinsPhi,minPhi,maxPhi)
    inclusive = ROOT.TH2D("Inclusive","",nBinsEta,minEta,maxEta,nBinsPhi,minPhi,maxPhi) #=Had+EM
    inclusive.SetTitle("Distribution of total number of sums per tower;eta;phi;required number of sums")

    HadSilic = fillTowersIncl(parMtxHadSilic, HadSilic, DoNumOfSums=True, N_div=8) #N_div is ignored when DoNumOfSums=True
    HadScint = fillTowersIncl(parMtxHadScint, HadScint, DoNumOfSums=True, N_div=16)
    Had.Add(HadScint)
    Had.Add(HadSilic)
    EM = fillTowersIncl(parMtxEM, EM, DoNumOfSums=True, N_div=8)
    inclusive.Add(Had)
    inclusive.Add(EM)

    ####Printing
    fileType = 'pdf'
    if not Path(outputdir).exists():
        Path(outputdir).mkdir(parents=True)
        
    SaveHist(HadSilic, outputdir, 'HadSilic', fileType, AddGrid=False)
    SaveHist(HadScint, outputdir, 'HadScint', fileType, AddGrid=False)
    SaveHist(Had, outputdir, 'Had', fileType, AddGrid=False)
    SaveHist(EM, outputdir, 'EM', fileType, AddGrid=False)
    SaveHist(inclusive, outputdir, 'inclusive', fileType, AddGrid=False)

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
                        param_mtx_had_name_silic=config['param_mtx']['param_mtx_had_name_silic'],\
                        param_mtx_had_name_scint=config['param_mtx']['param_mtx_had_name_scint'],\
                        outputdir=config['towersPerStage1']['outputdir'],\
                        histRange=config['histRange']\
                        )
    
    if (config['plotterFuncs']['plotSCsOverTower']):
        plotSCsOverTower(SC_file=config['param_mtx']['inputdir']+config['param_mtx']['SC_position_file'],\
                        outputdir=config['plotSCsOverTower']['outputdir']\
                        )

    if (config['plotterFuncs']['plotSCsOverTower_singleModule']):
        plotSCsOverTower_singleModule(SC_file=config['param_mtx']['inputdir']+config['param_mtx']['SC_position_file'],\
                        outputdir=config['plotSCsOverTower_singleModule']['outputdir']\
                        )

    if (config['plotterFuncs']['TotNumOfSumsPerTower']):
        TotNumOfSumsPerTower(outputdir=config['TotNumOfSumsPerTower']['outputdir'],\
                             histRange=config['histRange']\
                            )
        

if __name__ == "__main__":
    main()
