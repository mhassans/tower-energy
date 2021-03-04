import sys
import time
import math
import yaml
from funcs import getParMtxPerBundle, getModulesPerBundle, SaveHist
from funcs_plotter import find_eta, find_phi
import ROOT

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()

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

    perStage1_towerFit_EM = {} #would be like {0:2D Hist of towerFit of 1st Stage1 FPGA towers, 1:...2nd...,..., 23:...24th..}
    perStage1_towerFit_Had = {}
    
    with open(lpgbtMappingsFile) as f:
        lines = [line.rstrip('\n') for line in f]
        f.close()
    bundles = getModulesPerBundle(lines)
    parMtxEM_PerBundle, parMtxHad_PerBundle = getParMtxPerBundle(bundles, inputdir_paramMtx, param_mtx_em_name, param_mtx_had_name)
    for bundle in parMtxEM_PerBundle: #bundle = 0,1,2,...,23
        sumPerTower = parMtxEM_PerBundle[bundle][(parMtxEM_PerBundle[bundle]!=0).any(axis=1)].sum(axis=1)
        perStage1_towerFit_EM[bundle] = ROOT.TH2D("towerFit_EM_Stage1FPGA_num"+str(bundle),""\
                                                  ,nBinsEta,minEta,maxEta,nBinsPhi,minPhi,maxPhi)
        for tower in sumPerTower.index:
            perStage1_towerFit_EM[bundle].SetBinContent(find_eta(tower)+2, find_phi(tower)+8, sumPerTower[tower])#2 and 8 are just offset
    for bundle in parMtxHad_PerBundle: #bundle = 0,1,2,...,23
        sumPerTower = parMtxHad_PerBundle[bundle][(parMtxHad_PerBundle[bundle]!=0).any(axis=1)].sum(axis=1)
        perStage1_towerFit_Had[bundle] = ROOT.TH2D("towerFit_Had_Stage1FPGA_num"+str(bundle),""\
                                                   ,nBinsEta,minEta,maxEta,nBinsPhi,minPhi,maxPhi)
        for tower in sumPerTower.index:
            perStage1_towerFit_Had[bundle].SetBinContent(find_eta(tower)+2, find_phi(tower)+8, sumPerTower[tower])#2 and 8 are just offset
    for bundle in parMtxEM_PerBundle:
        perStage1_towerFit_EM[bundle].GetXaxis().SetTitle('eta')
        perStage1_towerFit_EM[bundle].GetYaxis().SetTitle('phi')
        SaveHist(perStage1_towerFit_EM[bundle], outputdir, 'towerFit_EM_Stage1FPGA_num'+str(bundle), 'root')
        perStage1_towerFit_Had[bundle].GetXaxis().SetTitle('eta')
        perStage1_towerFit_Had[bundle].GetYaxis().SetTitle('phi')
        SaveHist(perStage1_towerFit_Had[bundle], outputdir, 'towerFit_Had_Stage1FPGA_num'+str(bundle), 'root')
    

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

if __name__ == "__main__":
    start = time.time()
    main()
    print('The program ran in', time.time() - start, 'seconds!')
