import itertools
import math
import ROOT
from  root_numpy import hist2array, array2hist
from scipy import ndimage
import sys
import numpy as np
import pandas as pd

def partitions(n, k): 
    """
    This is a generator function. It can be used for iteration. For example, if a=partitions(3, 2) then "a" can be used only once.
    
    partitions(n,k) yields ALL arrays of non-negative integers with len=k and sum=n
    
    e.g. patitions(3,2) are [0 3], [1 2], [2 1], [3 0]
    """
    for c in itertools.combinations(range(n+k-1), k-1):
        yield np.array([b-a-1 for a, b in zip((-1,)+c, c+(n+k-1,))])

def chisquare(vec1, vec2):
    if len(vec1)!=len(vec2):
        print(20*'*' + 'ERROR: Different vector length for Chi2-Squared' + 20*'*')
        sys.exit(1)
    chi2 = 0
    for i in range(len(vec1)):
        chi2 += np.power(vec1[i]-vec2[i], 2)
    return chi2

def calc_kernel(n): #n should be odd
    kernel = np.zeros((n,n))
    center = n//2
    kernel[center, center] = 1.
    for i in range(n):
        for j in range(n):
            if not (i==center and j==center):
                kernel[i,j] = 0.001 / ( (i-center)**2 + (j-center)**2 )
    
    kernel /= kernel.sum() 
    return kernel


def getModulesPerBundle(lines):

    bundles = {} #will be like {0:[], 1:[],..., 23:[]} with lists being module IDs
    bundleID = -1 #goes from 0 to 23
    for l in lines[1:]: # 1st line (column names) skipped
        if len(l)==0:
            print(20*'*'+'ERROR: Found an empty line in the file!'+20*'*')
            sys.exit(1)
        elif len(l.split())==1: #Lines with just the bundle number, i.e. 0, 1,2,...,23
            bundleID += 1
            bundles[bundleID]=[]
            continue
        else:
            bundles[bundleID] += ['type'+l.split()[5*i+2]+'-l'+l.split()[5*i+3]+'-u'+l.split()[5*i+4]+'-v'+l.split()[5*i+5] for i in range(int(l.split()[1]))]
    
    for bundle in bundles:
        bundles[bundle] = list(dict.fromkeys(bundles[bundle])) #remove repeated modules (some modules are connected to >1 lpgbt)
        for module in reversed(bundles[bundle]): #must be in reverse. Otherwise skips the next member when a member is removed (python trick)
            if(len(module))<5: #check for possible bugs
                print(20*'*'+'ERROR: Check the module name: ' + module + 20*'*')
                sys.exit(1)
            if(module[:5]=='type1'):
                bundles[bundle].remove(module) # remove scintillator modules as data is not yet included 
            elif(module[:5]!='type0'):
                print(20*'*'+'ERROR: Check the module name: ' + module + 20*'*')
                sys.exit(1)

    for bundle in bundles:#This loop could be merged with the previous one, but safer not to!
        for index, module in enumerate(bundles[bundle]):
            if(module[:5]=='type0'):
                bundles[bundle][index] = module[6:] #remove 'type0-'
            else:
                print(20*'*'+'ERROR: Check the module name: ' + module + 20*'*')
                sys.exit(1)

    return bundles

def getParMtxPerBundle(bundles, inputdir_paramMtx, param_mtx_em_name, param_mtx_had_name):
    parMtxEM = pd.read_pickle(inputdir_paramMtx + param_mtx_em_name).astype('int')
    parMtxHad = pd.read_pickle(inputdir_paramMtx + param_mtx_had_name).astype('int')

    
    parMtxEM_PerBundle = {}
    parMtxHad_PerBundle = {}

    for i in range(len(bundles)):
        parMtxEM_PerBundle[i] = parMtxEM[parMtxEM.columns.intersection(bundles[i])]
        parMtxHad_PerBundle[i] = parMtxHad[parMtxHad.columns.intersection(bundles[i])]
    
    return parMtxEM_PerBundle, parMtxHad_PerBundle

def writeParMtxPerBundleToFile(outputdir, parMtx, name):
    for i in parMtx:
        with open(outputdir + name + '_' + str(i) + '.txt', "w") as f:
            f.write('tower numModules moduleID fraction\n')
            for idx, tower in parMtx[i].iterrows():
                towerModuleOverlap = tuple(zip(tower.loc[tower!=0].index,\
                                        tower.loc[tower!=0].values))
                                        #becomes sth like (('l29-u0-v2', 7), ('l31-u0-v2', 3), ('l33-u0-v2', 1))
                f.write('{} {}'.format(tower.name, len(towerModuleOverlap)))
                for module in towerModuleOverlap:
                    f.write(' {} {}'.format(parMtx[i].columns.get_loc(module[0]), module[1]))
                f.write('\n')
        f.close()

def writeTowerPerModuleToFile(outputdir, parMtxEM, parMtxHad):
    with open(outputdir + 'tower_per_module.txt', 'w') as f:
        f.write('layer waferu waferv numTowers binEta binPhi fraction\n')
        for parMtx in [parMtxEM, parMtxHad]:
            for col in parMtx.columns:
                f.write('{} {} {} '.format(col[col.find('l')+1 : col.find('-u')],\
                                           col[col.find('u')+1 : col.find('-v')],\
                                           col[col.find('v')+1 : ]))
                towersInModule = parMtx[col].loc[parMtx[col]!=0]
                f.write('{} '.format(len(towersInModule)))
                for tower, frac in towersInModule.items(): #tower like 'had-eta2-phi23'. frac is 1,2,3,..
                   f.write('{} {} {} '.format(tower[tower.find('eta')+3 : tower.find('-phi')], tower[tower.find('-phi')+4:], frac))
                f.write('\n')
    f.close()


def getModulesWithTC(bundlesFileFullPath):
    with open(bundlesFileFullPath) as f:
        lines = [line.rstrip('\n') for line in f]
    f.close()
    bundles = getModulesPerBundle(lines)
    modules =[]
    for bundle in bundles:
        modules += bundles[bundle]
    return modules

def applyKernel(tower, kernel):
    tower_array = hist2array(tower)#convert to array
    towerSmoothed = ndimage.correlate(tower_array, kernel, mode='constant', cval = 0.0)
    for i in range(tower_array.shape[0]): #smoothing should not turn zeros to non-zeros
        for j in range(tower_array.shape[1]):
            if tower_array[i][j] == 0:
                towerSmoothed[i][j] = 0
    return towerSmoothed

def sortAndNormalize(tower, N_div):
    tower_flat = tower.flatten()
    tower_flat.sort()
    
    if (N_div<len(tower_flat)):
        if (tower_flat[-1 * N_div] == tower_flat[-1 * N_div -1] and tower_flat[-1 * N_div]!=0): #check degeneracy before removing low values
            print('check: degenerate!')
            print('cut value=', tower_flat[-1 * N_div], ', max=', tower_flat[-1], ' ratio=', tower_flat[-1 * N_div]/tower_flat[-1])
            print(20*'-')    
        tower_flat = tower_flat[-1 * N_div:] #keep N_div highest values
    
    tower_flat = tower_flat[tower_flat!=0] #remove zeros to prevent energy allocation to towers with zero overlap after fit.
    
    tower_flat = (tower_flat/tower_flat.sum())*N_div #normalize s.t. sum=N_div

    if len(tower_flat)==0:
        print(20*'*' + 'ERROR: Empty array! Should not be the case?: ' + 20*'*')
        sys.exit(1)
    return tower_flat

def findBestFit(towerSortedNormed, N_div):
    l = [] #check for degeneracy
    isDegenerate = False
    chi2_min = N_div * 1000000.0 #a very large number
    for p in partitions(N_div, len(towerSortedNormed)): #partitions(A,B) yields ALL arrays of non-negative integers with len=B and sum=A
        chi2_temp = chisquare(p,towerSortedNormed)
        if (chi2_temp <= chi2_min):
            chi2_min = chi2_temp
            bestFit = p 
            l.append(chi2_min)
    if(l.count(np.min(l))>1):#more than one best solution exist
        isDegenerate = True
        
    return bestFit, isDegenerate

def SaveHist(hist, outputdir, name, fileType):
    c = ROOT.TCanvas("canvas")
    hist.Draw('colz')
    nominal_low = ROOT.TLine(1.4,0,3.2,0)
    nominal_high = ROOT.TLine(1.4,2*math.pi/3,3.2,2*math.pi/3)
    nominal_low.Draw('same')
    nominal_low.SetLineStyle(7)
    nominal_high.Draw('same')
    nominal_high.SetLineStyle(7)
    c.SetRightMargin(0.15)
    c.Print(outputdir+name+'.'+fileType)
    #f = ROOT.TFile.Open(outputdir+name+'.'+whatType, 'RECREATE')
    #hist.Write()
    #f.Close()

def find_eta(name):
    return int(name[name.find('eta')+3 : name.find('-phi')])

def find_phi(name):
    return int(name[name.find('-phi')+4:])

def getPerStage1TowerHists(parMtx_PerBundle, coord, name):
    below30 = {}
    above30 = {}
    for bundle in parMtx_PerBundle: #bundle = 0,1,2,...,23
        sumPerTower = parMtx_PerBundle[bundle][(parMtx_PerBundle[bundle]!=0).any(axis=1)]
        sumPerTower = sumPerTower.astype(bool).astype(int).sum(axis=1)
        below30[bundle] = ROOT.TH2D("towerFit_"+name+"_below30deg_Stage1FPGA_num"+str(bundle),""\
                                    ,coord[0],coord[1],coord[2],coord[3],coord[4],coord[5])
        below30[bundle].SetTitle(";eta;phi;required number of sums")
        
        above30[bundle] = ROOT.TH2D("towerFit_"+name+"_above30deg_Stage1FPGA_num"+str(bundle),""\
                                    ,coord[0],coord[1],coord[2],coord[3],coord[4],coord[5])
        above30[bundle].SetTitle(";eta;phi;required number of sums")
        
        for tower in sumPerTower.index:
            if(find_phi(tower)<=5):
                below30[bundle].SetBinContent(find_eta(tower)+2, find_phi(tower)+8, sumPerTower[tower])#2 and 8 are just offset
            else:
                above30[bundle].SetBinContent(find_eta(tower)+2, find_phi(tower)+8, sumPerTower[tower])#2 and 8 are just offset
    return below30, above30



