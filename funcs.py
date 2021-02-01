import itertools
import numpy as np

def partitions(n, k):
    for c in itertools.combinations(range(n+k-1), k-1):
        yield np.array([b-a-1 for a, b in zip((-1,)+c, c+(n+k-1,))])

def chisquare(vec1, vec2):
    chi2 = 0
    for i in range(len(vec1)):
        chi2 += np.power(vec1[i]-vec2[i], 2)
    return chi2

def calc_kernel(n):
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

    for bundle in bundles:
        for index, module in enumerate(bundles[bundle]):
            if(module[:5]=='type0'):
                bundles[bundle][index] = module[6:]
            else:
                print(20*'*'+'ERROR: Check the module name: ' + module + 20*'*')
                sys.exit(1)

    return bundles
