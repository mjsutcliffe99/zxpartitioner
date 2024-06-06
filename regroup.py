import pyzx as zx
from .partition import get_unique_params

def decomp(g):
    gDecomp = zx.simulate.find_stabilizer_decomp(g)
    #print("No. of terms:",len(gDecomp))
    total = (0+0j)
    for i in gDecomp: total += i.scalar.to_number()
    return total

def prepParams(g):
    pVerts = [] # parameterised vertices
    pVertsParams = [] # the parameters of the parameterised vertices
    localParams = set() # the SET of unique local parameters
    for v in g.vertices():
        if len(g.get_params(v)) > 0:
            pVerts.append(v)
            p = list(g.get_params(v))[0] # TEMP - currently only works for verts with one param
            pVertsParams.append(p)
            localParams.add(p)
            g.set_params(v,set())
    
    localParams = list(localParams)
    localParams.sort()
    return g,pVerts,pVertsParams,localParams

def intToBin(x,n_params): return format(x, '#0'+str(n_params+2)+'b')[2:] # TODO - improve this!

# TODO - this is very messy and inefficient right now...
def iterateCuts(g):
    g = g.copy()
    g,pVerts,pVertsParams,localParams = prepParams(g)
    gOrig = g.copy()
    n_loc_params = len(localParams)
    scalars = [0]*(2**n_loc_params)
    for i in range(2**n_loc_params):
        g = gOrig.copy()
        strBin = intToBin(i,n_loc_params)
        for j,v in enumerate(pVerts):
            p_index = localParams.index(pVertsParams[j]) # e.g. 'a' -> 0
            phase = int(strBin[p_index]) # e.g. if a=0 let phase=0, if a=1 let phase=1
            g.set_phase(v,phase)
        for pp in g.scalar.phasepairs:
            pp_B_index = localParams.index(list(pp.paramsB)[0])
            B = int(strBin[pp_B_index]) # B = 0 or 1
            if B == 1: g.scalar.add_phase(pp.alpha/4) # TODO - should really also remove the param phasepair data here (but it shouldn't matter)
        zx.simplify.full_reduce(g)
        scalars[i] = decomp(g)
    return localParams,scalars

def globalToLocalBits(globalParams,localParams,globBitstr):
    locBitstr = ''
    for i,p in enumerate(globalParams):
        if p in localParams: locBitstr += globBitstr[i]
    return locBitstr

class Segment: # TODO - should be in HVert
    localParams = []
    scalars = []
    
    def __init__(self,g=None):
        if not g == None:
            localParams,scalars = iterateCuts(g)
            self.localParams = localParams
            self.scalars = scalars
            return
        self.localParams = []
        self.scalars = []

def getExclusivelyCommonParams(segs,A,B):
    commonParams = set(segs[A].localParams).intersection(set(segs[B].localParams))
    nonexclusives = set()
    for p in commonParams:
        for i,seg in enumerate(segs):
            if i in {A,B}: continue
            if p in seg.localParams: nonexclusives.add(p)
    return commonParams.symmetric_difference(nonexclusives) # Return params EXCLUSIVELY common to A and B

def regroupPair(segs,A,B):
    segA = segs[A]
    segB = segs[B]
    
    pairParams = list(set(segA.localParams).union(set(segB.localParams))) # all params common to A and B
    pairParams.sort()
    n_pair_params = len(pairParams)
    
    noncommonParams = list(set(pairParams).symmetric_difference(getExclusivelyCommonParams(segs,A,B)))
    noncommonParams.sort()
    newScalars = [0]*(2**len(noncommonParams))
    
    for i in range(2**n_pair_params):
        bitstr = intToBin(i,n_pair_params)
        
        ab = globalToLocalBits(pairParams,segA.localParams,bitstr)
        if ab == '': ab = '0' # Special case (free nodes)
        A_ab = segA.scalars[int(ab,2)]
        
        bc = globalToLocalBits(pairParams,segB.localParams,bitstr)
        if bc == '': bc = '0' # Special case (free nodes)
        B_bc = segB.scalars[int(bc,2)]
    
        AB_abc = A_ab*B_bc
        ac = globalToLocalBits(pairParams,noncommonParams,bitstr)
        if ac == '': ac = '0' # If all segments have combined and no params remain
        newScalars[int(ac,2)] += AB_abc
    
    segAB = Segment()
    segAB.localParams = noncommonParams
    segAB.scalars = newScalars

    segs[A] = segAB
    segs[B] = Segment()
    
    return True
    
#----------------
    
# Estimate runtime for precomp...
def estimateCostPrecomp(gs,k):
    calcs_tot = 0
    for i in range(k):
        calcs_tot += 2**(len(get_unique_params(gs[i])))
    print(calcs_tot)
    #print(" ~",max([calcs_A,calcs_B,calcs_C,calcs_D]))
    return calcs_tot

# PRECOMPILE SEGMENTS...
def precompSegments(gs,k):
    segs = [None]*k
    for i in range(k): segs[i] = Segment(gs[i])
    return segs

# Estimate runtime for cross-ref...
def estimateCostCrossref(hNet):
    hNetCopy = hNet.copy()
    calcs_tot = 0
    steps = 0
    while len(hNetCopy.getActiveVerts()) > 1:
        fusePair,cost = hNetCopy.stepReduce(False)
        calcs_tot += 2**cost
        steps += 1
    
    hNetCopy = hNet.copy()
    max_crossref_calcs = hNetCopy.fullReduce(True,False)
    print("\n",calcs_tot)
    print("  ~",2**max_crossref_calcs)
    return calcs_tot#,steps
    
# CROSS-REFERENCE SEGMENTS...
def crossrefSegments(segs,hNet):
    #calcs,steps = estimateCostCrossref() # TEMP
    steps = len(segs)-1

    hNetCopy = hNet.copy()
    for i in range(steps):
        fusePair,cost = hNetCopy.stepReduce(False)
        X = fusePair[0]
        Y = fusePair[1]
        #print(X,Y) #TEMP
        regroupPair(segs,X,Y)
    
    #assert len(segs[0].scalars) == 1
    result = segs[0].scalars[0]
    return result