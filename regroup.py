# MATTHEW SUTCLIFFE, 2024

import pyzx as zx
from .utils import *
import cuda
import cupy
import math

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
    newScalars = [(0+0j)]*(2**len(noncommonParams))
    
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
    
def regroupPairGPU(segs,A,B):
    segA = segs[A]
    segB = segs[B]
    
    pairParams = list(set(segA.localParams).union(set(segB.localParams))) # all params common to A and B
    pairParams.sort()
    n_pair_params = len(pairParams)
    size = 2**n_pair_params
    
    noncommonParams = list(set(pairParams).symmetric_difference(getExclusivelyCommonParams(segs,A,B)))
    noncommonParams.sort()
    n_ex_params = len(noncommonParams) # exclusive param pairs (i.e. those of the future grouped AB segment)
    
    #--
    
    #                  e.g.     a b c
    A_params     = 0   #  6  =  1 1 0
    B_params     = 0   #  3  =  0 1 1
    AB_ex_params = 0   #  5  =  1 0 1   (be careful to not collapse (n>2)-edged hyderedges
    
    for i,p in enumerate(pairParams):
        if p in segA.localParams: A_params += 2**(n_pair_params-i-1)
        if p in segB.localParams: B_params += 2**(n_pair_params-i-1)
        if p in noncommonParams:  AB_ex_params += 2**(n_pair_params-i-1)
    
    A_re  = cupy.zeros(2**len(segA.localParams), dtype=cupy.float32)
    A_im  = cupy.zeros(2**len(segA.localParams), dtype=cupy.float32)
    B_re  = cupy.zeros(2**len(segB.localParams), dtype=cupy.float32)
    B_im  = cupy.zeros(2**len(segB.localParams), dtype=cupy.float32)
    AB_re = cupy.zeros(2**n_ex_params, dtype=cupy.float32)
    AB_im = cupy.zeros(2**n_ex_params, dtype=cupy.float32)
    
    # separate the complex scalars into real and imag components...
    for i,s in enumerate(segA.scalars):
        A_re[i] = s.real
        A_im[i] = s.imag
    for i,s in enumerate(segB.scalars):
        B_re[i] = s.real
        B_im[i] = s.imag
    
    # CUDA code for regrouping pair...
    cuda_code = r'''
    extern "C"
    __global__ void regroup_pair_gpu(int paramsA, int paramsB, int paramsC, float * A_re, float * A_im, float * B_re, float * B_im, float * AB_re, float * AB_im, const int N_params, const int size)
    {
        int index = blockIdx.x * blockDim.x + threadIdx.x;
        
        // LOCALLY INDEX...
        
        int ab = 0;
        int bc = 0;
        int ac = 0;
        int abc = index;
        int x = 0; // current length of ab
        int y = 0; // current length of bc
        int z = 0; // current length of ac
        
        for (int i=0; i<N_params; ++i)
        {
            if (paramsA & 1) ab = ((abc & 1) << x++) | ab;
            if (paramsB & 1) bc = ((abc & 1) << y++) | bc;
            if (paramsC & 1) ac = ((abc & 1) << z++) | ac;
            abc >>= 1;
            paramsA >>= 1;
            paramsB >>= 1;
            paramsC >>= 1;
        }
        
        // MULTIPLY SCALARS (A_ab * B_bc -> AB_abc) ...
        // (A+ai)(B+bi) = (AB-ab) + (Ab+aB)i
        
        float A = A_re[ab];
        float a = A_im[ab];
        float B = B_re[bc];
        float b = B_im[bc];

        atomicAdd(&AB_re[ac], (A*B) - (a*b));  //AB_re[index] = (A*B) - (a*b);
        atomicAdd(&AB_im[ac], (A*b) + (a*B));  //AB_im[index] = (A*b) + (a*B);
        __syncthreads();
    }
    '''
    regroup_pair_gpu = cupy.RawKernel(cuda_code, "regroup_pair_gpu")
    
    threads_per_block = 1024
    grid_size  = (1,1,1)
    block_size = (size,1,1)
    if size > 1024:
        grid_size  = (math.ceil(size/threads_per_block),1,1)
        block_size = (threads_per_block,1,1)
    
    #%%time
    regroup_pair_gpu(grid_size, block_size, (A_params, B_params, AB_ex_params, A_re, A_im, B_re, B_im, AB_re, AB_im, n_pair_params, size)) #regroup_pair_gpu(grid_size, block_size, (A_params, B_params, AB_ex_params, A_re, A_im, B_re, B_im, AB_re, AB_im, n_pair_params, size))
    cupy.cuda.stream.get_current_stream().synchronize()
    
    # recollect the real and imag components into a scalar...
    newScalars = AB_re + (AB_im*1j)
    if type(newScalars) == cupy.ndarray: newScalars = newScalars.tolist()
    
    #--
    
    segAB = Segment()
    segAB.localParams = noncommonParams
    segAB.scalars = newScalars

    segs[A] = segAB
    segs[B] = Segment()
    
    return True
    
#----------------
    
# Estimate runtime for precomp...
def estimateCostPrecomp(gs,doPrint=False):
    calcs_tot = 0
    for g in gs:
        calcs_tot += 2**(len(get_unique_params(g)) + (EST_ALPHA*zx.tcount(g)))
    if doPrint: print(calcs_tot)
    #print(" ~",max([calcs_A,calcs_B,calcs_C,calcs_D]))
    return calcs_tot

# PRECOMPILE SEGMENTS...
def precompSegments(gs):
    k = len(gs)
    segs = [None]*k
    for i in range(k): segs[i] = Segment(gs[i])
    return segs

# Estimate runtime for cross-ref...
def estimateCostCrossref(hNet,doPrint=False):
    hNetCopy = hNet.copy()
    calcs_tot = 0
    steps = 0
    while len(hNetCopy.getActiveVerts()) > 1:
        fusePair,cost = hNetCopy.stepReduce(False)
        calcs_tot += 2**cost
        steps += 1
    
    hNetCopy = hNet.copy()
    max_crossref_calcs = hNetCopy.fullReduce(doPrint,False)
    if doPrint:
        print("\n",calcs_tot)
        print("  ~",2**max_crossref_calcs)
    return calcs_tot#,steps
    
# CROSS-REFERENCE SEGMENTS...
def crossrefSegments(segs,hNet,useGPU=False):
    #calcs,steps = estimateCostCrossref() # TEMP
    steps = len(segs)-1

    hNetCopy = hNet.copy()
    for i in range(steps):
        fusePair,cost = hNetCopy.stepReduce(False)
        X = fusePair[0]
        Y = fusePair[1]
        #print(X,Y) #TEMP
        if useGPU: regroupPairGPU(segs,X,Y)
        else: regroupPair(segs,X,Y)
    
    #assert len(segs[0].scalars) == 1
    result = segs[0].scalars[0]
    return result