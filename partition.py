# MATTHEW SUTCLIFFE, 2024

import pyzx as zx
import kahypar as kahypar
from .hypergraph import *
from .utils import *

def getK(): return 6 #TEMP

def weighHyperverts(g):
    edges = []
    for e in g.edges(): edges.append(e)
    eWeights = [0]*len(edges)
    #print(edges)
    
    for v in g.vertices():
        if g.phase(v) in (0.25,0.75,1.25,1.75):
            for neigh in g.neighbors(v):
                e = (min(v,neigh), max(v,neigh))
                weight = 1/len(g.neighbors(v))
                eWeights[edges.index(e)] += weight
    return eWeights

def getHyperData(g):
    hypernodes = [] # i.e. the edges of the normal graph
    nets = []
    
    for e in g.edges(): hypernodes.append(e)
    
    for v in g.vertices():
        net = []
        for e in g.neighbors(v):
            hnode = (min(v,e),max(v,e))
            hnodeIndx = hypernodes.index(hnode)
            net.append(hnodeIndx)
        nets.append(net)
    
    #--
    
    eptr = [0]
    eind = []
    
    for net in nets:
        eind.extend(net)
        eptr.append(len(net)+eptr[-1])

    #--

    e_weights = [1]*len(nets)
    v_weights = weighHyperverts(g)

    # upscale to integers...
    weight_factor = 100
    for i in range(len(v_weights)):
        v_weights[i] = int(v_weights[i]*weight_factor)
        if v_weights[i] == 0: v_weights[i] = 1
    
    num_nodes = len(hypernodes)
    num_nets  = len(nets)
    
    return (num_nodes,num_nets,eptr,eind,e_weights,v_weights)

def kahyparPartition(g):
    num_nodes,num_nets,eptr,eind,e_weights,v_weights = getHyperData(g)
    
    k = getK() # 10
    epsilon = 0.5 #0.5 # 0.25 # 0.8 #0.03
    
    hypergraph = kahypar.Hypergraph(num_nodes, num_nets, eptr, eind, k, e_weights, v_weights)
    
    context = kahypar.Context()
    context.loadINIconfiguration(GET_MODULE_PATH()+"/km1_kKaHyPar_sea20.ini") #TEMP
    
    context.setK(k)
    context.setEpsilon(epsilon)
    
    context.setPartitionFileName("PARTITION_OUTPUT.txt")
    context.writePartitionFile(True)
    
    kahypar.partition(hypergraph, context)
    return k
    
def readPartitionData():
    f = open("PARTITION_OUTPUT.txt", "r")
    strMem = f.read()[0:-1]
    f.close()
    mem = [int(ln) if ln.isdigit() else -1 for ln in strMem.split('\n')] # partition membership data of the vertices
    return mem
    
def memToConnectivityGraphs(g, mem):
    #for v in g.vertices(): g.set_phase(v,0)
    
    chars = 'abcdefghijklmnopqrstuvwxyz' #TEMP
    edges = []
    for e in g.edges(): edges.append(e)

    vmems = [set() for i in range(len(g.vertices()))]
    
    parts = max(mem)+1
    gs = []
    for p in range(parts):
        h = g.copy()
        for e in range(len(mem)):
            if mem[e]!=p: h.remove_edge(edges[e])
            else:
                #h.set_phase(edges[e][0],1); h.set_phase(edges[e][1],1) #h.set_phase(edges[e][0],chars[p]); h.set_phase(edges[e][1],chars[p])
                vmems[edges[e][0]].add(chars[p]); vmems[edges[e][1]].add(chars[p])
        gs.append(h)
    
    return gs,vmems
    
def getCutsFromVmems(vmems,nParts,asParams=False):
    cuts = dict()
    for v in range(len(vmems)):
        if len(vmems[v])>1: cuts[v] = vmems[v]    
    if asParams: return cuts
    
    vDict = getVdict(nParts)
    # Convert 'a'->0, 'b'->1, etc....
    vcuts = dict()
    for c in cuts:
        pSet = []
        for p in list(cuts[c]):
            pSet.append(vDict[p])
        vcuts[c] = set(pSet)
    return vcuts
    
def ave_pos(a,b,perc=1/2): return (abs(a-b))*(perc) + min(a,b) # Return position 'perc'%-distance between 2 points

def cut_vertex_param(g,v,ph):
    g  = g.copy()
    gBuffer = g.copy()

    n = len(gBuffer.neighbors(v))
    g.scalar.add_power(-n)
    g.scalar.add_phase_pair(g.phase(v),0,{},{ph}) # account for e^(i*pi*alpha) on right branch
    
    vtype = -1
    match gBuffer.type(v):
        case 1: vtype = 2
        case 2: vtype = 1
        case _: raise ValueError("Attempted illegal cut on boundary vertex "+str(v))

    for i in gBuffer.neighbors(v):
        etype = gBuffer.edge_type((v,i)) # maintain edge type
        qubit = ave_pos(gBuffer.qubit(v),gBuffer.qubit(i),1/2)
        row   = ave_pos(gBuffer.row(v),gBuffer.row(i),1/2)

        if '*' in g.get_params(i): # don't count null-verts as legit neighbors
            g.scalar.add_power(1) # correct for counting it as a neighbour earlier
            continue
        newV = g.add_vertex(vtype,qubit,row,ph) # add and connect the new vertices
        g.add_edge((newV,i),etype)

    g.set_phase(v,'*') # mark for removal #TEMP

    return g
    
def cutAndClean(g,gs,cuts,vmems):
    pChars = 'abcdefghijklmnopqrstuvwxyz'#ABCDEFGHIJKLMNOPQRSTUVWXYZαβγδεϝͷϛζͱηθικλμξϻρσͼφχψωϡͳϸ' #TEMP
    
    # Apply the cuts...
    for i,c in enumerate(cuts):
        for j in range(len(gs)):
            gs[j] = cut_vertex_param(gs[j],c,pChars[i])
    
    # Mark for removal the verts that aren't members of each partition...
    for v in g.vertices():
        for i in range(len(gs)):
            if not pChars[i] in vmems[v]:
                gs[i].set_phase(v,'*')
    
    # Remove the marked verts...
    for i in range(len(gs)):
        for v in g.vertices():
            if '*' in gs[i].get_params(v): gs[i].remove_vertex(v)
    
    # Display partitioned graphs...
    for i in range(len(gs)):
        zx.draw(gs[i],scale=20,labels=True)
        print("T-count:",zx.tcount(gs[i]))

    return gs
    
def accountParamScalars(gs):
    pps = []
    for g in gs:
        for pp in g.scalar.phasepairs:
            pps.append(pp)
    for g in gs: g.scalar.phasepairs = []
    
    ticked_params = set()
    for pp in pps:
        param = list(pp.paramsB)[0]
        if param in ticked_params: continue
        cutphase = pp.alpha/4
        for i in range(len(gs)):
            if param in get_unique_params(gs[i]): # find the first graph-part to contain param p
                gs[i].scalar.add_phase_pair(cutphase,0,{},{param}) # account for e^(i*a*alpha)'s #TEMP (CUTPHASES, CHARS)
                ticked_params.add(param)
                break # go to next p
    return True
    
def get_unique_params(g): #TODO - move this into graph_s.py
    unique_params = set()
    allParams = g.get_all_params()
    for v in allParams:
        unique_params = unique_params.union(allParams[v])
    return unique_params
    
def partition(g):
    k = kahyparPartition(g)
    mem = readPartitionData()
    gs,vmems = memToConnectivityGraphs(g,mem)
    for i in range(1,len(gs)): gs[i].scalar = zx.Graph().scalar # ensure the parent scalar is only counted once
    cuts = getCutsFromVmems(vmems,k)
    hNet = genHNet(cuts,k)
    gs = cutAndClean(g,gs,cuts,vmems)
    accountParamScalars(gs)
    return hNet,gs