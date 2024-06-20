# MATTHEW SUTCLIFFE, 2024

import os
import igraph as ig
import matplotlib.pyplot as plt
import math
import random
import ipywidgets as widgets
from ipywidgets import interact, interactive, fixed, interact_manual
from PIL import Image
from .utils import *

DIR_MODULE = os.path.dirname(__file__) + '/'
#DIR_WORKING = os.path.abspath('') + '/'

#TEMP...
vDict = {
    'A':0, 'B':1, 'C':2, 'D':3, 'E':4, 'F':5, 'G':6, 'H':7, 'I':8, 'J':9, 'K':10, 'L':11, 'M':12, 'N':13, 
    'O':14, 'P':15, 'Q':16, 'R':17, 'S':18, 'T':19, 'U':20, 'V':21, 'W':22, 'X':23, 'Y':24, 'Z':25
}

class hgraph:
    graph = ig.Graph()
    hyperCount = 0 # number of (hyper)edges
    
class HVert:
    label = ''
    active = True
    tcount = 0 # or weight, i.e. T-count of the graph-part (for runtime estimation purposes)
    #scalars = []
    def __init__(self,label,active=True,tcount=0,scalars=[]):
        self.label = label
        self.active = active
        self.tcount = tcount
        #self.scalars = []

class HEdge:
    verts = set()
    weight = 0 # i.e. number of hyperedges connecting these verts
    def __init__(self,verts,weight):
        self.verts = set()
        for v in verts:
            if type(v) == int: self.verts.add(v)
            else: self.verts.add(vDict[v]) #TODO - this should not refer to global vDict!
        self.weight = weight

class HNet:
    hEdges = []
    hVerts = []
    vDict = dict()
    color_dict = dict()
    frame_color_dict = dict()
    
    def __init__(self):
        self.hEdges = []
        self.hVerts = []
        self.vDict = vDict #TODO - this should not refer to global vDict!
        self.color_dict,self.frame_color_dict = initColours()

    def pairWeight(self,a,b,doPrint=False):
        if type(a) != int: a = self.vDict[a]
        if type(b) != int: b = self.vDict[b]
        
        totWeight = 0
        for he in self.hEdges:
            if a in he.verts or b in he.verts:
                totWeight += he.weight
                if doPrint: print(he.verts," = ",he.weight)
        return totWeight

    def clean(self):
        bufferEdgeSets = []
        bufferWeights = []
        for he in self.hEdges:
            if len(he.verts) < 2: # if self-connection
                continue
            elif not he.verts in bufferEdgeSets: # if unique
                bufferEdgeSets.append(he.verts)
                bufferWeights.append(he.weight)
            else: # if repeat
                #print("REPEAT @",he.verts) #TEMP
                bufferWeights[bufferEdgeSets.index(he.verts)] += he.weight
        self.hEdges = []
        for i in range(len(bufferEdgeSets)):
            self.hEdges.append(HEdge(bufferEdgeSets[i],bufferWeights[i]))
    
    def fusePair(self,a,b):
        if type(a) != int: a = self.vDict[a]
        if type(b) != int: b = self.vDict[b]
        
        for he in self.hEdges:
            if b in he.verts:
                he.verts.remove(b)
                he.verts.add(a)
        self.clean()
        self.hVerts[a].label += self.hVerts[b].label

        # Sort the string: e.g. CDAEB -> ABCDE...
        lbl = ''
        buffer = list(self.hVerts[a].label)
        buffer.sort()
        for i in buffer: lbl += i
        self.hVerts[a].label = lbl
        self.hVerts[a].tcount += self.hVerts[b].tcount
        self.hVerts[b].active = False

    def findAdjPairs(self):
        uniqueAdjPairs = []
        for he in self.hEdges:
            for a in he.verts:
                for b in he.verts:
                    newPair = {a,b}
                    if (not newPair in uniqueAdjPairs) and (a != b): uniqueAdjPairs.append({a,b})
        return uniqueAdjPairs
        
    def findCheapestFusePair(self):
        # IF THERE ARE FREE (DISCONNECTED) HVERTS...
        freeNodes = self.findFreeNodes()
        if len(freeNodes) > 1: # if multiple free nodes, return a pair of them
            return sorted({freeNodes[0],freeNodes[1]})
        elif len(freeNodes) > 0: # if only one free node, pair it with the cheapest non-free node
            a = freeNodes[0]
            best_b = -1
            cheapest = 999999
            for b,hv in enumerate(self.hVerts):
                if not hv.active: continue
                if b==a: continue
                cost = self.pairWeight(a,b)
                if cost < cheapest:
                    cheapest = cost
                    best_b = b
            return sorted({a,best_b})
            
        # IF THERE ARE NO FREE (DISCONNECTED) HVERTS...  
        adjPairs = self.findAdjPairs()
        best_w = 999999
        best_ab = {}
        for pair in adjPairs:
            a = list(pair)[0]
            b = list(pair)[1]
            w = self.pairWeight(a,b)
            #print(a,b," => ",w)
            if w < best_w:
                best_w = w
                best_ab = {a,b}
        return sorted(best_ab)#,best_w

    def getActiveVerts(self):
        activeVerts = []
        for i,hv in enumerate(self.hVerts):
            if hv.active: activeVerts.append(i)
        return activeVerts

    def stepReduce(self,output=True,draw=False):
        cheapestPair = self.findCheapestFusePair()
        w = self.pairWeight(cheapestPair[0],cheapestPair[1],doPrint=False)
        if output: print("Fuse",self.hVerts[cheapestPair[0]].label,"and",self.hVerts[cheapestPair[1]].label,"\t=> cost:",w)
        self.fusePair(cheapestPair[0],cheapestPair[1])
        if draw: self.draw()
        return (cheapestPair,w)
    
    def fullReduce(self,output=True,drawSteps=False,doHighlight=True):
        topCost = 0
        while len(self.getActiveVerts()) > 1:
            highlightedPair = None
            if doHighlight: highlightedPair = self.findCheapestFusePair()
            if drawSteps: self.draw(highlightedVerts=highlightedPair)
            cp,w = self.stepReduce(output,draw=False)
            if w > topCost: topCost = w
        if drawSteps: self.draw()
        return topCost

    def fullReduceWidget(self,doHighlight=True):
        outFigs = []
        outNextMoves = []
        n = 0
        while True: #len(self.getActiveVerts()) > 1:
            if n>0: self.stepReduce(False)
            if len(self.getActiveVerts()) == 1:
                g,vs = drawHypergraph(self)
                strNextMove = "[No more moves]"
                outFigs.append(ig.plot(g, **vs,target='outfig'+str(n)+'.png'))
                outNextMoves.append(strNextMove)
                break
            cheapestPair = self.findCheapestFusePair()
            highlightPair = None
            if doHighlight: highlightPair = cheapestPair
            g,vs = drawHypergraph(self,highlightPair)
            strNextMove = "Fuse "+str(self.hVerts[cheapestPair[0]].label)+" and "+str(self.hVerts[cheapestPair[1]].label)
            strNextMove += "\t=> cost: "+str(self.pairWeight(cheapestPair[0],cheapestPair[1],doPrint=False))
            outFigs.append(ig.plot(g, **vs,target='outfig'+str(n)+'.png'))
            outNextMoves.append(strNextMove)
            n += 1
        #return outFigs,outNextMoves
        # OUTPUT AS WIDGET...
        a = widgets.IntSlider(min=0,max=len(outFigs)-1)
        ui = a
        def f_wid(a):
            display(outFigs[a])
            print(outNextMoves[a]) # TODO - print the highlighted weight (for next move) here
            #[probably also maintain consistent colourings on edges]
            #Resize the labels too, based on no. of chars
        out = widgets.interactive_output(f_wid, {'a': a})
        display(ui, out)

    def getEdgesUnpacked(self):
        es = []
        for he in self.hEdges:
            for i in range(he.weight):
                es.append(tuple(he.verts))
        return es
    
    def getVertId(self,strVert):
        for i,hv in enumerate(self.hVerts):
            if hv.label == strVert: return i
        return -1

    def copy(self):
        newNet = HNet()
        newNet.vDict  = self.vDict.copy()
        for he in self.hEdges: newNet.hEdges.append(HEdge(he.verts.copy(),he.weight))
        for hv in self.hVerts: newNet.hVerts.append(HVert(hv.label,hv.active,hv.tcount))
        return newNet

    def draw(self,highlightedVerts=None):
        hg,visual_style = drawHypergraph(self,highlightedVerts)
        display(ig.plot(hg, **visual_style))

    def save(self,fileName):
        f = open(fileName+".txt", "w")
        strLine = ''
        for hv in self.hVerts:
            strLine += hv.label+','
            strLine += str(int(hv.active))+','
            strLine += str(hv.tcount)+';'
        f.write(strLine[:-1])
        for he in self.hEdges:
            strLine = '\n'
            for v in he.verts: strLine += str(v)+','
            f.write(strLine[:-1])
            f.write(';'+str(he.weight))
        f.close()
        print("Saved to:",os.path.realpath(f.name))
        return True
        
    def findFreeNodes(self):
        adjPairs = self.findAdjPairs()
        freeNodes = []
        for i,hv in enumerate(self.hVerts):
            i_is_free = True
            if not hv.active: continue
            for pair in adjPairs:
                if i in pair:
                    i_is_free = False
                    break
            if i_is_free: freeNodes.append(i)
        return freeNodes
    
def getVdict(nParts, inCaps = False): #TEMP (just use charsUpper.find()?)
    charsUpper = CHARS_SEGMENTS
    charsLower = CHARS_PARAMS
    vDictUpper = dict()
    vDictLower = dict()
    
    for i in range(nParts):
        vDictUpper[charsUpper[i]] = i
        vDictLower[charsLower[i]] = i
        
    if inCaps: return vDictUpper
    else: return vDictLower
    
def addEdge(hg,e,w=0):
    g = hg.graph
    hg.hyperCount += 1
    g.add_edge(e[0],e[1])
    hyperedgeIndexes = g.es["hyperedgeIndx"]; hyperedgeIndexes[-1] = hg.hyperCount; g.es["hyperedgeIndx"] = hyperedgeIndexes
    if w>0: weights = g.es["weight"]; weights[-1] = w; g.es["weight"] = weights

# Insert the "dummy" nodes accordingly...
def addHyperedge(hg, nodes, layout, w=0):
    g = hg.graph
    hg.hyperCount += 1
    g.add_vertex()
    labels = g.vs["label"]; labels[-1] = ""; g.vs["label"] = labels
    isDummy = g.vs["isDummy"]; isDummy[-1] = 1; g.vs["isDummy"] = isDummy
    #if w>0: labels = g.vs["label"]; labels[-1] = w; g.vs["label"] = labels #TEMP
    newDummy = len(g.vs)-1
    
    xs = [layout[i][0] for i in (nodes)]
    ys = [layout[i][1] for i in (nodes)]
    x = (max(xs)-min(xs))/2 + min(xs)
    y = (max(ys)-min(ys))/2 + min(ys)
    layout.append((x,y))
    for i in nodes:
        g.add_edge(newDummy,i)
        hyperedgeIndexes = g.es["hyperedgeIndx"]; hyperedgeIndexes[-1] = hg.hyperCount; g.es["hyperedgeIndx"] = hyperedgeIndexes
        hypernodeIndexes = g.vs["hypernodeIndx"]; hypernodeIndexes[-1] = hg.hyperCount; g.vs["hypernodeIndx"] = hypernodeIndexes
    
    if w>0: weights = g.es["weight"]; weights[-1] = w; g.es["weight"] = weights
    
def initIg(n_verts): #TEMP
    g = ig.Graph()
    g.add_vertices(n_verts)
    
    g.vs["label"] = ["A", "B", "C", "D", "E", "F", "G", "H"]
    g.vs["isDummy"] = [0,0,0,0,0,0,0,0]
    g.vs["hypernodeIndx"] = [0,0,0,0,0,0,0,0]
    g.vs["dead"] = [0,0,0,0,0,0,0,0]
    g.es["hyperedgeIndx"] = []
    g.es["weight"] = []
    #g.es["label"] = []

    hg = hgraph()
    hg.graph = g
    
    return hg
    
def bundleEdges(edges):
    checkedEdges = []
    edgeWeights = []
    
    # Sort the node-sets in the edge list...
    for i in range(len(edges)):
        e = edges[i]
        e = list(e)
        e.sort()
        edges[i] = e
    
    for e in edges:
        if e in checkedEdges: continue
        checkedEdges.append(e)
        edgeWeights.append(edges.count(e))
    
    edges = checkedEdges
    return edges,edgeWeights
    
def constructEdgeBundles(hg,edges,layout,edgeWeights):
    for i in range(len(edges)):
        e = edges[i]
        w = edgeWeights[i]
        if len(e) < 3: addEdge(hg,e,w)
        else: addHyperedge(hg,e,layout,w)
    return hg
    
def randCol(): return '#{:02x}{:02x}{:02x}'.format(random.randrange(0,256), random.randrange(0,256), random.randrange(0,256))

# The colour palette is randomly generated (but overwritten by palette.png with as many colours as provided - any thereafter will still be random):
# [Note the first colour is reserved for the main vertices interior, and the second for the main vertices border]
def initColours():
    colCount = 256
    color_dict = ["#FFB8D7"]*colCount
    frame_color_dict = ["#FF0000"]*colCount
    for i in range(2,colCount):
        col = randCol()
        color_dict[i] = col
        frame_color_dict[i] = col
    
    #--
    
    im = Image.open(DIR_MODULE+'palette.png')
    width, height = im.size
    pix = im.load()
    
    colParts = pix[0,0]; color_dict[0]       = '#{:02x}{:02x}{:02x}'.format(colParts[0],colParts[1],colParts[2])
    colParts = pix[1,0]; frame_color_dict[0] = '#{:02x}{:02x}{:02x}'.format(colParts[0],colParts[1],colParts[2])
    
    for x in range(2,width):
        colParts = pix[x,0]
        color_dict[x-1] = '#{:02x}{:02x}{:02x}'.format(colParts[0],colParts[1],colParts[2])
        frame_color_dict[x-1] = '#{:02x}{:02x}{:02x}'.format(colParts[0],colParts[1],colParts[2])
    
    im.close()
    return (color_dict,frame_color_dict)
    
def prepHighlight(g,highlightVerts,doHighlightNeighbours=False):
    vs = [0]*len(g.vs)
    es = [0]*len(g.es)
    
    for i in range(len(g.vs)):
        if g.vs["isDummy"][i] == 0: vs[i] = 1
    
    for vert in highlightVerts:
        for v in g.vs[vert].neighbors():
            i = v.index
            es[g.get_eid(vert,i)] = 1
            if g.vs["isDummy"][i]:
                vs[i] = 2
                for u in g.vs[i].neighbors():
                    j = u.index
                    vs[j] = 3 if doHighlightNeighbours else 1
                    es[g.get_eid(i,j)] = 1
            else: vs[i] = 3 if doHighlightNeighbours else 1
    
    for vert in highlightVerts: vs[vert] = 3 # 0 = background dummy, 1 = background genuine, 2 = highlight dummy, 3 = highlight genuine
    for vert in range(len(g.vs)):
        if g.vs["dead"][vert] == 1: vs[vert] = -1
    
    return vs,es
    
def resetStyle(g,layout,color_dict,frame_color_dict):
    #layout = g.layout("circle")
    size_dict  = {0: 30, 1: 6, 2: 0}
    
    visual_style = {}
    visual_style["vertex_size"]  = [size_dict[isDummy]  for isDummy in g.vs["isDummy"]]
    visual_style["vertex_color"] = [color_dict[hypernodeIndx] for hypernodeIndx in g.vs["hypernodeIndx"]]
    visual_style["vertex_label_color"] = ["#000000" for i in range(len(g.vs))]
    visual_style["vertex_frame_color"] = [frame_color_dict[hypernodeIndx] for hypernodeIndx in g.vs["hypernodeIndx"]]
    visual_style["vertex_label"] = [lbl for lbl in g.vs["label"]]
    visual_style["edge_width"] = [2 for i in range(len(g.es))]
    visual_style["edge_color"] = [color_dict[hyperedgeIndx] for hyperedgeIndx in g.es["hyperedgeIndx"]]
    visual_style["edge_label_color"] = [color_dict[hyperedgeIndx] for hyperedgeIndx in g.es["hyperedgeIndx"]] #TEMP
    visual_style["edge_label"] = [w for w in g.es["weight"]]
    
    visual_style["layout"] = layout
    visual_style["bbox"] = (300, 300)
    visual_style["margin"] = 20
    
    for i in range(len(g.vs)):
        if g.vs["dead"][i] == 1:
            visual_style["vertex_size"][i] = 0
            visual_style["vertex_label"][i] = " "
    
    return visual_style
    
# vs: 0 = background dummy, 1 = background genuine, 2 = highlight dummy, 3 = highlight genuine
# es: 0 = background edge,  1 = highlight edge
def percToHex(perc): return '{:02x}'.format(int(perc*255))
    
def drawHighlight(g,selectedVerts):
    vs,es = prepHighlight(g,selectedVerts)
    visual_style = resetStyle(g,layout)
    
    opacityMid = 0.2
    opacityLow = 0.1
    hexOpacMid = percToHex(opacityMid)
    hexOpacLow = percToHex(opacityLow)
    
    for i in range(len(g.vs)):
        match vs[i]:
            case 0: # background dummy
                visual_style["vertex_size"][i]         = 1
                visual_style["vertex_color"][i]       += hexOpacLow
                visual_style["vertex_frame_color"][i] += hexOpacLow
                visual_style["vertex_label_color"][i] += hexOpacLow
            case 1: # background genuine
                visual_style["vertex_color"][i]       += hexOpacMid
                visual_style["vertex_frame_color"][i] += hexOpacMid
                visual_style["vertex_label_color"][i] += hexOpacMid
    
    for i in range(len(g.es)):
        if es[i] == 0: # background edge
            visual_style["edge_width"][i] = 0.1
            visual_style["edge_label_color"][i] += hexOpacLow
        if es[i] == 1: # highlight edge
            visual_style["edge_width"][i] = 2.5

    display(ig.plot(g, **visual_style))
    
# vs: 0 = background dummy, 1 = background genuine, 2 = highlight dummy, 3 = highlight genuine, -1 = dead
# es: 0 = background edge,  1 = highlight edge
def percToHex(perc): return '{:02x}'.format(int(perc*255))

def prepHightlightStyle(g,vs,es,visual_style):
    opacityMid = 0.2
    opacityLow = 0.1
    hexOpacMid = percToHex(opacityMid)
    hexOpacLow = percToHex(opacityLow)
    
    for i in range(len(g.vs)):
        match vs[i]:
            case 0: # background dummy
                visual_style["vertex_size"][i]         = 1
                visual_style["vertex_color"][i]       += hexOpacLow
                visual_style["vertex_frame_color"][i] += hexOpacLow
                visual_style["vertex_label_color"][i] += hexOpacLow
            case 1: # background genuine
                visual_style["vertex_color"][i]       += hexOpacMid
                visual_style["vertex_frame_color"][i] += hexOpacMid
                visual_style["vertex_label_color"][i] += hexOpacMid
            case -1: # dead vertex
                visual_style["vertex_size"][i]         = 0
                visual_style["vertex_color"][i]       += "00"
                visual_style["vertex_frame_color"][i] += "00"
                visual_style["vertex_label_color"][i] += "00"
    
    for i in range(len(g.es)):
        if es[i] == 0: # background edge
            visual_style["edge_width"][i] = 0.1
            visual_style["edge_label_color"][i] += hexOpacLow
        if es[i] == 1: # highlight edge
            visual_style["edge_width"][i] = 2.5

    return visual_style
    
def loadHNet(fileName):
    hNet = HNet()
    f = open(fileName+".txt", "r")
    strData = f.read()
    f.close()
    lines = strData.split('\n')
    vData = lines[0].split(';')
    for v in vData:
        vParts = v.split(',')
        hVert = HVert(vParts[0],bool(int(vParts[1])),int(vParts[2]))
        hNet.hVerts.append(hVert)
    for i in range(1,len(lines)):
        line = lines[i]
        eData = line.split(';')
        eSet = eData[0].split(',')
        verts = []
        for e in eSet: verts.append(int(e))
        hEdge = HEdge(verts,int(eData[1]))
        hNet.hEdges.append(hEdge)
    return hNet
    
def drawHypergraph(hNet,highlightedVerts=None):
    n_verts = len(hNet.hVerts)
    
    # Arrange the "genuine" nodes into a circle...
    layout = [(-1,-1) for i in range(n_verts)]
    for i in range(n_verts):
        theta = i/n_verts
        layout[i] = (math.cos(theta*2*math.pi), math.sin(theta*2*math.pi))

    #--

    edges = hNet.getEdgesUnpacked() # [(0,1),(0,1,2),(0,2),(0,2),(0,2),(0,2),(0,2),(0,2),(0,2),(0,2,3),(0,2,3),(0,2,3),(2,3),(2,3),(2,3),(2,3)]

    hg = initIg(n_verts)
    edges,edgeWeights = bundleEdges(edges)
    hg = constructEdgeBundles(hg,edges,layout,edgeWeights)
    g = hg.graph
    
    labels = ['' for i in range(len(g.vs))]
    for i in range(len(hNet.hVerts)): labels[i] = labels[i] = hNet.hVerts[i].label
    g.vs["label"] = labels
    
    deads = g.vs["dead"]
    for i,hv in enumerate(hNet.hVerts):
            deads[i] = int(not hv.active)
    g.vs["dead"] = deads
    #deads = g.vs["dead"]; deads[a] = 1; g.vs["dead"] = deads #TEMP
    #labels = g.vs["label"]; labset = [labels[a],labels[b]]; labset.sort(); labels[b] = ''.join(labset); labels[a] = ''; g.vs["label"] = labels #TEMP
    
    visual_style = resetStyle(g,layout,hNet.color_dict,hNet.frame_color_dict)
    if not highlightedVerts == None:
        vs,es = prepHighlight(g,highlightedVerts,False)
        visual_style = prepHightlightStyle(g,vs,es,visual_style)
    
    #ig.plot(g, **visual_style)
    return g,visual_style
    
def genHNet(cuts,nParts,gs=[]):
    vDict = getVdict(nParts,True)
    
    hEdges = []
    for c in cuts:
        pSet = cuts[c]
        isRepeat = False
        for he in hEdges:
            if he.verts == pSet:
                he.weight += 1
                isRepeat = True
                break
        if not isRepeat:
            hEdge = HEdge(pSet,1)
            hEdges.append(hEdge)
    
    hVerts = []
    for i,v in enumerate(vDict):
        tc = 0
        if len(gs) > 0: tc = zx.tcount(gs[i])
        hVerts.append(HVert(v,tcount=tc))
    
    hNet = HNet()
    hNet.vDict = vDict
    hNet.hEdges = hEdges
    hNet.hVerts = hVerts
    return hNet