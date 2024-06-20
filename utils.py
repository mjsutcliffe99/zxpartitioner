# MATTHEW SUTCLIFFE, 2024

import os

EST_ALPHA = 0.5 # Estimate of the typical value of alpha for the decomposition(s) used (e.g. BSS has alpha~=0.47, but in practise is closer to ~0.42, so should cite the latter here)

CHARS_SEGMENTS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyzαβγδεϝͷϛζͱηθικλμξϻρσͼφχψωϡͳϸ' #TEMP (should probably avoid using the same chars between segments and params)
CHARS_PARAMS   = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZαβγδεϝͷϛζͱηθικλμξϻρσͼφχψωϡͳϸ'

def GET_MODULE_PATH():
    return os.path.dirname(__file__)
    
def get_unique_params(g): #TODO - move this into pyzx/graph_s.py
    unique_params = set()
    allParams = g.get_all_params()
    for v in allParams:
        unique_params = unique_params.union(allParams[v])
    return unique_params