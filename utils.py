# MATTHEW SUTCLIFFE, 2024

import os

def GET_MODULE_PATH():
    return os.path.dirname(__file__)
    
def get_unique_params(g): #TODO - move this into pyzx/graph_s.py
    unique_params = set()
    allParams = g.get_all_params()
    for v in allParams:
        unique_params = unique_params.union(allParams[v])
    return unique_params