#!/usr/env/bin python3
import warnings
import pandas as pd
from datetime import datetime
import pickle 
from tabulate import tabulate
from arcs.analysis_functions import *
from arcs.setup_functions import *
from monty.serialization import loadfn

if __name__ == "__main__":
    warnings.simplefilter('ignore')
    
    print('''   \   _ \  __|  __| 
  _ \    / (   \__ \ 
_/  _\_|_\\___|____/''' )

    functional = 'pbe'

    print('\n loading reactions...')

    reactions = pickle.load(open('../../../reaction_data/{}_reactions.p'.format(functional),'rb'))
    
    trange = [300]
    prange = [100]
    
    gs = GraphGenerator(reactions).multidigraph_from_t_and_p_range(trange,prange)
    
    ic = loadfn('initial_concentrations.json')

    headers = ['compound','conc /ppm']
    
    print(tabulate(ic.items(),headers=headers))
    
    sample_length = 10000 
    no2_concs = np.linspace(1e-3,0,8) 
    
    data = {}
    for i in no2_concs:
        ic['SO2'] = i
        sets = Traversal(gs,reactions,ic,trange,prange)
        t_and_p_data = sets.graph_sampling_processes(sample_length=sample_length,path_depth=200,nprocs=32,random_path_depth=True)
        data[i] = t_and_p_data

    print("Writing to file...")
    
    pickle.dump(data,open('concs_and_stats_{}_{}.p'.format(functional,sample_length),'wb'))

    print("Data written to concs_and_stats_{}_{}.p".format(functional,sample_length)) 
