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

    functional = 'hse06'

    print('loading reactions...')

    reactions = pickle.load(open('../reaction_data/{}_reactions.p'.format(functional),'rb'))
    
    trange = [300]
    prange = [100]
    
    gs = GraphGenerator(reactions).multidigraph_from_t_and_p_range(trange,prange)
    
    ic = loadfn('initial_concentrations.json')
    
    headers = ['compound','conc /ppm']
    
    print(tabulate(ic.items(),headers=headers))

    samples = [1,10,100,1000,5000,10000]

    start = datetime.now()
    
    sets = Traversal(gs,reactions,ic,trange,prange)

    convergence_data = {}
    with tqdm(total=len(samples),bar_format='{desc:<5.5}{percentage:3.0f}%|{bar:10}{r_bar}') as pbar:
        for sample in samples:
            convergence_data[sample] = sets.graph_sampling_processes(sample_length=sample,path_depth=200,nprocs=32,random_path_depth=True)
            pbar.update(1)
    
    pickle.dump(convergence_data,open('convergence_data_{}.p'.format(functional),'wb'))
