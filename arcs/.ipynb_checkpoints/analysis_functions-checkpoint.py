# tidy up these imports
import random
from chempy.equilibria import Equilibrium,EqSystem
from chempy import Substance
import copy
import networkx as nx
import itertools as it
#from tqdm.notebook import tqdm
from tqdm import tqdm
from collections import defaultdict
from numpy.random import choice

import os
import warnings

import pathos.multiprocessing as multiprocessing
from pathos.pools import ProcessPool
from pathos.pp import ParallelPool
from pathos.serial import SerialPool
from datetime import datetime
import random
import math
import numpy as np
import pandas as pd 
from pathos.helpers import mp as pmp
import queue
from collections import Counter
   
    
def get_reaction_statistics(t_and_p_data):
    #trange = list(t_and_p_data.keys()) #test
    #prange = list(t_and_p_data[trange[0]].keys())
    equations = {}
    for T in t_and_p_data:
        eqs_p  = {}
        for P in t_and_p_data[T]:
            eqs = []
            for x in t_and_p_data[T][P]:
                if t_and_p_data[T][P][x]['equation_statistics']:
                    eqs.append(t_and_p_data[T][P][x]['equation_statistics'])
            eqs_p[P] = eqs
        equations[T] = eqs_p
    
    def get_dataframes(list_of_equations):
        appearances = defaultdict(int)
        for sample in list_of_equations:
            for i in sample:
                appearances[i] += 1
    
        equation_statistics = {}
        for equation,frequency in appearances.items():
            eq,k = equation.split(';')
            equation_statistics[eq] = {'k':k.split('\n')[0],'frequency':frequency}
        d = pd.DataFrame(equation_statistics).T.sort_values(by='frequency',ascending=False)
        return(d)

    dict_of_dataframes = {}
    for T in equations:
        dfs = {}
        for P in equations[T]:
            try:
                dfs[P] = get_dataframes(equations[T][P])
            except:
                dfs[P] = []
        dict_of_dataframes[T] = dfs

    return(dict_of_dataframes)

def str_to_int_dict(data):
    new_dict = {}
    for t in data:
        tdata = {}
        for p in data[t]:
            pdata = {}
            for x in data[t][p]:
                pdata[int(x)] = {'data':data[t][p][x]['data']}
            tdata[int(p)] = pdata
        new_dict[int(t)] = tdata
    return(new_dict)

def get_mean_change_in_data(t_and_p_data,percentage=False):

    tr = list(t_and_p_data.keys())
    pr = list(t_and_p_data[tr[0]].keys())
    xr = list(t_and_p_data[tr[0]][pr[0]].keys())
    
    if isinstance(xr[0],str):
        t_and_p_data = str_to_int_dict(t_and_p_data)
        

    final = {T:
                {P:
                np.mean(pd.DataFrame({x:t_and_p_data[T][P][x]['data'] for x in t_and_p_data[T][P].keys()})[1:].T) 
                for P in t_and_p_data[T].keys()}
                for T in t_and_p_data.keys()}
        
    init = {T:
            {P:
                pd.DataFrame({x:t_and_p_data[T][P][x]['data'] for x in t_and_p_data[T][P].keys()})[0]
                for P in t_and_p_data[T].keys()}
            for T in t_and_p_data.keys()}
        
    if not percentage == True: 
        md = {T:
              {P:
               final[T][P] - init[T][P]
               for P in t_and_p_data[T].keys()}
              for T in t_and_p_data.keys()}
    else:
        md = {T:
              {P:
               ((final[T][P] - init[T][P])/init[T][P]) * 100
               for P in t_and_p_data[T].keys()}
              for T in t_and_p_data.keys()}
        
    return(pd.DataFrame(md)) 

def reaction_paths(data,T,P,max_rows=10):
    '''currently chooses the top reaction, and only does what comes after'''
    
    stats = {int(x):{y:d.split(';')[0] for y,d in enumerate(data[T][P][x]['equation_statistics']) if d} for x in data[T][P]}
    top_reaction = str(get_reaction_statistics(data)[T][P].head(1).index[0])
    
    valid_samples = []
    for x in stats:
        for y in stats[x]:
            if top_reaction in stats[x][y]:
                valid_samples.append(x)

    paths_2_length = []
    for x in valid_samples:
        if len(stats[x]) > 1:
            for y in stats[x]:
                if stats[x][y] == top_reaction:
                    try:
                        paths_2_length.append(stats[x][y]+':'+stats[x][y+1])
                    except:
                        paths_2_length.append(stats[x][y-1]+':'+stats[x][y])
    frequencies = Counter(paths_2_length)
    freq_sort = {frequencies[f]:{x:d for x,d in enumerate(f.split(':'))} for x,f in enumerate(frequencies)}
    df = pd.DataFrame(dict(reversed(sorted(freq_sort.items())))).T.head(max_rows).reset_index()
    
    df.columns = 'frequency','reaction 1','reaction 2'
    df.set_index('frequency')
    return({int(T):{int(P):df}})
from chempy.equilibria import EqSystem,Equilibrium

def get_eqsys_from_top_paths(top_paths,applied_reactions,T,P):
    eqsys = {}
    for i in top_paths[T][P].T:
        r_1 = Equilibrium.from_string(top_paths[T][P].T[i]['reaction 1'])
        r_2 = Equilibrium.from_string(top_paths[T][P].T[i]['reaction 2'])
        for x in applied_reactions[T][P]:
            if applied_reactions[T][P][x]['e'].string() == r_1.string():
                k_1 = applied_reactions[T][P][x]['k']
            elif applied_reactions[T][P][x]['e'].string() == r_2.string():
                k_2 = applied_reactions[T][P][x]['k']
        eqsys[i] = EqSystem([Equilibrium(r_1.reac,r_1.prod,k_1),Equilibrium(r_2.reac,r_2.prod,k_2)])
    return(eqsys)

            
