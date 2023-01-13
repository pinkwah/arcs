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
   
def get_reaction_statistics(data):
        
    def _get_stats(list_of_equations):
        appearances = defaultdict(int)
        for sample in list_of_equations:
            for i in sample:
                appearances[i] += 1

        equation_statistics = {}
        for equation,frequency in appearances.items():
            eq,k = equation.split(';')
            equation_statistics[eq] = {'k':k.split('\n')[0],'frequency':frequency}
        try:
            d = pd.DataFrame(equation_statistics).T.sort_values(by='frequency',ascending=False)
            d = d.reset_index()
            d.T['index'] = 'reaction'
            d = d.to_dict()
        except:
            d = {}
        return(d)
    
    eqt = {}
    for T in data:
        eqp = {}
        for P in data[T]:
            equations = []
            for x in data[T][P]:
                eqs = data[T][P][x]['equation_statistics']
                if eqs:
                    equations.append(eqs)
            try:
                eqp[float(P)] = _get_stats(equations)
            except:
                eqp[float(P)] = []
        eqt[float(T)] = eqp
   
    return(eqt)

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

def get_mean_change_in_data(data):
    tr = list(data.keys())
    pr = list(data[tr[0]].keys())
    xr = list(data[tr[0]][pr[0]].keys())
    if isinstance(xr[0],str):
        data = str_to_int_dict(data)
        final = {float(T):
                 {float(P):
                 np.mean(pd.DataFrame({x:data[T][P][x]['data'] 
                                       for x in data[T][P].keys()})[1:].T)
                 for P in data[T].keys()}
                 for T in data.keys()}
        init = {float(T):
                 {float(P):
                     pd.DataFrame({x:data[T][P][x]['data'] 
                                   for x in data[T][P].keys()})[0]
                     for P in data[T].keys()}
                 for T in data.keys()}

    md = {float(T):
          {float(P):
           final[T][P] - init[T][P]
           for P in data[T].keys()}
          for T in data.keys()}
            
    return(md)# needs to be a dict

def reaction_paths(data,max_rows=10):
    '''currently chooses the top reaction, and only does what comes after'''

    def eqsys_from_path(pathsstats):
        new  = {}
        new['paths'] = {}
        new['k'] = {}
        new['frequency'] = pathsstats['frequency']
        for i in pathsstats['frequency']:
            r_1,k_1 = pathsstats['reaction 1'][i].split(';')
            k_1 = float(k_1.split('k=')[1])
            r_2,k_2 = pathsstats['reaction 2'][i].split(';')
            k_2 = float(k_2.split('k=')[1])
    
            string_1 = r_1 + '\n' + r_2 
            string_2 = str(k_1) + '\n' + str(k_2)
            new['paths'][i] = string_1
            new['k'][i] = string_2
        return(new)

    df1 = {}
    for T in data:
        df2 = {}
        for P in data[T]:
            stats = {int(x):{y:{'reaction':d.split(';')[0],'k':d.split(';')[1]} for y,d in enumerate(data[T][P][x]['equation_statistics']) if d} for x in data[T][P]}
            try:
                top_reaction = str(get_reaction_statistics(data)[float(T)][float(P)]['index'][0])
            except:
                top_reaction = None
    
            valid_samples = []
            for x in stats:
                if stats[x]:
                    for y in stats[x]:
                        if top_reaction in stats[x][y]['reaction']:
                             valid_samples.append(x)

            p2l = []
            for x in valid_samples:
                if len(stats[x]) > 1:
                    for y in stats[x]:
                        if stats[x][y]['reaction'] == top_reaction:
                            try:
                                p2l.append(stats[x][y]['reaction']+' ; k='+stats[x][y]['k'].split('\n')[0]+':'+stats[x][y+1]['reaction']+' ; k='+stats[x][y+1]['k'].split('\n')[0])
                            except:
                                p2l.append(stats[x][y-1]['reaction']+' ; k='+stats[x][y-1]['k'].split('\n')[0]+':'+stats[x][y]['reaction']+' ; k='+stats[x][y]['k'].split('\n')[0])
            try:
                frequencies = Counter(p2l)
                freq_sort = {frequencies[f]:{x:d for x,d in enumerate(f.split(':'))} for x,f in enumerate(frequencies)}
                df = pd.DataFrame(dict(reversed(sorted(freq_sort.items())))).T.head(max_rows).reset_index()
    
                df.columns = 'frequency','reaction 1','reaction 2'
                df.set_index('frequency')
                dict_ = df.to_dict()
                df2[float(P)] = eqsys_from_path(dict_)
            except:
                df2[float(P)] = {'frequency':[None],'paths':[None],'k':[None]}

        df1[float(T)] = df2
    return(df1)
