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
   
    
def get_reaction_statistics(t_and_p_data):
    trange = list(t_and_p_data.keys()) #test
    prange = list(t_and_p_data[trange[0]].keys())
    equations = {}
    for T in trange:
        eqs_p  = {}
        for P in prange:
            eqs = []
            for x in t_and_p_data[T][P].keys():
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

    try:
        dict_of_dataframes = {T:{P:get_dataframes(equations[T][P]) 
                             for P in prange} 
                          for T in trange}
    except:
        dict_of_dataframes = {T:{P:[] for P in prange}
                              for T in trange}
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

def get_mean_change_in_data(t_and_p_data,percentage=True):
    
    if isinstance(list(t_and_p_data)[0],str):
        t_and_p_data = str_to_int_dict(t_and_p_data)
        

    final = {T:
                {P:
                np.mean(pd.DataFrame({x:t_and_p_data[T][P][x]['data'] for x in t_and_p_data[T][P]})[1:].T) 
                for P in t_and_p_data[T]}
                for T in t_and_p_data}
        
    init = {T:
            {P:
                pd.DataFrame({x:t_and_p_data[T][P][x]['data'] for x in t_and_p_data[T][P]})[0]
                for P in t_and_p_data[T]}
            for T in t_and_p_data}
        
    if not percentage == True: 
        md = {T:
              {P:
               final[T][P] - init[T][P]
               for P in t_and_p_data[T]}
              for T in t_and_p_data}
    else:
        md = {T:
              {P:
               ((final[T][P] - init[T][P])/init[T][P]) * 100
               for P in t_and_p_data[T]}
              for T in t_and_p_data}
        
    return(pd.DataFrame(md)) 


            
