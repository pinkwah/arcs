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

class Traversal:
    def __init__(self,graph,reactions,concs,trange,prange,**kwargs):
        self.graph = graph
        self.reactions = reactions
        self.concs = copy.deepcopy(concs)
        #self.progress = kwargs['progress']
        self.trange = trange
        self.prange = prange
        pass
        
    def _get_weighted_random_compound(self,T,P,co2=False,force_selection=None):     # doesnt' include looping

        nodes = [n for n in self.graph[T][P].nodes() if isinstance(n,str)]
        temp_concs = copy.deepcopy(self.concs)            

        if not force_selection == None:
            new_temp_concs = {}
            force_selection = force_selection + ['CO2']
            for c in temp_concs:
                if c in force_selection:
                    new_temp_concs[c] = temp_concs[c]

            temp_concs = new_temp_concs

        if co2 == False:
            del temp_concs['CO2']
            
        available = [choice(list(temp_concs.keys()),
                len(temp_concs),
                p=[x/sum(temp_concs.values()) for x in temp_concs.values()])[0]]

        source = random.choice([n for n in nodes if n in available])
        return(source)

    def _random_choice_unconnected(self,T,P,force_direct=False,co2=False): # currently randomly disjointed reactions that are weighted
        nodes = [n for n in self.graph[T][P].nodes() if isinstance(n,str)]
        if force_direct == True:
            pstring = [0,1,2]
            while len(pstring) > 2:
                source = self._get_weighted_random_compound(T,P,co2=co2,force_selection=None) 
                target = random.choice(nodes)
                p = nx.shortest_path(self.graph[T][P],source,target)
                pstring = [n for n in p if isinstance(p,str)]
        else:
                source = self._get_weighted_random_compound(T,P)
                target = random.choice(nodes)
                p = nx.shortest_path(self.graph[T][P],source,target)
        return(p)

    def _random_choice_connected(self,T,P,force_direct=False,previous_index=None,co2=False): # this will be joined - I think we can also make a ranking of potential reactions based upon components in the stream as well 
        if previous_index == None:
            raise ValueError('no previous compound selected')
        nodes = [n for n in self.graph[T][P].nodes() if isinstance(n,str)]
        if force_direct == True:
            pstring = [0,1,2]
            while len(pstring) > 2:
                present = [c for c in list(self.reactions[T][P][previous_index]['e'].reac) + list(self.reactions[T][P][previous_index]['e'].prod) ] # this should probably be weighted according to stoichiometry i.e. 2CO2 + H2O = [CO2, CO2, H2O]
                source = self._get_weighted_random_compound(T,P,co2=co2,force_selection=present)
                target = random.choice(nodes) # the next path will be random 
                p = nx.shortest_path(self.graph[T][P],source,target)
                pstring = [n for n in p if isinstance(p,str)]
        else:
                source = self._get_weighted_random_compound(T,P)
                target = random.choice(nodes)
                p = nx.shortest_path(self.graph[T][P],source,target)
        return(p)


    def random_walk(self,T,P,path_depth=10):
        nodes = [n for n in self.graph[T][P].nodes() if isinstance(n,str)]
        temp_concs = copy.deepcopy(self.concs)
        del temp_concs['CO2']
        available = [choice(list(temp_concs.keys()),
                len(temp_concs),
                p=[x/sum(temp_concs.values()) for x in temp_concs.values()])[0]]
        source = random.choice([n for n in nodes if n in available]) # this needs to be better
        path = []
        for i in range(path_depth):
            target = random.choice(nodes)
            p = nx.shortest_path(self.graph[T][P],source,target)
            source = target
            if not len(p) == 1:
                path.append(p)
        return(path)   
    
    def generate_eqsystems_from_path(self,path,T,P):
        charged_species = {'CO3H':-1,'NH4':+1,'NH2CO2':-1} # this needs to be added to the arguments
        r_i = [i for i in path if isinstance(i,int)]
        rs = [self.reactions[T][P][i] for i in r_i]
        eqsystems = []
        for i in rs:
            r = i['e'].reac
            p = i['e'].prod
            k = i['k']
            substances = {}
            for n in list(it.chain(*[list(r)+list(p)])):
                if n in list(charged_species.keys()):
                    s = Substance.from_formula(n,charge=charged_species[n])
                    substances[s.name] = s
                    #substances[s.charge] = charged_species[n]
                else:
                    s = Substance.from_formula(n,charge=0)
                    substances[s.name] = s 
                    #substances[s.charge] = 0 
            #substances = {s.name:s for s in [Substance.from_formula(n) 
            #                                     for n in list(it.chain(*[list(r) + list(p)]))]}
            try:
                eqsystems.append(EqSystem([Equilibrium(r,p,k)],substances))
            except:
                eqsystems.append([])
        return(eqsystems[0])
    
    def equilibrium_concentrations_from_walk(self,eqwalk):
        fc = copy.deepcopy(self.concs)
        concs = {0:copy.deepcopy(fc)}
        successful_equations = []
        for i,eqsys in enumerate(eqwalk):
            fc = copy.deepcopy(fc)
            try:
                x,sol,sane = eqsys.root(fc)
                assert sol['success'] and sane
                for n,c in enumerate(x):
                    fc[eqsys.substance_names()[n]] = c
                    
                concs[i+1] = fc
                successful_equations.append(eqsys.string())
            except:
                concs[i+1] = fc
        return(concs,successful_equations)
    
    def graph_sampling_serial(self,sample_length=100,path_depth=50):
        with tqdm(total=len(self.trange)*len(self.prange),bar_format='{desc:<20}{percentage:3.0f}%|{bar:20}{r_bar}',position=0,leave=True) as pbar1:
            temperature_data = {}
            for T in self.trange:
                pressure_data = {}
                for P in self.prange:
                    samples = list(range(0,sample_length,1))
                    pbar1.set_description('T = {},P = {}'.format(T,P))
                    pressure_data[P] = self.sampling_function_no_queue(samples,T,P,path_depth)
                    pbar1.update(1)
                temperature_data[T] = pressure_data
        return(temperature_data)
    
    def sampling_function_no_queue(self,samples,T,P,path_depth):
        sample_data = {}
        with tqdm(total=len(samples),bar_format='{desc:<20}|{bar:10}|',position=1,leave=False) as pbar2:
            for sample in samples:
                pbar2.set_description('Sample {}'.format(sample))
                random_path = self.random_walk(T,P,path_depth)
                walk = [self.generate_eqsystems_from_path(step,T,P) for step in random_path]
                concs,equations = self.equilibrium_concentrations_from_walk(walk)
                sample_data[sample] = {'data':concs[len(concs)-1],'equation_statistics':equations}
                pbar2.update(1)
            pbar2.reset()
        return(sample_data)
    
    def sampling_function_queue(self,samples,T,P,path_depth,out_q):
        sample_data = {}
        with tqdm(total=len(samples),bar_format='{desc:<20}|{bar:10}|',position=1,leave=False) as pbar2:
            for sample in samples:
                pbar2.set_description('Samples'.format(sample))
                random_path = self.random_walk(T,P,path_depth) # here we predefine a random walk
                walk = [self.generate_eqsystems_from_path(step,T,P) for step in random_path] 
                concs,equations = self.equilibrium_concentrations_from_walk(walk)
                sample_data[sample] = {'data':concs[len(concs)-1],'equation_statistics':equations}
                pbar2.update(1)
            pbar2.reset()
        out_q.put(sample_data)
    
    def graph_sampling_processes(self,sample_length=100,path_depth=50,nprocs=4,random_path_depth=False):
        init_concs = copy.deepcopy(self.concs)
        with tqdm(total=len(self.trange)*len(self.prange),bar_format='{desc:<20}{percentage:3.0f}%|{bar:20}{r_bar}',position=0,leave=True) as pbar1:
            temperature_data = {}
            for T in self.trange:
                pressure_data = {}
                for P in self.prange:
                    pbar1.set_description('T = {},P = {}'.format(T,P))
                    result_dict = {0:{'data':init_concs,'equation_statistics':[]}}
                    out_queue = pmp.Queue()
                    samples = list(range(1,sample_length+1,1))
                    data_chunks = [samples[chunksize*i:chunksize*(i+1)] 
                            for i in range(nprocs) 
                            for chunksize in [int(math.ceil(len(samples)/float(nprocs)))]]
                    jobs = []
                    for chunk in data_chunks:
                        if random_path_depth == True:
                            path_depth = random.randint(1,100) #implementing a random function
                        process = pmp.Process(target=self.sampling_function_queue,
                                    args=(chunk,T,P,path_depth,out_queue))
                        jobs.append(process)
                        process.start()


                    for proc in jobs:
                        result_dict.update(out_queue.get())

                    for proc in jobs:
                        proc.join()

                    out_queue.close()

                    pressure_data[P] = result_dict
                    
                    pbar1.update(1)
                temperature_data[T] = pressure_data
        return(temperature_data)
            
