import random
from chempy.equilibria import Equilibrium,EqSystem
from chempy import Substance
import copy
import networkx as nx
import itertools as it
#from tqdm.notebook import tqdm
from tqdm import tqdm
from collections import defaultdict

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
        
    
    def random_walk(self,T,P,path_depth=10):
        nodes = [n for n in self.graph[T][P].nodes() if isinstance(n,str)]
        available = [a for a,i in self.concs.items() if i !=0]
        source = random.choice([n for n in nodes if n in available])
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
                random_path = self.random_walk(T,P,path_depth)
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
    
#    def graph_sampling_pool_apply_async(self,sample_length=100,path_depth=50,nprocs=4):
#        with tqdm(total=len(self.trange)*len(self.prange),bar_format='{desc:<20}{percentage:3.0f}%|{bar:20}{r_bar}',position=0,leave=True) as pbar1:
#            temperature_data = {}
#            for T in self.trange:
#                pressure_data = {}
#                for P in self.prange:
#                    pbar1.set_description('T = {},P = {}'.format(T,P))
#                    samples = list(range(1,sample_length+1,1))
#                    data_chunks = [samples[chunksize*i:chunksize*(i+1)] 
#                            for i in range(nprocs) 
#                            for chunksize in [int(math.ceil(len(samples)/float(nprocs)))]]
#
#                    pool =  multiprocessing.Pool(nprocs)
#                    jobs = []
#                    for chunk in data_chunks:
#                        jobs.append(pool.apply_async(self.sampling_function_no_queue,args=(chunk,T,P,path_depth)))
#
#                    pressure_data[P] = [result.get() for result in jobs]
#                    pool.close()
#                    pbar1.update(1)
#                temperature_data[T] = pressure_data
#        return(temperature_data)
#    
#    def graph_sampling_pool_imap(self,sample_length=100,path_depth=50,nprocs=4):
#        with tqdm(total=len(self.trange)*len(self.prange),bar_format='{desc:<20}{percentage:3.0f}%|{bar:20}{r_bar}',position=0,leave=True) as pbar1:
#            temperature_data = {}
#            for T in self.trange:
#                pressure_data = {}
#                for P in self.prange:
#                    samples = list(range(0,sample_length,1))
#                    pool = multiprocessing.Pool(nprocs)
#                    pbar1.set_description('T = {},P = {}'.format(T,P))
#                    pressure_data[P]  = pool.imap(self.sampling_function_no_queue,(samples,T,P,path_depth))
#                    pbar1.update(1)
#                temperature_data[T] = pressure_data
#        return(temperature_data)
    
    
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

    dict_of_dataframes = {T:{P:get_dataframes(equations[T][P]) 
                             for P in prange} 
                          for T in trange}
    return(dict_of_dataframes)

def get_mean_change_in_data(t_and_p_data,percentage=True):
    if not percentage == True:
        data = {T:
                {P:
                  {x:t_and_p_data[T][P][x]['data']
                   for x in t_and_p_data[T][P].keys()} 
                 for P in t_and_p_data[T].keys()} 
                for T in t_and_p_data.keys()}
        d = pd.DataFrame(data) 
        mean_dataframe = pd.DataFrame({T:
                                       {P:pd.DataFrame(d[T][P]).T.mean().drop('CO2') - pd.DataFrame(d[T][P])[0].drop('CO2')
                                        for P in d.index} 
                                       for T in d.columns})
    else:
        data = {T:
                {P:
                  {x:t_and_p_data[T][P][x]['data']
                   for x in t_and_p_data[T][P].keys()}
                 for P in t_and_p_data[T].keys()}
                for T in t_and_p_data.keys()}
        d = pd.DataFrame(data)
        mean_dataframe = pd.DataFrame({T:
                                       {P:((pd.DataFrame(d[T][P]).T.mean().drop('CO2') - pd.DataFrame(d[T][P])[0].drop('CO2')) / pd.DataFrame(d[T][P])[0].drop('CO2'))*100
                                        for P in d.index}
                                       for T in d.columns})
    return(mean_dataframe)
    

class PrettyPlot:
    
    def __init__(self,graph,concs,path,index,directory):
        self.graph = graph
        self.concs = concs
        self.path = path
        self.index = index
        self.directory = directory
        
        
    def make_options(self,p):
        node_sizes = [1 if isinstance(n,str) else 0 for n in self.graph.nodes()]
        node_colours = []
        alphas = []
        for n in self.graph.nodes:
             if isinstance(n,str):
                    if n in list(self.concs.keys()):
                        node_colours.append((0.8,0.0,0.0))
                        alphas.append(1.0)
                    else:
                        node_colours.append((0.6,0.4,0.9))
                        alphas.append(0.5)
             else:
                 node_colours.append((0.0,0.4,0.8))
                 alphas.append(0.2)
            
        node_options = {'node_color': node_colours,
                   'alpha': alphas,
                   'node_size': node_sizes}
        
        edge_colours = []
        for e in self.graph.edges:
            if p[0] == e[0] and p[1] == e[1]:
                edge_colours.append((0.9,0.0,0.0,1.0))
            else:
                edge_colours.append((0.2,np.random.random(),np.random.random(),0.1))
        
        edge_options = {'connectionstyle':'arc3,rad=0.9',
                        'width':1,
                       'edge_color':edge_colours}
        
        return([node_options,edge_options])
    
    def make_directory(self):
        np = os.path.join(self.directory,str(self.index))
        try:
            os.mkdir(np)  
        except:
            pass
        return(np)
        
    def plot_graph(self):
        np = self.make_directory()
        for i,pt in enumerate(self.path):
            node_options,edge_options = self.make_options(pt)        
            fig,ax = plt.subplots(figsize=(20,20),dpi=100)
            pos = nx.kamada_kawai_layout(self.graph)
            n = nx.draw_networkx_nodes(self.graph,pos,ax=ax,
                                       node_color=node_options['node_color'],
                                       alpha=node_options['alpha'],
                                      node_size=node_options['node_size'])
            e = nx.draw_networkx_edges(self.graph,pos,ax=ax,**edge_options)
            ax.set_facecolor((0.1, 0.1, 0.1))
            f = os.path.join(np,'nl_{}_{}.png'.format(self.index,i))
            plt.savefig(f,transparent=False)        
            
