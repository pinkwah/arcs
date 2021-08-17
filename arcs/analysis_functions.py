import random
from chempy.equilibria import Equilibrium,EqSystem
from chempy import Substance
import copy
import networkx as nx
import itertools as it
from tqdm import tqdm
from collections import defaultdict


import os
import warnings

from multiprocessing import Process, Queue
from datetime import datetime
import random
import math
import threading
import multiprocessing
import numpy as np
import pandas as pd 
from pathos.helpers import mp as pmp


class Traversal:
    def __init__(self,graph,reactions,concs,**kwargs):
        self.graph = graph
        self.reactions = reactions
        self.concs = copy.deepcopy(concs)
        self.progress = kwargs['progress']
        
    
    def random_walk(self,T,P,depth=10):
        nodes = [n for n in self.graph[T][P].nodes() if isinstance(n,str)]
        available = [a for a,i in self.concs.items() if i !=0]
        source = random.choice([n for n in nodes if n in available])
        path = []
        for i in range(depth):
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
                successful_equations.append(eqsys)
            except:
                concs[i+1] = fc
        return(concs,successful_equations)
    
    def sample_graph_serial(self,T,P,sample_length=100,path_depth=100): # to be used when there are problems with sample_graph_mp
        total_data = {0:self.concs}
        s_equs = []
        with tqdm(total=sample_length*path_depth,disable=self.progress) as pbar:
            for i in range(1,sample_length,1):
                p = self.random_walk(T,P,path_depth)
                walk = [] 
                for path in p:
                    walk.append(self.generate_eqsystems_from_path(path,T,P))
                    pbar.update(1)    
                da,eq = self.equilibrium_concentrations_from_walk(walk)
                total_data[i] = da[len(da)-1]
                s_equs.append(eq)
        return(total_data,s_equs)
    
    def sample_graph_mp(self,T,P,sample_length=100,path_depth=100,nprocs=4):
        
        def mp_function(samples,out_q):
            data = {}
            for sample in samples:
                p = self.random_walk(T,P,path_depth)
                walk = []
                for path in p:
                    walk.append(self.generate_eqsystems_from_path(path,T,P))
                da,eq = self.equilibrium_concentrations_from_walk(walk)
                data[sample] = {'data':da[len(da)-1],'equation_statistics':eq}
            out_q.put(data)
            
            
        resultdict = {0:{'data':self.concs,'equation_statistics':[]}}
        #manager = pmp.Manager()
        queue = pmp.Queue()
        total_samples = range(1,sample_length,1)
        chunksize = int(math.ceil(len(total_samples)/float(nprocs)))
        procs = []
        for i in range(nprocs):
            p = pmp.Process(target=mp_function,
                        args=(total_samples[chunksize*i:chunksize*(i+1)],queue))
            procs.append(p)
            p.start() 
    
        for i in range(nprocs):
            #resultdict.update(queue.get(timeout=1800))
<<<<<<< HEAD
            try:
                resultdict.update(queue.get(block=True, timeout=1800))
            except queue.empty():
                continue
            #resultdict.update(queue.get())
=======
            resultdict.update(queue.get())
>>>>>>> 216b2c73a86f669e764f37a3008faf0726169ec5
            
        queue.close()
            
        for p in procs:
            p.join()
            
<<<<<<< HEAD
        queue.close()
=======
        #queue.join_thread()
>>>>>>> 216b2c73a86f669e764f37a3008faf0726169ec5
            
            
        return(resultdict)
    
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
                appearances[i.string()] += 1
    
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

def get_mean_change_in_data(t_and_p_data):
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
            
