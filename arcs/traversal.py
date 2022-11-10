import random
from chempy.equilibria import Equilibrium,EqSystem
from chempy import Substance
import copy
import networkx as nx
import itertools as it
from tqdm.notebook import tqdm
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
import warnings
import time

class Traversal:
    def __init__(self,graph,reactions,concs,trange,prange,co2=False,**kwargs):
        self.graph = graph
        self.reactions = reactions
        self.concs = copy.deepcopy(concs)
        #self.progress = kwargs['progress']
        self.trange = trange
        self.prange = prange
        if co2==False:
            self.concs['CO2'] = 0.0
        pass
    
#    ######################################################################################################################################################        
    def _get_weighted_random_compounds(self,T,P,concs=None,co2=False,max_compounds=2,probability_threshold=0.05):     # doesnt' include looping
        # 1. choose a compound
        # 2. choose another compound
        
        nodes = [n for n in self.graph[T][P].nodes() if isinstance(n,str)]
        temp_concs = copy.deepcopy(concs)            

        #if not force_selection == None:
        #    new_temp_concs = {}
        #    force_selection = force_selection + ['CO2']
        #    for c in temp_concs:
        #        if c in force_selection:
        #            new_temp_concs[c] = temp_concs[c]

        #    temp_concs = new_temp_concs

        if co2 == False:
            del temp_concs['CO2']
            
        probabilities = {k:v/sum(temp_concs.values()) for k,v in temp_concs.items()}
        available = [choice(list(temp_concs.keys()),
                len(temp_concs),
                p=list(probabilities.values()))][0]
        num_species = len([x for x in temp_concs.values() if x > 0])
        if max_compounds > num_species:
            #print('Warning: max_compounds {} > {} -> now {}'.format(max_compounds,num_species,num_species)) 
            max_compounds = num_species
            
            
        choices = {}
        for c in range(max_compounds): 
            if c == 0:
                c1 = random.choice([n for n in nodes if n in available])
                if probabilities[c1] >= probability_threshold:
                    choices[c1] = probabilities[c1]
            else:
                #c2 = random.choice([x for x in available if not x in choices])
                c2 = random.choice([x for x in available])

                if probabilities[c2] >= probability_threshold:
                    choices[c2] = probabilities[c2]
        return(choices)
    
    def _get_weighted_reaction_rankings(self,T,P,choices,max_rank=5,method='Bellman-Ford'):
        rankings = {}
        if len(choices) > 1:
            possibilities = list(nx.shortest_paths.all_shortest_paths(self.graph[T][P],list(choices)[0],list(choices)[1],method=method)) # Bellman-Ford
            for x in possibilities:
                candidates = list(self.graph[T][P][x[1]])
                if len(choices) > 2:
                    for c in list(choices)[2:]:
                        if c in candidates:
                            weight = self.graph[T][P].get_edge_data(x[0],x[1])[0]['weight']
                            rankings[x[1]] = {'candidates':candidates,'weight':weight}
                else:
                    weight = self.graph[T][P].get_edge_data(x[0],x[1])[0]['weight']
                    rankings[x[1]] = {'candidates':candidates,'weight':weight}
        if rankings:
            topranks = [x for i,x in enumerate(rankings) if i<=max_rank]
            rankings = {x:rankings[x] for x in topranks}
            return(pd.DataFrame(rankings).sort_values(by='weight',axis=1).to_dict()) # sort them according to lowest weight first
        else:
            return(None)
        

    def _random_choice_unconnected(self,T,P,force_direct=False,co2=False): # currently randomly disjointed reactions that are weighted
        nodes = [n for n in self.graph[T][P].nodes() if isinstance(n,str)]
        if force_direct == True:
            pstring = [0,1,2]
            while len(pstring) > 2:
                source = self._get_weighted_random_compound(T,P,co2=co2,force_selection=None) 
                target = random.choice(nodes)
                p = nx.shortest_path(self.graph[T][P],source,target,weight='weight')
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
                p = nx.shortest_path(self.graph[T][P],source,target,weight='weight')
                pstring = [n for n in p if isinstance(p,str)]
        else:
                source = self._get_weighted_random_compound(T,P)
                target = random.choice(nodes)
                p = nx.shortest_path(self.graph[T][P],source,target)
        return(p)
    
    def generate_eqsystem(self,index,T,P):
        charged_species = {'CO3H':-1,'NH4':+1,'NH2CO2':-1} # this needs to be added to the arguments
        rs = self.reactions[T][P][index]
        r = rs['e'].reac
        p = rs['e'].prod
        k = rs['k']
        substances = {}
        for n in list(it.chain(*[list(r)+list(p)])):
            if n in list(charged_species.keys()):
                s = Substance.from_formula(n,charge=charged_species[n])
                substances[s.name] = s
            else:
                s = Substance.from_formula(n,charge=0)
                substances[s.name] = s 
        try:
            return(EqSystem([Equilibrium(r,p,k)],substances)) # might not just be able to try a return...
        except:
            return(None)
        
        
    def equilibrium_concentrations(self,concs,eq):
        fc = copy.deepcopy(concs)
        try:
            x,sol,sane = eq.root(fc)
            assert sol['success'] and sane
            for n,c in enumerate(x):
                fc[eq.substance_names()[n]] = c
                    
            concs = fc
            eq = eq.string()
        except:
            concs = fc
            eq = None
        return(concs,eq)
    
    def random_path_scenario3(self,T,P,
                              probability_threshold=0.05,
                              path_depth=50,
                              #concs=None,
                              max_compounds=5,
                              max_rank=5):
    
        final_concs = {0:copy.deepcopy(self.concs)} 
        reactionstats = {0:None}

        for ip in range(1,path_depth+1):
            fcs = copy.deepcopy(final_concs[ip-1])
            try:
                choices = self._get_weighted_random_compounds(T,P,concs=fcs,
                                               max_compounds=max_compounds,
                                               probability_threshold=probability_threshold)
            except:
                path_depth = ip+1
                break
            if len(choices) <= 1: # not sure this is necessary....
                path_depth = ip+1
                break
        
            rankings = self._get_weighted_reaction_rankings(T,P,choices,max_rank=max_rank,method='bellman-ford')
        
            if not rankings:
                break
        
            weights = {k:1/rankings[k]['weight'] for k in rankings}
            probabilities = {k:v/sum(weights.values()) for k,v in weights.items()}
            chosen_reaction = random.choice([choice(list(probabilities.keys()),
                                len(probabilities),
                                p=list(probabilities.values()))][0])

            eqsyst = self.generate_eqsystem(chosen_reaction,T,P)
            # if reaction was previous reaction then break
            path_available = [r for r in reactionstats.values() if not r==None]
            if path_available:
                if eqsyst.string() == path_available[-1] and eqsyst.string() == path_available[-1]:
                    break   # extra break
                    
            final_concs[ip],reactionstats[ip] = self.equilibrium_concentrations(fcs,eqsyst)
        return({'data':final_concs[list(final_concs)[-1]],
                'equation_statistics':[r for r in reactionstats.values() if not r==None],
                'path_length':len([r for r in reactionstats.values() if not r==None])})  
