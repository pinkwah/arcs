import os,gzip
from ase.io import read,write
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
from ase.thermochemistry import IdealGasThermo
from scipy.constants import Boltzmann, e
from monty.serialization import loadfn
import numpy as np 
import itertools as it
import numpy as np 
from chempy import balance_stoichiometry
from chempy import Reaction
from chempy.equilibria import Equilibrium
from chempy.reactionsystem import Substance
from tqdm import tqdm
import networkx as nx
from pathos.helpers import mp as pmp
import math
import copy
import warnings

class ReactionsGenerator:
    ''' once you have used this it is highly recommended that you  run class RemoveDuplicateReactions'''
    
    def __init__(self,compounds,max_length=5):
        self.compounds = compounds
        self.max_length = max_length
        if max_length < 2:
            raise ValueError('supplied max_length is too low. max_length should be >= 2')
    
    def _get_combinations(self,length):
        combinations = []
        init = list(it.permutations(self.compounds,length))
        for i in range(1,length):
            for j in init:
                combinations.append([j[0:i],j[i:len(j)]])
        return(combinations)
    
    def _screen_combinations(self):
        screened = []
        for l in range(2,self.max_length+1):
            with tqdm(total=len(self._get_combinations(l)),bar_format='{desc:<5.5}{percentage:3.0f}%|{bar:10}{r_bar}') as pbar:
                for poss in self._get_combinations(l):
                    r,p = poss
                    pbar.set_description('{}'.format(l))
                    pbar.update(1)
                    try:
                        screened.append(balance_stoichiometry(list(r),list(p)))
                    except:
                        pass
        return(screened)
       
    def get_reactions(self):
        reactions = []
        for i in self._screen_combinations():
            r,p = i
            try:
                #substances = list(r) + list(p)
                reactions.append(Equilibrium(r,p))  
            except:
                pass
        return(reactions)


    def _screen_combinations_write(self,filename): # this is test - should be more memory efficient to write to an outfile then load later
        with open('{}'.format(filename),'w') as outfile:
            for l in range(2,self.max_length+1):
                with tqdm(total=len(self._get_combinations(l)),bar_format='{desc:<5.5}{percentage:3.0f}%|{bar:10}{r_bar}') as pbar:
                    for poss in self._get_combinations(l):
                        r,p = poss
                        pbar.set_description('{}'.format(l))
                        pbar.update(1)
                        try:
                            re,pr = balance_stoichiometry(list(r),list(p))
                            re = Reaction(re,pr)
                            outfile.write('\n')
                            outfile.write(re.string())
                            outfile.flush()
                        except:
                            pass
        outfile.close()

    def get_reactions_from_file(self,filename):
        reactions = []
        self._screen_combinations_write(filename)
        infile = open(filename,'r').readlines()
        for i in self._screen_combinations_write():
            r,p = i
            try:
                reactions.append(Equilibrium(r,p))  
            except:
                pass
        return(reactions)
    
class RemoveDuplicateReactions: # this needs checking
    def __init__(self,reactions,divisor):
        self.reactions = copy.deepcopy(reactions)
        self.divisor = divisor
    
    def check_same(self,i,j):
        r1 = list(self.reactions[i].reac)
        p1 = list(self.reactions[i].prod)
    
        r2 = list(self.reactions[j].reac)
        p2 = list(self.reactions[j].prod)
    
        action = 0
    
        if all(x in r2 for x in r1):
            if all(y in p2 for y in p1):
                action = 1
        elif all(x in r2 for x in p1):
            if all(y in p2 for y in r1):
                action = 1
        
        return(action)

    def run_through_permutations(self,l):
        for i,j in l:
            act = self.check_same(i,j)
            if act == 1:
                self.reactions.pop(j)
                return(0)
                break
            
    def iterate(self):
        value = 0
        while value == 0:
            l = list(it.combinations(list(range(len(self.reactions))),2))
            value = self.run_through_permutations(l)
        return(self.reactions)
    
    def clean(self,split):
        splitted_reactions = np.array_split(list(it.chain(*[self.reactions])),split)
        cleaned = []
        for r_split in splitted_reactions:
            self.reactions = copy.deepcopy(r_split.tolist()) # this should probably be multiprocessed
            d = self.iterate()
            cleaned.append(d)
        return(list(it.chain(*cleaned)))
    
    def whileclean(self):
        for div in tqdm(range(self.divisor,0,-1)):        
            finished = False
            while finished == False:
                init = len(list(it.chain(*[self.reactions])))
                self.reactions = self.clean(div)
                final = len(list(it.chain(*[self.reactions])))
                if init == final:
                    finished = True
        return(self.reactions)
