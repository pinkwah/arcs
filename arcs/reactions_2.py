import numpy as np
import itertools as it
import tqdm
from datetime import datetime
from monty.serialization import loadfn,dumpfn
import gzip,os
from copy import deepcopy

class ReactionsDictionaryGenerator:
    ''' a class that creates the initial reference dictionary for solving all permutations of reactions between N compounds ( in this case defaults to 30 compounds) 
    
    This should initially be run at the start as it takes a long time but acts as a reference dictionary for all compounds N length or lower '''
    
    def __init__(self,no_compounds=30,path='.'):
        self.nc = no_compounds
        self.path= os.path.abspath(path)
        
        
    def reaction_run_and_clean(self,indexes=None,reactants_length=None,products_length=None):
        rs = list(it.combinations(indexes,reactants_length)) # can this be an iterator?
        ps = list(it.combinations(indexes,products_length))
        
        filtered_list = []
        with tqdm.tqdm(total=len(rs)*len(ps),bar_format='{desc:<5.5}{percentage:3.0f}%|{bar:10}{r_bar}') as pbar:
             for i in rs:
                for j in ps:
                    pbar.update(1)
                    if not set(j)  ==  set(i):
                        if not any([x in j for x in i]):
                            if filtered_list:
                                if not (tuple(set(i)),tuple(set(j))) in filtered_list:
                                    if not (tuple(set(j)),tuple(set(i))) in filtered_list:
                                        filtered_list.append((tuple(set(i)),tuple(set((j)))))
                            else:
                                filtered_list.append((tuple(set(i)),tuple(set((j)))))
            
        return(filtered_list)


    def datawriter(self,data,name):
        with open(name,'w') as f:
            [f.write(str(i)+'\n') for i in data]
            
    
    def runner(self,reaction_length=4,split_files=False):
        sizing = [x for x in it.combinations_with_replacement(np.arange(1,reaction_length),2) 
                  if np.sum(x) == reaction_length] 
        
        print('running_runner for {} component reactions...'.format(reaction_length))
        print('{} possibilities : {}'.format(len(sizing),sizing))
        
        if len(sizing) == 1 or split_files==False:
            for i,size in enumerate(sizing):
                if i == 0:
                    l = self.reaction_run_and_clean(indexes=list(np.arange(self.nc+1)),
                                                    reactants_length=size[0],
                                                    products_length=size[1])
                else:
                    l.extend(self.reaction_run_and_clean(indexes=list(np.arange(self.nc+1)),
                                                    reactants_length=size[0],
                                                    products_length=size[1]))
                filename=os.path.join(self.path,'data_{}.dat'.format(reaction_length))
                self.datawriter(data=l,name=filename)
                             
        else:
            for i,size in enumerate(sizing):
                l = self.reaction_run_and_clean(indexes=list(np.arange(self.nc+1)),
                                                reactants_length=size[0],
                                                products_length=size[1])
                filename=os.path.join(self.path,'data_{}-{}.dat'.format(reaction_length,i))
                self.datawriter(data=l,name=filename)
                
                

class ReactionsDictionaryToEquation:
    ''' a class that takes a premade reactions reference dictionary and generates a list of reactions with it that are then further balanced and filtered'''
    
    def __init__(path,file