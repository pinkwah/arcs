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


    def _screen_combinations_write(self,filename): #this is test - should be more memory efficient to write to an outfile then load later
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


# backup from official reactions     
import numpy as np
import itertools as it
import tqdm
from datetime import datetime
from monty.serialization import loadfn,dumpfn
import gzip,os,math
from copy import deepcopy
from pathos.helpers import mp as pmp 
from chempy import balance_stoichiometry,Reaction
from chempy.equilibria import Equilibrium
from chempy.reactionsystem import Substance

'''
1. Starts with a ReactionDictionaryGenerator - which generates a combination of numbers
2. Filters based on compounds present
3. Filters based on balanced reactions
'''


class ReactionsDictionaryGenerator:
    ''' a class that creates the initial reference dictionary for solving all permutations of reactions between N compounds ( in this case defaults to 30 compounds) 
    
    This should initially be run at the start as it takes a long time but acts as a reference dictionary for all compounds N length or lower '''
    
    def __init__(self,no_compounds=30,path='.'):
        self.nc = no_compounds
        self.path= os.path.abspath(path)
        

    def _reaction_run_and_clean_DEPRECATED(self,indexes=None,reactants_length=None,products_length=None,pbar=True):
        rs = tuple(it.combinations(indexes,reactants_length)) # can this be an iterator?
        ps = tuple(it.combinations(indexes,products_length))
        

        def _filter(rs,ps):
            filtered_list = []
            for i in rs:
                si = sorted(i)
                for j in ps:
                    sj = sorted(j)
                    if not si  ==  sj:
                        if not any([x in sj for x in si]):
                            if filtered_list:
                                if not tuple((si,sj)) or not tuple((sj,si)) in filtered_list:
                                    filtered_list.append(tuple((si,sj)))
                            else:
                                filtered_list.append(tuple((si,sj)))
            
            return(filtered_list)
        
        
    def _reaction_run_and_clean_speedup_pbar_DEPRECATED(self,indexes=None,reactants_length=None,products_length=None):
        rs = tuple(it.combinations(indexes,reactants_length)) # can this be an iterator?
        ps = tuple(it.combinations(indexes,products_length))
        
        filtered_list = []
        with tqdm.tqdm(total=len(rs)*len(ps),bar_format='{desc:<5.5}{percentage:3.0f}%|{bar:10}{r_bar}') as pbar:
            for i in rs:
                si = sorted(i)
                for j in ps:
                    sj = sorted(j)
                    pbar.update(1)
                    if not si  ==  sj:
                        if not any([x in sj for x in si]):
                            if filtered_list:
                                if not tuple((si,sj)) or not tuple((sj,si)) in filtered_list:
                                    filtered_list.append(tuple((si,sj)))
                            else:
                                filtered_list.append(tuple((si,sj)))
            
        return(filtered_list)
    
    
    def _runner_2_pbar_DEPRECATED(self,reaction_length=4,split_files=False):
        sizing = [x for x in it.combinations_with_replacement(np.arange(1,reaction_length),2) 
                  if np.sum(x) == reaction_length] 
        
        print('running_runner for {} component reactions...'.format(reaction_length))
        print('{} possibilities : {}'.format(len(sizing),sizing))
        
        if len(sizing) == 1 or split_files==False:
            for i,size in enumerate(sizing):
                if i == 0:
                    l = self._reaction_run_and_clean_speedup_pbar_DEPRECATED(indexes=np.arange(self.nc+1),
                                                    reactants_length=size[0],
                                                    products_length=size[1])
                else:
                    l.extend(self._reaction_run_and_clean_speedup_pbar_DEPRECATED(indexes=np.arange(self.nc+1),
                                                    reactants_length=size[0],
                                                    products_length=size[1]))
                filename=os.path.join(self.path,'data_{}.dat'.format(reaction_length))
                self.datawriter(data=l,name=filename)
                             
        else:
            for i,size in enumerate(sizing):
                l = self.reaction_run_and_clean_speedup_pbar_DEPRECATED(indexes=np.arange(self.nc+1),
                                                reactants_length=size[0],
                                                products_length=size[1])
                filename=os.path.join(self.path,'data_{}-{}.dat'.format(reaction_length,i))
                self.datawriter(data=l,name=filename)
        
        
        
    def reaction_filter_serial(self,tl): # for serial applications
        fl = []
        for k in tl:
            si = sorted(k[0])
            sj = sorted(k[1])
            if not si == sj:
                if not any([x in si for x in sj]):
                    if fl:
                        if not tuple((si,sj)) or not tuple((sj,si)) in fl:
                            fl.append(tuple((si,sj)))
                    else:
                        fl.append(tuple((si,sj)))
        return(fl)
    
    def reaction_filter_mp(self,tl,out_q): # for multiprocessing 
        fl = []
        for k in tl:
            si = sorted(k[0])
            sj = sorted(k[1])
            if not si == sj:
                if not any([x in si for x in sj]):
                    if fl:
                        if not tuple((si,sj)) or not tuple((sj,si)) in fl:
                            fl.append(tuple((si,sj)))
                    else:
                        fl.append(tuple((si,sj)))
        out_q.put(fl)
    

    def datawriter(self,data,name):
        with open(name,'w') as f:
            [f.write(str(i)+'\n') for i in data]
            
            
    def runner_serial(self,reaction_length=4,split_files=False):
        sizing = [x for x in it.combinations_with_replacement(np.arange(1,reaction_length),2) 
                  if np.sum(x) == reaction_length] 
        
        print('running_runner for {} component reactions...'.format(reaction_length),end=' ')
        print('{} possibilities : {}'.format(len(sizing),sizing))
        
        if len(sizing) == 1 or split_files==False:
            for i,size in enumerate(sizing):
                if i == 0:
                    rs = it.combinations([x for x in range(self.nc+1)],size[0])
                    ps = it.combinations([x for x in range(self.nc+1)],size[1])
                    tl = tuple(it.product(rs,ps))
                    print(size,':',len(tl),end='->')
                    l = self.reaction_filter_serial(tl)
                    print(len(l))
                else:
                    rs = it.combinations([x for x in range(self.nc+1)],size[0])
                    ps = it.combinations([x for x in range(self.nc+1)],size[1])
                    tl = tuple(it.product(rs,ps))
                    print(size,':',len(tl),end='->')
                    l2  = self.reaction_filter_serial(tl)
                    l.extend(l2)
                    print(len(l2))
                
                filename=os.path.join(self.path,'data_{}.dat'.format(reaction_length))
                self.datawriter(data=l,name=filename)
        else:
            for i,size in enumerate(sizing):
                rs = it.combinations([x for x in range(self.nc+1)],size[0])
                ps = it.combinations([x for x in range(self.nc+1)],size[1])
                tl = tuple(it.product(rs,ps))
                print(size,':',len(tl),end='->')
                l = self.reaction_filter_serial(tl)
                print(len(l))
                filename=os.path.join(self.path,'data_{}-{}.dat'.format(reaction_length,i))
                self.datawriter(data=l,name=filename)  
    
    def _mp_run(self,tl,nprocs):
        queue = pmp.Queue()
        chunksize = int(math.ceil(len(tl)/float(nprocs)))
        
        processes = []
        for i in range(nprocs):
            pr = pmp.Process(target=self.reaction_filter_mp,
                               args=(tl[chunksize*i:chunksize*(i+1)],queue)
                               )
            processes.append(pr)
            pr.start()
        
        data = []
        for i in range(nprocs):
            data.append(queue.get())
    
        for pr in processes:
            pr.join()    
        
        data_chain = tuple(it.chain(*data))
        # should the serial part come here?
        return(data_chain)
    
    
    def _while_procs(self,tl,nprocs):
        dat = []
        while nprocs>1:
            print(nprocs,end=':')
            nprocs = int(nprocs/2)
            tl = self._mp_run(tl,nprocs)
            if len(dat) == 3:
                if dat[-2] == dat[-1] == len(tl):
                    nprocs=1
            dat.append(len(tl))
            print(len(tl),end='->')
                
        return(tl)
    
            
        
        
                
    def runner_mp(self,reaction_length=4,split_files=False,nprocs=4):
        import math
        '''protect behind if __name__ == '__main__':'''
        sizing = [x for x in it.combinations_with_replacement(np.arange(1,reaction_length),2) 
                  if np.sum(x) == reaction_length] 
        
        print('running_runner for {} component reactions...'.format(reaction_length),end=' ')
        print('{} possibilities : {}'.format(len(sizing),sizing))
        
        if len(sizing) == 1 or split_files==False:
            for i,size in enumerate(sizing):
                if i == 0:
                    rs = it.combinations([x for x in range(self.nc+1)],size[0])
                    ps = it.combinations([x for x in range(self.nc+1)],size[1])
                    tl = tuple(it.product(rs,ps))
                    print(size,':',len(tl),end='->')
                    l = self._mp_run(tl,nprocs)
                    print(len(l))
                    #l_pre = self._mp_run(tl,nprocs)
                    #print(len(l_pre),end='->')
                    #l = self.reaction_filter_serial(l_pre)
                    #print(len(l))
                else:
                    rs = it.combinations([x for x in range(self.nc+1)],size[0])
                    ps = it.combinations([x for x in range(self.nc+1)],size[1])
                    tl = tuple(it.product(rs,ps))
                    print(size,':',len(tl),end='->')
                    l2 = self._mp_run(tl,nprocs)
                    l.extend(l2)
                    print(len(l2))
                    #l_pre = self._mp_run(tl,nprocs)
                    #print(len(l_pre),end='->')
                    #l2 = self.reaction_filter_serial(l_pre)
                    #print(len(l2))
                    #l.extend(l2)
                
                filename=os.path.join(self.path,'data_{}.dat'.format(reaction_length))
                self.datawriter(data=l,name=filename)
        else:
            for i,size in enumerate(sizing):
                rs = it.combinations([x for x in range(self.nc+1)],size[0])
                ps = it.combinations([x for x in range(self.nc+1)],size[1])
                tl = tuple(it.product(rs,ps))
                print(size,':',len(tl),end='->')
                l = self._mp_run(tl,nprocs)
                print(len(l))
               # l_pre = self._mp_run(tl,nprocs)
               #print(len(l_pre),end='->')
               # l = self.reaction_filter_serial(l_pre)
               # print(len(l))
                filename=os.path.join(self.path,'data_{}-{}.dat'.format(reaction_length,i))
                self.datawriter(data=l,name=filename)  
                
    def runner_mp_while(self,reaction_length=4,split_files=False,nprocs=4):
        import math
        '''protect behind if __name__ == '__main__':'''
        sizing = [x for x in it.combinations_with_replacement(np.arange(1,reaction_length),2) 
                  if np.sum(x) == reaction_length] 
        
        print('running_runner for {} component reactions...'.format(reaction_length),end=' ')
        print('{} possibilities : {}'.format(len(sizing),sizing))
        
        if len(sizing) == 1 or split_files==False:
            for i,size in enumerate(sizing):
                if i == 0:
                    rs = it.combinations([x for x in range(self.nc+1)],size[0])
                    ps = it.combinations([x for x in range(self.nc+1)],size[1])
                    tl = tuple(it.product(rs,ps))
                    print(size,':',len(tl),end='->')                    
                    l_pre = self._while_procs(tl,nprocs)
                    l = self.reaction_filter_serial(l_pre)
                    print(len(l))
                else:
                    rs = it.combinations([x for x in range(self.nc+1)],size[0])
                    ps = it.combinations([x for x in range(self.nc+1)],size[1])
                    tl = tuple(it.product(rs,ps))
                    print(size,':',len(tl),end='->')                    
                    l_pre = self._while_procs(tl,nprocs)
                    l2 = self.reaction_filter_serial(l_pre)
                    print(len(l2))
                    l.extend(l2)
                
                filename=os.path.join(self.path,'data_{}.dat'.format(reaction_length))
                self.datawriter(data=l,name=filename)
        else:
            for i,size in enumerate(sizing):
                rs = it.combinations([x for x in range(self.nc+1)],size[0])
                ps = it.combinations([x for x in range(self.nc+1)],size[1])
                tl = tuple(it.product(rs,ps))
                print(size,':',len(tl),end='->')
                l_pre = self._while_procs(tl,nprocs)
                l = self.reaction_filter_serial(l_pre)
                print(len(l))
                filename=os.path.join(self.path,'data_{}-{}.dat'.format(reaction_length,i))
                self.datawriter(data=l,name=filename)  
                    

    def clean_from_file(self,filename,nprocs=4):
        import gzip
        '''file is a gzip'''
        def _convert_line(file):
            for line in file:
                r,p = line.strip().split('],')
                r = tuple(int(x) for x in r.split('([')[1].split(',') if x)
                p = tuple(int(x) for x in p.split('[')[1].split('])')[0].split(',') if x)
                yield ((r,p))
        print('reading from file...',end=' ')
        f = tuple(_convert_line(gzip.open(filename,'rt')))
        print(len(f),end='->')
        l = self._mp_run(f,nprocs)
        print(len(l))
        filename = filename.replace('.dat.gz','_reloaded.dat')
        self.datawriter(data=l,name=filename)




class MappingtoReaction:
    ''' a class that takes a premade reactions reference dictionary and generates a list of reactions with it that are then further balanced and filtered'''
    
    def __init__(self,filename,compounds):
        self.filename = filename
        self.compounds = {i:c for i,c in enumerate(compounds)}
        
    def convert_file(self,file): # iterator that reads in the file (file is a gzip)
        for line in file:
            r,p = line.strip().split('],')
            r = tuple(int(x) for x in r.split('([')[1].split(',') if x)
            p = tuple(int(x) for x in p.split('[')[1].split('])')[0].split(',') if x)    
            yield ((r,p))
            
    def remove_indices(self,reaction_indexes):
        indexes =  tuple(self.compounds)
        approved = []
        for r in reaction_indexes:
            if not [x for x in it.chain(*r) if not x in indexes]:
                approved.append(r)
        return(approved)
    
    def convert_to_string(self,approved_list):
    
        def _screen_string(r):
            re,pr = r
            c = []
            for i in re:
                rs = [x for x in i]
                for j in rs:
                    try:
                        int(j)
                    except:
                        c.append(j)
        
            dr = set(dict.fromkeys(c))
        
            c = []
            for i in pr:
                ps = [x for x in i]
                for j in ps:
                    try:
                        int(j)
                    except:
                        c.append(j)
        
            dp = set(dict.fromkeys(c))
        
            if dr == dp:
                return(r)
    
        converted = []
        for r in approved_list:
            d = [[self.compounds[i] for i in r[0]],[self.compounds[i] for i in r[1]]]
            ds = _screen_string(d)
            if ds:
                converted.append(ds)
        return(converted)
    
    def convert_to_equation(self,converted_strings):
        converted = []
        for r in tqdm.tqdm(converted_strings):
            re,pr = r
            try:
                converted.append(balance_stoichiometry(list(re),list(pr),underdetermined=None))
            except:
                try:
                    converted.append(balance_stoichiometry(list(re),list(pr)))
                except:
                    pass                            
        return(converted) 
    
    def screen_converted(self,converted_reactions): # this should be multiprocessed - > perhaps a mp.Pool? as it only needs the big list
        def _convert_ord_to_dict(r):
            re,pr =  r
            try:
                reacs = {k:int(re[k]) for k in re}
                prods = {k:int(pr[k]) for k in pr}
            except:
                print('\n error with {}'.format(r))
            return(reacs,prods)
    
        screened = []
        for r in converted_reactions:
            try:
                re,pr = _convert_ord_to_dict(r)
                try:
                    screened.append(Equilibrium(re,pr))
                except:
                    pass
            except:
                pass
        return(screened)    
    
    def run_all(self):
        s = datetime.now()
        loi = tuple(self.convert_file(gzip.open(self.filename,'rt')))
        print('orig =',len(loi),end='...')
        approved = self.remove_indices(loi)
        print(' approved = ',len(approved),end='...')
        strings = self.convert_to_string(approved)
        print(' prescreening = ',len(strings))
        equations = self.convert_to_equation(strings)
        screened = self.screen_converted(equations)
        print(' final = ',len(screened),end='...')
        f = datetime.now() - s 
        print(' time = ',f)
        return(screened)
