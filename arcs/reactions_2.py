import numpy as np
import itertools as it
import tqdm
from datetime import datetime
from monty.serialization import loadfn,dumpfn
import gzip,os,math
from copy import deepcopy
from pathos.helpers import mp as pmp 

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
            data.append(queue.get(timeout=3600))
    
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
            dat.append(tl)
            print(len(tl),end='->')
            if dat[-2] == dat[-1]:
                break
                
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
                    l_pre = self._mp_run(tl,nprocs)
                    print(len(l_pre),end='->')
                    l = self.reaction_filter_serial(l_pre)
                    print(len(l))
                else:
                    rs = it.combinations([x for x in range(self.nc+1)],size[0])
                    ps = it.combinations([x for x in range(self.nc+1)],size[1])
                    tl = tuple(it.product(rs,ps))
                    print(size,':',len(tl),end='->')
                    l_pre = self._mp_run(tl,nprocs)
                    print(len(l_pre),end='->')
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
                l_pre = self._mp_run(tl,nprocs)
                print(len(l_pre),end='->')
                l = self.reaction_filter_serial(l_pre)
                print(len(l))
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
                l_pre = self._while_procs(tl,nprocs)
                l = self.reaction_filter_serial(l_pre)
                print(len(l))
                filename=os.path.join(self.path,'data_{}-{}.dat'.format(reaction_length,i))
                self.datawriter(data=l,name=filename)  
                    
                
###blach                

                
                

#class ReactionsDictionaryToEquation:
#    ''' a class that takes a premade reactions reference dictionary and generates a list of reactions with it that are then further balanced and filtered'''
#    
#    def __init__(path,file