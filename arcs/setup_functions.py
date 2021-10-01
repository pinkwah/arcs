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

def get_compound_directory(base,compound,size):
    return(os.path.join(base,compound,size))

class GetEnergyandVibrations:
    '''Class to get the Total Energy and Vibrations from a directory containing a calculations'''
    def __init__(self,directory):
        self.directory = directory
        
    def atoms(self):
        structure = read('{}/POSCAR.gz'.format(self.directory))
        return(structure)
        
    def energy(self):
        outcar = gzip.open('{}/OUTCAR.gz'.format(self.directory),'tr').readlines()
        for line in outcar:
            if 'y=' in line:
                energy = float(line.split()[-1])
        if len(list(dict.fromkeys(self.atoms().get_atomic_numbers()))) == 1:
            if not any(x in self.atoms().symbols.get_chemical_formula() for x in ['H2','O2','N2']):
                energy = energy / self.atoms().get_global_number_of_atoms()
    
        return(energy)
    
    def spin(self):
        outcar = gzip.open('{}/OUTCAR.gz'.format(self.directory),'tr').readlines()
        for line in outcar:
            if 'NELECT' in line:
                nelect = float(line.split()[2])
                return([0 if nelect %2 == 0 else 1][0])
            
    def pointgroup(self):
        atoms = self.atoms()
        pg = PointGroupAnalyzer(AseAtomsAdaptor.get_molecule(atoms)).get_pointgroup()
        return(pg.sch_symbol)
    
    def islinear(self):
        num_at = self.atoms()
        if num_at.get_global_number_of_atoms() == 1:
            return('monatomic')
        else:
            pg = self.pointgroup()
            if '*' in pg:
                return('linear')
            else:
                return('nonlinear')

    def rotation_num(self):
        pg = [x for x in self.pointgroup()]
        if pg[0] == 'C':
            if pg[-1] == 'i':
                rot = 1
            elif pg[-1] == 's':
                rot = 1
            elif pg[-1] == 'h':
                rot = int(pg[1])
            elif pg[-1] == 'v':
                if pg[1] == '*':
                    rot = 1
                else:
                    rot = int(pg[1])
            elif len(pg) == 2:
                rot = int(pg[-1])
                
        elif pg[0] == 'D':
            if pg[-1] == 'h':
                if pg[1] == '*':
                    rot = 2
                else:
                    rot = 2*int(pg[1])
            elif pg[-1] == 'd':
                rot = 2*int(pg[1])
            elif len(pg) == 2:
                rot = 2*int(pg[1])
            
        elif pg[0] == 'T':
            rot = 12
        
        elif pg[0] == 'O':
            rot = 24
        
        elif pg[0] == 'I':
            rot = 60        
                    
        return(rot)          
        
    def vibrations(self):
        new_directory = os.path.join(self.directory,'ibrion_6')
        outcar = gzip.open('{}/OUTCAR.gz'.format(new_directory),'tr').readlines()
        frequencies = []
        for line in outcar:
            if 'THz' in line:
                if not 'f/i' in line: # we ignore imaginary modes
                    ev = float(line.split()[-2]) / 1000
                    frequencies.append(ev)
        return(frequencies)      
    
    def as_dict(self):
        return({'atoms':self.atoms(),
                'spin':self.spin(),
                'rotation_num':self.rotation_num(),
                'islinear':self.islinear(),
                'energy':self.energy(),
                'vibrations':self.vibrations()})
    
    
class ReactionGibbsandEquilibrium:
    
    def __init__(self,reaction,temperature,pressure,reaction_input):
        self.reaction = reaction
        self.temperature = temperature
        self.pressure = pressure*100000 #pressure in bar
        self.reaction_input = reaction_input
        
    def Gibbs(self,c):
        data = self.reaction_input[c]
        igt = IdealGasThermo(vib_energies=data.vibrations(),
                                        geometry=data.islinear(),
                                        potentialenergy=data.energy(),
                                        atoms=data.atoms(),
                                        symmetrynumber=data.rotation_num(),
                                        spin=data.spin(),
                                        natoms=data.atoms().get_global_number_of_atoms())
        G = igt.get_gibbs_energy(self.temperature,self.pressure,verbose=False)
        H = igt.get_enthalpy(self.temperature,verbose=False)
        S = igt.get_entropy(self.temperature,self.pressure,verbose=False) 
        Z = igt.get_entropy(self.temperature,self.pressure,verbose=False)
        return({'G':G,'H':H,'S':S,'Z':Z})
    
    def reaction_energy(self):
        prod = self.reaction.prod
        reac = self.reaction.reac
        reaction_compounds = list(prod)+list(reac)
        # need to add a charge neutrality condition and mass balance violation
        gibbs = {c:self.Gibbs(c) for c in reaction_compounds}
        prod_sum = np.sum([self.Gibbs(c)['G']*prod[c] for c in gibbs if c in prod])
        reac_sum = np.sum([self.Gibbs(c)['G']*reac[c] for c in gibbs if c in reac])
        return(float(prod_sum - reac_sum))
    
    def equilibrium_constant(self):
        K = np.exp(-(self.reaction_energy()*e)/(Boltzmann*self.temperature))
        return(K)
    
    def as_dict(self):
        return({'G_react':self.reaction_energy(),'K_react':self.equilibrium_constant()})


class ReactionsGenerator:
    ''' once you have used this it is highly recommended that you  run class RemoveDuplicateReactions'''
    
    def __init__(self,compounds,max_length):
        self.compounds = compounds
        self.max_length = max_length
    

    def get_combinations(self,length):
        combinations = []
        init = list(it.permutations(self.compounds,length))
        for i in range(1,length):
            for j in init:
                combinations.append([j[0:i],j[i:len(j)]])
        return(combinations)
    
    def screen_combinations(self):
        screened = []
        for l in range(2,self.max_length+1):
            with tqdm(total=len(self.get_combinations(l)),bar_format='{desc:<5.5}{percentage:3.0f}%|{bar:10}{r_bar}') as pbar:
                for poss in self.get_combinations(l):
                    r,p = poss
                    pbar.update(1)
                    try:
                        screened.append(balance_stoichiometry(list(r),list(p)))
                    except:
                        pass
        return(screened)
       
    def get_reactions(self):
        reactions = []
        for i in self.screen_combinations():
            r,p = i
            try:
                #substances = list(r) + list(p)
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

class ApplyDataToReaction:
    ''' this class applies the gibbs data to a specific reaction'''
    
    def __init__(self,trange,prange,reactions,compound_data,nprocs):
        self.trange = trange
        self.prange = prange
        self.reactions = reactions
        self.compound_data = compound_data
        self.nprocs = nprocs
        
    def get_t_p_data(self,t,p): #serial
        reactions = {i:{'e':r,
            'k':ReactionGibbsandEquilibrium(r,t,p,self.compound_data).equilibrium_constant(),
            'g':ReactionGibbsandEquilibrium(r,t,p,self.compound_data).reaction_energy()} 
                     for i,r in tqdm(enumerate(self.reactions))}
        return(reactions)

    def get_t_p_data_mp(self,t,p): #multiprocessed

        manager = pmp.Manager()
        queue = manager.Queue()
        
        def mp_function(reaction_keys,out_q):

            data = {}
            for r in reaction_keys:
                rge = ReactionGibbsandEquilibrium(self.reactions[r],t,p,self.compound_data)
                data[r] = {'e':self.reactions[r],
                        'k':rge.equilibrium_constant(),
                        'g':rge.reaction_energy()}
            out_q.put(data)

        resultdict = {}
        r_keys = list(self.reactions.keys())
        chunksize = int(math.ceil(len(self.reactions)/float(self.nprocs)))
        processes = []

        for i in range(self.nprocs):
            pr = pmp.Process(target=mp_function,
                            args=(r_keys[chunksize*i:chunksize*(i+1)],queue))
            processes.append(pr)
            pr.start()

        for i in range(self.nprocs):
            resultdict.update(queue.get(timeout=1800))

        for pr in processes:
            pr.join()

        return(resultdict)
        
    def as_dict(self):
        data = {}
        with tqdm(total=len(self.trange)*len(self.prange),bar_format='{desc:<5.5}{percentage:3.0f}%|{bar:10}{r_bar}') as pbar:
            for t in self.trange:
                pdat = {}
                for p in self.prange:
                    pdat[p] = self.get_t_p_data_mp(t,p)
                    pbar.update(1)
                data[t] = pdat
        return(data) 

        
class GraphGenerator:
    
    def __init__(self,preloaded_data):
        self.preloaded_data = preloaded_data

    def _cost_function(self,gibbs,T,reactants):
        '''takes the cost function that is used in https://www.nature.com/articles/s41467-021-23339-x.pdf'''

        comps = []
        for r,n in reactants.items():
            for i in range(n):
                comps.append(r)

        num_atoms = np.sum([np.sum([y 
                             for x,y in Substance.from_formula(c).composition.items()]) 
                     for c in comps]) 

        return(np.log(1+(273/T)*np.exp(gibbs/num_atoms)))

    def multidigraph_cost(self,T,P):
        ''' this will weight the graph in terms of a cost function which makes it better for a Djikstra algorithm to work'''
        t = nx.MultiDiGraph(directed=True)
        for i,reac in self.preloaded_data[T][P].items():
            f_cost = self._cost_function(reac['g'],T,reac['e'].reac) #forward cost
            b_cost = self._cost_function(-reac['g'],T,reac['e'].prod) #backward cost
            r = list(reac['e'].reac)
            p = list(reac['e'].prod)
            t.add_weighted_edges_from([c,i,f_cost] for c in r) #reactants -> reaction
            t.add_weighted_edges_from([i,c,b_cost] for c in r) #reaction -> reactants
            t.add_weighted_edges_from([i,c,f_cost] for c in p) #reaction -> products
            t.add_weighted_edges_from([c,i,b_cost] for c in p) #products -> reaction
        return(t) #probably don't need something that has if and elif

    def multidigraph(self,T,P):
        t = nx.MultiDiGraph(directed=True)
        for i,reac in self.preloaded_data[T][P].items():
            r = list(reac['e'].reac)
            p = list(reac['e'].prod)
            k = reac['k'] # maybe check equilibrium.as_reactions ( gives forward and backward reactions!)
            if k <= 1: #favours reactants
                t.add_weighted_edges_from([c,i,1/k] for c in r)
                t.add_weighted_edges_from([i,c,k] for c in r)
                t.add_weighted_edges_from([i,c,1/k] for c in p)
                t.add_weighted_edges_from([c,i,k] for c in p)
            elif k >= 1: #favours products
                t.add_weighted_edges_from([c,i,1/k] for c in r)
                t.add_weighted_edges_from([i,c,k] for c in r)
                t.add_weighted_edges_from([i,c,1/k] for c in p)
                t.add_weighted_edges_from([c,i,k] for c in p)
        return(t) #need to check shortest pathways
    
    
    def multidigraph_from_t_and_p_range(self,trange,prange):
        graphs = {}
        with tqdm(total=len(trange)*len(prange),bar_format='{desc:<20}{percentage:3.0f}%|{bar:10}{r_bar}') as pbar:
            pbar.set_description('generating graph')
            for T in trange:
                pdict = {}
                for P in prange:
                    pdict[P] = self.multidigraph(T,P)
                    pbar.update(1)
                graphs[T] = pdict
        return(graphs)

    def multidigraph_from_t_and_p_range_cost(self,trange,prange):
        graphs = {}
        with tqdm(total=len(trange)*len(prange),bar_format='{desc:<20}{percentage:3.0f}%|{bar:10}{r_bar}') as pbar:
            pbar.set_description('generating graph')
            for T in trange:
                pdict = {}
                for P in prange:
                    pdict[P] = self.multidigraph_cost(T,P)
                    pbar.update(1)
                graphs[T] = pdict
        return(graphs)

class GenerateInitialConcentrations:

    def __init__(self,graph,T,P):
        self.graph = graph
        self.T = T
        self.P = P

    def all_random(self):
        compounds = [node for node in self.graph[self.T][self.P].nodes() if isinstance(node,str)]
        ic = {c:np.random.random()/1e6 for c in compounds}
        ic['CO2'] = 1
        return(ic)

    def specific_random(self,compounds=None):
        full_list = [n for n in self.graph[self.T][self.P].nodes() if isinstance(n,str)]
        ic = {}
        for c in full_list:
            if c  in self.compounds:
                ic[c] = np.random.random()/1e6
            else:
                ic[c] = 0 
        ic['CO2'] = 1
        return(ic)

    def from_file(self,file_name):
        nodes = [n for n in self.graph[self.T][self.P].nodes() if isinstance(n,str)]
        file_concentrations = loadfn(file_name)
        loaded_compounds = list(file_concentrations.keys())
        ic = {}
        for c in nodes:
            if c not in loaded_compounds:
                ic[c] = 0
            else:
                ic[c] = file_concentrations[c]
        ic['CO2'] = 1
        return(ic)




