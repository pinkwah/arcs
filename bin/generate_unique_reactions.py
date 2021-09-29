import pickle
import os,glob
import numpy as np 
import itertools as it
from arcs.setup_functions import ReactionsGenerator,RemoveDuplicateReactions

compounds = ['H2O',
 'H2S',
 'CH4N2O',
 'C2H6',
 'C3H6O2',
 'HNO2',
 'C2H5NO',
 'CH6N2O2',
 'H2',
 'H2SO4',
 'H2SO3',
 'HNO3',
 'NO',
 'C2H4O',
 'SO3',
 'O2',
 'SO2',
 'CO',
 'CO3H',
 'CH4O',
 'C4H10',
 'NH5CO3',
 'C4H8',
 'HSNO5',
 'CH4',
 'C2H4',
 'H2O2',
 'C2H4O2',
 'NH2CO2',
 'NO2',
 'N2',
 'C2H6S',
 'C3H6',
 'NH3',
 'NH4',
 'C2H6O',
 'CH2O',
 'CO2',
 'S']
print(compounds) 
try:
    reactions = pickle.load(open('uncleaned_reactions.p','rb'))
    print('reactions loaded!')
except:
    print('starting reaction generator...')
    reaction_generator = ReactionsGenerator(compounds=compounds,max_length=5)
    reactions = reaction_generator.get_reactions()
    pickle.dump(reactions,open('uncleaned_reactions.p','wb'))

print('total number of reactions = {}'.format(len(reactions)))
print('removing duplicate reactions...')
cl = RemoveDuplicateReactions(reactions,1000)
cleaned = cl.whileclean()
pickle.dump({i:j for i,j in enumerate(cleaned)},open('new_unique_reactions_2.p','wb'))
