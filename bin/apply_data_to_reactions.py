#!/usr/bin/env python3.8
import pickle
from arcs.setup_functions import *
import warnings
import os,glob

if __name__ == '__main__':
    functional = 'pbesol'
    trange = [250,300,325] #range(200,400,40)
    prange = [20,90,98,100] #range(7,300,40)
    dft_directory = '../../dft_data/{}/'.format(functional)
    compounds = [x.split('/')[-2] for x in  glob.glob('../../dft_data/{}/*/'.format(functional))]
    box_size = 10
    dft_scraped_data = {c:GetEnergyandVibrations(get_compound_directory(dft_directory,c,str(box_size))) for c in compounds}
    # this needs to be updated as in the reactions data this is 'S' and not 'S8'
    dft_scraped_data['S'] = dft_scraped_data.pop('S8')
    warnings.simplefilter('ignore')
    reactions = pickle.load(open('../../reaction_data/new_unique_reactions_2_100.p','rb'))
    nprocs = 32

    data = ApplyDataToReaction(trange,prange,reactions,dft_scraped_data,nprocs).as_dict()
    
    pickle.dump(data,open('{}_reactions.p'.format(functional),'wb'))
