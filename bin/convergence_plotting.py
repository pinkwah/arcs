import pickle
import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd
import os,glob
from arcs.analysis_functions import get_mean_change_in_data
functional = 'hse06'

datafiles = {functional:pickle.load(open('convergence_data_{}.p'.format(functional),'rb'))} 

samples = list(datafiles[functional].keys())

mean = get_mean_change_in_data(datafiles[functional])

def PD_plotter(df,ax,xlabel,title):
    df.plot.line(ax=ax,marker='o',cmap='tab20',legend=False)
    ax.set_yscale('symlog')
    ax.set_xlabel(xlabel)
    ax.set_title(title)

def get_legend(df,ax):
    legend = plt.legend(frameon=False,ncol=5,fontsize='xx-small',loc=(0,-0.5))
    label_names = []
    charged_species = {'CO3H':'-','NH4':'+','NH2CO2':'-'}
    separater = ''
    for label in list(df.keys()):
        data = []
        for i,x in enumerate(label):
            try:
                y = int(x)
                data.append('$_{}$'.format(y))
            except:
                data.append(x)
        if label in charged_species:
            data.append(str(charged_species[label]))
        label_names.append(separater.join(data))

    for l in range(len(df.keys())):
        legend.get_texts()[l].set_text(label_names[l])
    return(legend)


fig,axes = plt.subplots(2,3,figsize=(10,10),dpi=300,sharey=True)
compounds = list(mean[samples[0]].keys())
for i,S in enumerate(samples):
    df = pd.DataFrame({comp:{PD:mean[S][PD][comp] 
                        for PD in mean.index} 
                  for comp in compounds})
    if i == 0:
        ax = axes[0][0]
    elif i == 1:
        ax = axes[0][1]
    elif i == 2:
        ax = axes[0][2]
    elif i == 3:
        ax = axes[1][0]
    elif i == 4:
        ax = axes[1][1]
    elif i == 5:
        ax = axes[1][2]
        
    PD_plotter(df,ax,'Path Depth',S)
    
get_legend(df,axes[0][0])

axes[0][0].set_yticklabels([])
axes[0][0].set_yticks([])
    
plt.tight_layout()
plt.savefig('{}_convergence_samples.png'.format(functional),transparent=True)
plt.savefig('{}_convergence_samples.pdf'.format(functional))

