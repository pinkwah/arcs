import warnings
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from datetime import datetime
import pickle
import gzip
import math

import numpy as np 
from scipy.interpolate import interpn
import pandas as pd
from tqdm import tqdm
import matplotlib.colors as colors
import os
from monty.serialization import loadfn,dumpfn
import matplotlib




def plot_difference(axis,mean_data,temperature,pressure,threshold,ymin,ymax,mean_percentages=None):
    ax=axis
    T = temperature  
    P = pressure
    
    df_t = mean_data
    
    norm = colors.Normalize(vmin = ymin, vmax = ymax)
    cmap = matplotlib.cm.RdBu
    colours = [cmap(norm(y)) for x,y in df_t[T][P].items()]
    
    df_t[T][P].plot(kind='bar',ax=ax,color=colours)

    ax.set_yscale('symlog')
    label_names = []
    charged_species = {'CO3H':'-','NH4':'+','NH2CO2':'-'}
    separater = ''
    for label in list(df_t[T][P].keys()):
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
    
    ax.set_ylabel('$\Delta$ Concentrations (ppm/Arb. units)')    
    ax.set_xticklabels(label_names,rotation=80)
    ax.set_ylim(ymin,ymax)
    ax.set_yticks(np.linspace(ymin,ymax,5))
    ppm_labels = ['{:.1F}'.format(x/1e-6) for x in ax.get_yticks()]
    ax.set_yticklabels(ppm_labels)
    ax.axhline(0,color='k')
    rects = ax.patches
    percentages = []
    if isinstance(mean_percentages,pd.DataFrame):

        for comp,perc in mean_percentages[T][P].items():
            if not perc == np.inf:
                if not  math.isnan(perc) == True:
                    percentages.append('{:.0F}'.format(perc))
                else:
                    percentages.append(None)
            else:
                percentages.append(None)
    new_label_names = []
    if not len(percentages) == 0:
        for i in range(len(label_names)):
            if not percentages[i] == None:
                new_label_names.append('{}\n{}%'.format(label_names[i],percentages[i]))
            else:
                new_label_names.append(label_names[i])
    label_names = new_label_names
            
    for rect, label in zip(rects,label_names):
        y = rect.get_height()
        x = rect.get_x() + rect.get_width() / 2                
            
        space = 5
        va = 'bottom'
        if y < 0 :
            space *= -1
            va = 'top'
        if y > threshold or y < -threshold:    
            ax.annotate(
                label,
                (x,y),
                xytext=(0,space),
                textcoords='offset points',
                ha='center',
                va=va,
                rotation=60)
    plt.grid(axis = 'y',alpha=0.5)    

    threshold_rectangle = patches.Rectangle(
        (0-rects[0].get_width(),-threshold),len(label_names),threshold*2,alpha=0.5,color='grey',edgecolor=None
    )
    ax.add_patch(threshold_rectangle)
    legend_patch = patches.Patch(color='grey',alpha=0.5,label='detection threshold ({:.1F} ppm)'.format(threshold/1e-6))
    plt.legend(handles=[legend_patch],frameon=False)
