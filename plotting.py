# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 22:33:24 2018

@author: pnter
"""
import matplotlib.pyplot as plt
from math import ceil

def plotResults(axes,log,stats,title=None,legend=False):
    for i in range(4):
        axes.plot(range(len(log)),stats[:,i])
    if legend:
        axes.legend(['Susceptible','Exposed','Infected','Removed'],loc='upper right')
    plt.xlabel('t')
    plt.ylabel('Population')
    if type(title) == str:
        plt.title(title)

def plotExperiment(experiment,titles):
    f = plt.figure(figsize=(12,3))
    axesList = []
    numSims = len(experiment.stats)
    for i in range(numSims):
        axesList.append(f.add_subplot('{0}{1}{2}'.format(ceil(numSims/4),
                                      4,i)))
        plotResults(axesList[i],experiment.logs[i],experiment.stats[i],
                    titles[i])