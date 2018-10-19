# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 22:32:32 2018

@author: pnter
"""

# -*- coding: utf-8 -*-
"""
@author: Nick Terry

A test script for seir-sim that generates a random contact network and runs the
simulation
"""
from experiment import Experiment
from seirSim import SeirSim
from policies import vaccinateTopNTSH,vaccinateNAcquaintance
from tallyFuncs import numNodesInState
from networkalgs import generateRandomGraph
import numpy as np
import matplotlib.pyplot as pl
from scipy.interpolate import UnivariateSpline

numNodes = 1000

#Use parameters estimated from Ebola epidemics in Sierra Leone
#Source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4169395/
exposureRate = .45 #Beta
infectionRate = 1/5.3 #Sigma
recoveryRate = 1/5.61 #Gamma

numIterations = 100



def comparePolicies(k):
    
    numVaccinated = k
    
    TSHPolicy = vaccinateTopNTSH(numVaccinated)
    AcqPolicy = vaccinateNAcquaintance(numVaccinated)
    policies = [[TSHPolicy],[AcqPolicy]]
    
    #Store statistics on proportion of population infected (Average,Var)
    proportionInfected = np.zeros([numIterations,len(policies),2])
    #Store count of # nodes infected in each simulation
    numInfected = np.zeros([numIterations,len(policies)])
    
    for iteration in range(numIterations):
        print("Beginning Simulation {0}".format(iteration))
        G = generateRandomGraph(numNodes)
        
        simulation = SeirSim(G,exposureRate,infectionRate,recoveryRate,
                             logSim=True, 
                      tallyFuncs=[numNodesInState(0),numNodesInState(1),
                                  numNodesInState(2),numNodesInState(3)])
        
        exper = Experiment(simulation)
        
        logList,statsList = exper.compare(policies)

        for policy in range(len(policies)):
            propInfectedMean = np.mean(statsList[policy][:,2]/numNodes)
            propInfectedVar = np.var(statsList[policy][:,2]/numNodes)
            proportionInfected[iteration,policy,:] = np.array([propInfectedMean,propInfectedVar])
            
            numRemoved = np.sum(np.where(logList[policy][-1][0]==3,
                                       np.ones(numNodes),np.zeros(numNodes)))
            numInfected[iteration,policy] = numRemoved - numVaccinated
    
    #Compute mean, variance, 90% CI for each policy
    numInfectedStatsL = np.zeros([2,2])
    propInfectedStatsL = np.zeros([2,2])
    for i in range(len(policies)):
        mean = np.mean(proportionInfected[:,i,0])
        var= np.var(proportionInfected[:,i,0])
        mean2 = np.mean(numInfected[:,i])
        var2 = np.var(numInfected[:,i])
        
        propInfectedStatsL[i,:] = np.array([mean,var])
        numInfectedStatsL[i,:] = np.array([mean2,var2])
        
    return propInfectedStatsL,numInfectedStatsL

propInfectedStats = np.zeros([numNodes,2,2])
numInfectedStats = np.zeros([numNodes,2,2])

for k in range(numNodes):
    print('Beginning Simulation for K={0}...'.format(k))
    stats = comparePolicies(k)    
    propInfectedStats[k,:,:] = stats[0]
    numInfectedStats[k,:,:] = stats[1]

titles = ['Two-Step Heuristic Immunization','Acquaintance Immunization']
ax = pl.subplot(1,1,1)
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])
    
ax.plot(range(numNodes),propInfectedStats[:,0,0],' . ',
        linestyle='None')
ax.plot(range(numNodes),propInfectedStats[:,1,0],' . ',
        linestyle='None')
# Put a legend below current axis
ax.legend(titles,loc='upper center', bbox_to_anchor=(0.5, -0.15),
          fancybox=True, shadow=True, ncol=2)

pl.ylabel('Mean Proportion of Population Infected')
pl.xlabel('K')

pl.title('Comparison of TSH and Acquaintance Immunization for Varying K:\n Proportion of Population Infected')

pl.figure()
ax2 = pl.subplot(1,1,1)

box = ax2.get_position()
ax2.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])
    
ax2.plot(range(numNodes),numInfectedStats[:,0,0],' . ',
        linestyle='None')
ax2.plot(range(numNodes),numInfectedStats[:,1,0],' . ',
        linestyle='None')

ax2.legend(titles,loc='upper center', bbox_to_anchor=(0.5, -0.15),
          fancybox=True, shadow=True, ncol=2)

pl.ylabel('Mean Number of Infected Nodes')
pl.xlabel('K')

pl.title('Comparison of TSH and Acquaintance Immunization for Varying K:\nTotal Nodes Infected')

pl.figure()
diffNumInfected = numInfectedStats[:,1,0]-numInfectedStats[:,0,0]
pl.plot(range(numNodes),diffNumInfected,' . ',
        linestyle='None')
pl.ylabel('Difference in Total Nodes Infected')
pl.xlabel('K')
pl.title('Difference in Total Nodes Infected Using TSH vs. Acquaintance Immunization')
pl.show()

pl.figure()
ax2 = pl.subplot(1,1,1)

box = ax2.get_position()
ax2.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])
cvNumInfected = np.sqrt(numInfectedStats[:,:,1])
ax2.plot(range(numNodes),cvNumInfected[:,0],' . ',
        linestyle='None')
ax2.plot(range(numNodes),cvNumInfected[:,1],' . ',
        linestyle='None')

ax2.legend(titles,loc='upper center', bbox_to_anchor=(0.5, -0.15),
          fancybox=True, shadow=True, ncol=2)

pl.ylabel('Mean Number of Infected Nodes')
pl.xlabel('K')

pl.title('Comparison of TSH and Acquaintance Immunization for Varying K:\nStandard Deviation for Total Nodes Infected')