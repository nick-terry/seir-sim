# -*- coding: utf-8 -*-
"""
Created on Fri May 25 11:51:16 2018

@author: pnter
"""

from pyvis.network import Network

def showNetwork(G):
    
    net = Network("500px","500px")
    net.from_nx(G)
    net.show("graph.html")