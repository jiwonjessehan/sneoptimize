#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 10 17:01:04 2018

@author: jhan
"""
import numpy as np
import cPickle as pickle
import matplotlib.pyplot as plt
import gzip

#%% Define Directories
dir_param_root = "/Users/jhan/Desktop/UNITY/param_folder/"
dir_sne = "/Users/jhan/Desktop/UNITY/sim_sne/"
dir_home = "/Users/jhan/Desktop/UNITY/sim_sne_runs/"
dir_union = "/Users/jhan/Desktop/UNITY/union3/"

#%%


def getMBunc(NSNE, fit_params):    
    mbunc = np.std(fit_params['MB'])
    assert (mbunc>0), "Negative uncertainty"
    return mbunc

def getMuZunc(NSNE, fit_params):
    mbunc = np.std(fit_params['mu_zbins'])
    assert (mbunc>0), "Negative uncertainty"
    return mbunc

def getBetaB(NSNE, fit_params):
    alph = np.std(fit_params['beta_B'])
    assert (alph>0), "Negative uncertainty"
    return alph

def getNumber(NSNE, fit_params):
    size = len(fit_params['true_x1'][0])
    return size

def doStuff(Nlist):
    actNumlist = []
    muzlist = []
    betalist = []
    for num in Nlist:
        sample_pckl = dir_home+"%d/samples_sim_sne_list.pickle" %(num)
        fit_params = pickle.load(gzip.open(sample_pckl, 'rb'))
        actNum = getNumber(num, fit_params)
        actNumlist.append(actNum)
        
        muz = getMuZunc(num, fit_params)
        muzlist.append(muz)

        beta = getBetaB(num, fit_params)
        betalist.append(beta)
        
    fig = plt.figure(1, figsize=(10,10))
    ax = fig.add_subplot(111)
    ax.plot(Nlist, actNumlist, "x--")
    fig.savefig("N_actN.png")

    fig2 = plt.figure(2, figsize=(10,10))
    ax2 = fig2.add_subplot(111)
    ax2.plot(actNumlist, muzlist, "x--")
    fig2.savefig("actN_muZ.png")    

    fig3 = plt.figure(3, figsize=(10,10))
    ax3 = fig3.add_subplot(111)
    ax3.plot(actNumlist, betalist, "x--")
    fig3.savefig("actN_betaB.png")
    
    return 0


Nlist = [300, 500, 800, 900, 1000, 1200, 1400]

doStuff(Nlist)





