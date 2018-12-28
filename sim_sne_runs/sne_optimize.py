#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 2018

sne_optimize.py

@author: jhan
"""

import numpy as np
import copy
from sys import argv
from os import system as sys
from os import chdir
from os import environ
from scipy.stats import norm

#%% Pseudo-Code
"""
1. Params:
    1. True x1  ~ N(0,1)
    2. True c ~ N(0, 0.1^2)
    3. True E(B-V) = exp(0.1)
    4. True mag = -19.1 - 0.14*x1_true + 2*c_true, 4.1*E(B-V)_true
2. Observe
    1. 
3. Truncate an observed value (probabilistically)
    1. Use UNITY (error function for truncation)
4. Write light-curve fit results
5. Run UNITY
6. Repeat w/ more SNe

"""
#%% Define Directories
root = environ['ROOTDIR']
#print(root)
dir_param_root = root+"param_folder/"
dir_sne = root+"sim_sne/"
dir_home = root+"sim_sne_runs/"
dir_union = root+"union3/"

#%% Define Functions
def genParams(Nsne):
    x1mean, x1std = 0, 1 # Normal distribution
    cmean, cstd = 0, 0.1 # Normal distribution
    ex_scale = 0.1 # Exponential distribution
    dist_modulus = 40
    # Draw parameter samples from the population above. 
    t_x1s = [x1std*np.random.randn()+x1mean for i in range(Nsne)]
    t_cs = [cstd*np.random.randn()+cmean for i in range(Nsne)]
    t_exs = [(ex_scale)*np.random.exponential(1) for i in range(Nsne)]
    t_mags = [-19.1-0.14*t_x1s[i]+2*t_cs[i]+4.1*t_exs[i]+dist_modulus for i in range(Nsne)]
    ## Nomenclature: plist <-> param_list (to save typing)
    t_plist = [[t_x1s[i],t_cs[i],t_exs[i],t_mags[i]] for i in range(Nsne)]
    return t_plist

def obsErve(t_plist, siglist):
    ## Input order" [x1, c, ex, mag]
    Nsne = len(t_plist)
    xsig, csig, esig, msig = siglist
    ## Introduce observational uncertainties
    xs, cs, es, ms = xsig, csig, esig, msig
    o_x1s = [float(t_plist[i][0])+xs*np.random.randn() for i in range(Nsne)]
    o_cs = [float(t_plist[i][1])+cs*np.random.randn() for i in range(Nsne)]
    o_exs = [float(t_plist[i][2])+es*np.random.randn() for i in range(Nsne)]
    o_mags = [float(t_plist[i][3])+ms*np.random.randn() for i in range(Nsne)]
    o_plist = [[o_x1s[i],o_cs[i],o_exs[i],o_mags[i]] for i in range(Nsne)]
    return o_plist

def normProb(val, cut, scl):
    # Inputs: current value, threshold, decay scale
    prob = 100 # Not normalized, check value
    prob = norm.cdf(val, loc=cut, scale=scl)
    assert ((prob>=0)and(prob<=1)), "Probability is not normalized"
    return prob

def truncProb(prob): # Pass 1 (leave) or 0 (truncate) with probability
    rand = np.random.uniform(low=0.,high=1.)
    val = 100 # Not boolean, check value
    if rand>prob:
        val = 1
    else:
        val = 0
    assert ((val==0)or(val==1)), "Issue with truncProb function"
    return val

def truncAte(o_plist, cut, scl): # Return truncated list
    # List comprehensions are badass. 22, 0.5
    maglist = [o_plist[i][3] for i in range(len(o_plist))] # Observed mag list
    problist = [normProb(maglist[i], cut, scl) for i in range(len(o_plist))]
    deathnote = [bool(truncProb(problist[i])) for i in range(len(problist))]
    trunc_list = []
    for i in range(len(problist)):
        if deathnote[i]:
            trunc_list.append(o_plist[i])
    return trunc_list

def dx0(mag, dm):
    return 0.4*np.log(10)*10**(-0.4*mag)*dm

def paramToObs(plist, siglist): # Take in the truncated parameters and output observables
    ## Input: [x1, c, ex, mag], [x1_unc, c_unc, ex_unc, mag_unc]
    ## Output: [X0,X0_unc,X1,X1_unc,Color,Color_unc, mag, mag_unc]
    # Relations: x0_obs = 10^{-0.4*mag}, color = c + ex, dx0 = 0.4*ln(10)*10^-0.4m*dm
    obs = [[10**(-0.4*plist[i][3]),dx0(plist[i][3],siglist[3]), 
           plist[i][0],siglist[0],plist[i][1]+plist[i][2],
           0.04, plist[i][3], siglist[3]] for i in range(len(plist))]
    return obs

def mkFol(ID):    # Make folder with SNe ID
    # Inputs: ID. Pretty straightforward.
    ID = int(ID) # Just to make sure.
    sys("mkdir %s%d" %(dir_sne,ID))
    print("Folder for SNe %d Created." %ID)
    return 0

def writeFit(ID, observables): # Writes three parameter files to ID directory
    ## Inputs:
    # ID: Integer. Duh. Starts from 0.
    # Observables: [X0,X0_unc,X1,X1_unc,Color,Color_unc, m, m_unc] ... in that order!
    snedir = "%d/" %int(ID)     # In case you forgot to pass an integer... 
    sys("cp %s* %s" %(dir_param_root, dir_sne+snedir)) # Copy in the param files

    ## Re-write the result_salt2.dat file based on the simulated observables
    x0, x0_unc, x1, x1_unc, color, color_unc, mag, mag_unc = observables
    x0_cov, x1_cov, color_cov = x0_unc**2,x1_unc**2,color_unc**2
    # Value keys & Replacements
    val_keys = ["REPLACE_X0","REPLACE_UNC_X0","REPLACE_X1",
                "REPLACE_UNC_X1","REPLACE_COLOR","REPLACE_UNC_COLOR",
                "REPLACE_MAG", "REPLACE_UNC_MAG"]
    val_replace = [x0, x0_unc, x1, x1_unc, color, color_unc, mag, mag_unc]
    # Covariance keys & Replacements
    cov_keys = ["REPLACE_COV_X0","REPLACE_COV_X1","REPLACE_COV_COLOR"]
    cov_replace = [x0_cov, x1_cov, color_cov]
    # Salt2, for some reason, prints results twice...
    # But hey, at least the replacements are the same.
    # Combine everything into one dictionary
    total_keys = val_keys + cov_keys
    total_replace = val_replace + cov_replace
    find = dict(zip(total_keys, total_replace))
    # Write an updated resulst_salt2.dat file into sne directory
    with open(dir_param_root+'result_salt2.dat', "r") as oldf:
        with open(dir_sne+snedir+'result_salt2.dat', 'w') as newf:
            for line in oldf:
                check=0
                newl = copy.copy(line)
                for key in total_keys:
                    if key in line:
                        newl = newl.replace(key, str(find[key]))
                        check=1
                if check==1:
                    newf.write(newl)
                if check==0:
                    newf.write(line) # When nothing needs to be replaced.
    ## Assertion tests! Don't crash the computer when u cross the date line.
    with open(dir_param_root+'result_salt2.dat', "r") as f:
        og_numline = sum(1 for line in f)
    with open(dir_sne+snedir+'result_salt2.dat', 'r') as f:
        new_numline = sum(1 for line in f)
    assert (og_numline == new_numline), "New result_salt2.dat does not match original file length."
    junk_left_check = 0    
    with open(dir_sne+snedir+'result_salt2.dat', 'r') as f:
        for line in f:
            for key in total_keys:
                if key in line:
                    junk_left_check+=1
    assert (junk_left_check == 0), "Incomplete replacement of result_salt2.dat: %d remaining lines" %junk_left_check
    
    ## Now, re-write the lightfile
    key, value = "REPLACE_MASS", str(np.random.randn() + 10)
    with open(dir_param_root+'lightfile', "r") as oldf:
        with open(dir_sne+snedir+'lightfile', "w") as newf:
            for line in oldf:
                newl = copy.copy(line)
                if key in line:
                    newl = newl.replace(key, value)
                newf.write(newl)
    return 0

def makeParamFile(Nsne):
    # Write the correct paramfile.txt & sim_sne_list.txt to run UNITY
    # Rember that it's only ONE paramfile.txt for an entire sample of NSNE.
    sys("ls -d %s* > %ssim_sne_list.txt" %(dir_sne, dir_home))
    key = "REPLACE_LIST"
    snelist = "%ssim_sne_list.txt" %dir_home
    with open(dir_param_root+"paramfile.txt",'r') as f:
        with open(dir_home+"sim_sne_paramfile.txt",'w') as g:
            for line in f:
                newl = copy.copy(line)
                if key in line:
                    newl = newl.replace(key, snelist)
                    g.write(newl)
                else:
                    g.write(line)
    with open(dir_param_root+"paramfile.txt",'r') as f:
        og_numline = sum(1 for line in f)
    with open(dir_home+"sim_sne_paramfile.txt",'r') as f:
        new_numline = sum(1 for line in f)
    assert (og_numline == new_numline), "Paramfile length does not match original."
        
    return 0

def modMagCut(mag):
    key = "REPLACE_MAG"
    with open(dir_param_root+"mag_cuts.txt_mod", "r") as f:
        with open(dir_param_root+"mag_cuts.txt", "w") as g:
            for line in f:
                if key in line:
                    line = line.replace(key, str(mag))
                    g.write(line)
                else:
                    g.write(line)
    return 0

def NtoTexp(N):
    const = 300  # NSNe * texp = constant
    texp = const/N # Baseline: when N=300, uncertainties are at default values
    return texp

def NtoSig(N): # Relationship between NSNe and the uncertainties
    texp = NtoTexp(N)
    defaults = [0.5, 0.028, 0.028, 0.2]
    siglist = [defaults[i]/np.sqrt(texp) for i in range(len(defaults))]
    return siglist

def NtoLim(N): # Relationship between NSNe and the limiting magnitude
    texp = NtoTexp(N)
    c = 22  # Baseline: when N=100, t=1, and lim = 22
    lim = c + 1.25*np.log10(texp)
    return lim

def runUnity(Nsne): # Run UNITY given a number of SNe to simulate.
    Nsne = int(Nsne) # Just to make sure.
    rundir = "%s%d/" %(dir_home, Nsne)
    sys("mkdir %s" %rundir)
    chdir(rundir) # Run in the Nsne subdirectory
    
    siglist = NtoSig(Nsne) # Set observational uncertainties
    cut, scl = NtoLim(Nsne), 0.5 # Set magnitude cutoff
    
    t_plist = genParams(Nsne) # Generate true paramters
    o_plist = obsErve(t_plist, siglist) # Introduce observational uncertainties
    trunced = truncAte(o_plist, cut, scl) # Truncation effects
    obslist = paramToObs(trunced, siglist) # Generate final observed parameters
    
    trunc_SNe = len(obslist) # Truncated sample has less then original Nsne
    print("Number of runs after truncation: %d" %trunc_SNe)
    sys("rm -r %s*" %(dir_sne)) # Clear out previous run
    for i in range(trunc_SNe): # Prepare right files and folders for UNITY.
        mkFol(i) # Make folders
        writeFit(i, obslist[i]) # Write the necessary parameter files
    makeParamFile(trunc_SNe) # Make parameterfile
    modMagCut(cut) # Modify mag_cut for current NSNe
    
    scriptpath = dir_union + "scripts/read_and_sample.py"
    parampath = dir_home+"sim_sne_paramfile.txt"
    logpath = rundir + "log_%d.txt" %Nsne
    sys("python2 %s %s 2 > %s" %(scriptpath, parampath, logpath))

    return 1
#%%

assert (len(argv)==2), "Provide the following input(s): NSNe."
NSNE = int(argv[1])
checkval = 1
checkval = runUnity(NSNE)
assert (checkval==1), "Execution of runUNITY has failed."















