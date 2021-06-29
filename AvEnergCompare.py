#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  2 04:42:49 2021

@author: rajananderson

In this module we provide functions which may be useful to create functions that visualize differences
in models of interaction

"""
from accretiondisk import AccetionDisk 
import matplotlib.pyplot as plt
from RPAStructure.main import main as RPA
import os
import sys
import numpy as np
import datadictionaries as dicts # options are M3a0_datasetdict , timeY0_datasetdict , or 'NL3_datasetdict'


print(os.getcwd())
os.chdir('RPAStructure')
os.chdir('..')

print(os.getcwd())
        

def avgenergycompare(filenames,types,datasetdict,eflavs=['Mu and Tau'],limits=(0.01,2000),saveaspng=None):
    '''
    
    This function plots ratio of average energy for interacting nucleons of given mehtonds over the non-interacting case
    fir a given interval of optical depth.
    

    Parameters
    ----------
    filenames : Names of all hydrodynamical simulations investigating as organized in file system .
    types : Type of interaction as named in file ordering system.
    datasetdict : dictionary corresponding to particular hydrodynamic simulation.
    eflavs : TYPE, optional : The default is ['Mu and Tau']
        DESCRIPTION. For what flavour of neutrinos do we want to generate graphs for (a list of flavours will
            generate many different graphs). 
    limits : TYPE, optional : The default is (0.01,2000).
        DESCRIPTION. What should the limits of the optical depth be
    saveaspng : TYPE, optional : Do not save
        DESCRIPTION. If string is given then the function will save the figure as a pngs with that as the first part of the string

    Returns
    -------
    PlottedFigures.

    '''
    for filename in filenames:
        for eflav in eflavs:
    
            #more interactions
            filepro=[]
            acd=[]
            
            plt.clf()
            plt.figure(figsize=(13,7),dpi=700)
            plt.title('Average Energy at given Optical depth for '+ eflav + ' Neutrinos')

            for i,typ in enumerate(types):
                filepro.append('data/{}/{}/crossoptdata.dat'.format(filename,typ))   
                acd.append(AccetionDisk(dicts.preparam,filepro[i],datasetdict))
                crosssec,depth=acd[i].averageenergycompared(eflavour='Mu and Tau',surfdepthlimits=limits)
                plt.plot(depth,crosssec,label=typ)

            
            plt.xlabel('Optical Depth')
            plt.ylabel('Interacting/Non-Interacting Average Energy')
            plt.legend()
            if (saveaspng!=None):
                plt.savefig(saveaspng+filename+eflav+'.png')
            plt.show()
    
            plt.clf()
        
    
#avgenergycompare(filenames,types,dicts.timeY0_datasetdict,saveaspng=None)
def surfacecompared(filenames,types,datasetdict,eflav='Mu and Tau',elastic=True,saveaspng=None):
    '''
    This function plots neutrino surfaces and compares them

    Parameters
    ----------
    filenames : Names of all hydrodynamical simulations investigating as organized in file system .
    types : Type of interaction as named in file ordering system.
    datasetdict : dictionary corresponding to particular hydrodynamic simulation.
    eflav : TYPE, optional : The default is 'Mu and Tau'.
        DESCRIPTION. Choose betwen 'Mu and Tau','Electron' and 'Anti Electron'
    elastic : TYPE, optional : plot the non-interacting elastic neutrino surface
        DESCRIPTION. The default is True.
    saveaspng : TYPE, optional : Do not save
        DESCRIPTION. If string is given then the function will save the figure as a pngs with that as the first part of the string

    Returns
    -------
    Figure.
    
    '''   

 
    for filename in filenames:
        filepro=[]
        acd=[]
        
        if eflav=='Mu and Tau':
           itypes=[15]
        elif eflav=='Electron':
           itypes=[13]
        elif eflav=='Anti Electron':
           itypes=[14]
        plt.figure(figsize=(13,7),dpi=700)
        
        for i,typ in enumerate(types):
            filepro.append('data/{}/{}/crossoptdata.dat'.format(filename,typ))   
            acd.append(AccetionDisk(dicts.preparam,filepro[i],datasetdict))
            acd[i].label=typ
            acd[i].surface(itypes=itypes,surfacedepth=2/3)    
        if elastic:
            acd[0].surface(itypes=[12],surfacedepth=2/3) 
        plt.title(filename+' '+ eflav+'Neutrino Surfaces')

        plt.show()
        if (saveaspng!=None):
            plt.savefig(saveaspng+filename+eflav+'.png')
        
        plt.clf()

def probeopticaloncrosssurface(filename,types,datasetdict,eflav='Mu and Tau',elastic=False):
    '''
    This function plots a surface plot of a ratio of cross sectionss using different models of interaction.
    Then neutrino surfaces at varying optifcal depths are plotted along that surface plot. Neutrino surfaces
    are plotted based on the first item in the list of types.

    Parameters
    ----------
    filename : folder of hydrodynamical simulation of interest.
    
    types : input a list of the two models of interest. First item in the list will be in the denominator
        
    datasetdict : select appropriate dictionary for the hydrodynamic simulation
    
    eflav : TYPE, optional :The default is 'Mu and Tau'.
        DESCRIPTION. Choose betwen 'Mu and Tau','Electron' and 'Anti Electron'
    
    elastic : TYPE, optional : The default is False
        DESCRIPTION. If turned on the first interaction type in the list will be ignored and this function
        will create a surface plot using the non-interacting elastic cross section and the second item in 
        your list.

    Returns
    -------
    None.

    '''
    datatend='/crossoptdata.dat'
    filepros=['data/'+filename+'/'+t+datatend for t in types]
    if eflav=='Mu and Tau':
        m=12
    elif eflav=='Electron':
        m=11
    elif eflav=='Anti Electron':
        m=10
    plt.figure(figsize=(13,7),dpi=700)
    if elastic:
        n=m-3
    else:
        n=m
    acd=AccetionDisk(dicts.preparam,filepros[0],datasetdict)
    acd2=AccetionDisk(dicts.preparam,filepros[1],datasetdict)
    x=acd.fulldata[:,0]
    y=acd.fulldata[:,-1]
    z=acd.fulldata[:,n-3]/acd2.fulldata[:,m-3]
    z[z<0]=0
    z=z**-1
    nz=z
    #z[(z<1.03) & (z>0.97)]=1
    
    X=np.reshape(x,(len(list(set(list(x)))),len(list(set(list(y))))))
    Y=np.reshape(y,(len(list(set(list(x)))),len(list(set(list(y))))))
    Z=np.reshape(z,(len(list(set(list(x)))),len(list(set(list(y))))))
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #ax.plot_surface(X,Y,Z,cmap=cm.coolwarm)

    depths=[0.1,2/3,1.5,100]
    surflab=['1/10','2/3']+depths[2:]
    idx=[]
    for i,surf in enumerate(depths):
        idx=acd.surface(itypes=[m],surfacedepth=surf,plotyes=False)
        xx=x[idx]
        yy=y[idx]
        h=z[idx]
        ax.plot(xx,yy,zs=h,color=((0.7**i, 1-0.7**i, 1-0.7**i, 1)),label='optdepth={}'.format(surflab[i]))
    ax.plot_surface(X,Y,Z,color=(255/255, 255/255, 0/255, 0.3))


    ax.view_init(elev=20,azim=30)
    ax.set_zlabel('{} / {} Crosssections'.format(types[1],types[0]))
    ax.set_ylabel('z(km)')
    ax.set_xlabel('y(km)')
    ax.set_title('How far do we need to probe to see a difference ({})'.format(eflav))
    ax.legend(fontsize='small',loc='upper right')
    plt.show()

if __name__=='__main__': 
    eflavs=['Mu and Tau','Electron','Anti Electron']
    filenames=['time3500Y0','time3500Y0','time4500Y0','time5500Y0','M3a0.8t20','M3a0.8t60','NL3_t2.5ms_Y0']
    #filenames=['NL3_t2.5ms_Y0']#,'time4500Y0','time5500Y0']
    eflavs=['Mu and Tau']
    types=['RPABARE','RPASLY4','RPABARE','RPASLY4','SLY4T']#,'f0']
    
    probeopticaloncrosssurface(filenames[0],types,dicts.timeY0_datasetdict)
