#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  2 04:42:49 2021

@author: rajananderson
"""
from accretiondisk import AccetionDisk 
import matplotlib.pyplot as plt
from RPAStructure.main import main as RPA
import os


os.chdir('RPAStructure')
print(RPA(0.1,0.2,'SLY4',2,'NONE',True))
os.chdir('..')

print(os.getcwd())

M3a0_datasetdict={
    
    'axis_of_neutrino_emission': 1,
    'temp_idx':2,
    'density_idx':3,
    'electronfrac_idx':4,
    'addtional_spatial_dims':[0]

    }


timeY0_datasetdict={
    
    'axis_of_neutrino_emission': 2,
    'temp_idx':3,
    'density_idx':5,
    'electronfrac_idx':4,
    'addtional_spatial_dims':[0,1]

    }
        
NL3_datasetdict={
    
    'axis_of_neutrino_emission': 2,
    'temp_idx':3,
    'density_idx':5,
    'electronfrac_idx':4,
    'addtional_spatial_dims':[0,1]

    }
        
param={
       'preprocessed_crosssections':True,
       
       
       
       'crossparam':
           {
               'method':'viri',
               'viriparam':
                   {
                       'SVon':True,
                       'SAfit':True,
                       'SAlow':False,
                       
                       },
               
           }
       }

eflavs=['Mu and Tau','Electron','Anti Electron']
filenames=['time6500Y0','time3500Y0','time4500Y0','time5500Y0','M3a0.8t20','M3a0.8t60','NL3_t2.5ms_Y0']
filenames=['time3500Y0','time4500Y0','time5500Y0']
eflavs=['Electron']
'''
for filename in filenames:
    for eflav in eflavs:
        fileproLow='data/{}/SaLowSV/crossoptdata.dat'.format(filename)
        fileproAxFit='data/{}/AxFitSV/crossoptdata.dat'.format(filename)
        
        acdLow=AccetionDisk(param,fileproLow,timeY0_datasetdict)
        acdFit=AccetionDisk(param,fileproAxFit,timeY0_datasetdict)
        limits=(0.01,3)
        
        SaLow,depths=acdLow.averageenergycompared(eflavour=eflav,surfdepthlimits=limits)
        AxFit,depths=acdFit.averageenergycompared(eflavour=eflav,surfdepthlimits=limits)
        plt.clf()
        
        plt.figure(figsize=(13,7),dpi=200)

        plt.title(filename+' '+ eflav + ' Neutrinos')
        plt.plot(depths,SaLow,label='Pure Virial Expansion')
        plt.plot(depths,AxFit, label='Virial with Fit')
        plt.xlabel('Optical Depth')
        plt.ylabel('Interacting/Non-Interacting Average Energy')
        plt.legend()
        plt.show()
        plt.clf()
        
'''
'''
for filename in filenames:
    fileproLow='data/{}/SaLowSV/crossoptdata.dat'.format(filename)
    fileproAxFit='data/{}/AxFitSV/crossoptdata.dat'.format(filename)
    
    acd=AccetionDisk(param,fileproAxFit,timeY0_datasetdict)
    plt.clf()
    acd.surface(itypes=[12,15],surfacedepth=2/3)
    print(acd.luminocity(itype=[10]),acd.totflux(itype=[10]),'non',acd.averageenergy(itype=[10]))
    print(acd.luminocity(itype=[13]),acd.totflux(itype=[13]),'on',acd.averageenergy(itype=[13]))
'''