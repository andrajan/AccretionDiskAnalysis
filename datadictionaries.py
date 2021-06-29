#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 12:14:17 2021

@author: rajananderson
"""
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

# the following dictionary can be used after time consuming process of calculating all cross sections and optical depths has been acomplished 
preparam={
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
