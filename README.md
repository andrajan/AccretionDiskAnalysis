# AccretionDiskAnalysis

In this package various tools are provided to vissualize information about neutrinos in neutron star merger events and evaluate the effect of different models of nucleonic interaction on that information.

## Set Up

1. Use environment files to download appropriate conda environment for your system. Make sure you have gfortran installed. requirements.txt is also provided for the python packages.
2. First we compile edense.f90. To do this enter the following 

          f2py -m edense -c edense.f90
 
Proper compilation can be can be checked by running nointercross.py
3. Now navigate to the RPAStructure and type make. setup.py can be used to very the correct behaviour is observed.

## Index

data: This folder contains all data we wish to process. Opening this folder will give us a list folders named according to the hydrodynamical simulation that created them and the moment/slice of that simulation we are looking at. If we enter one of these folders, we see one file; the original data from the hydrodynamic simulation and then more folders labelled according to the model of nuclear reaction used to post-process them. In all of these folders there should be just one file called crossoptdata.dat
    
        -data
            -HydrodynamicSimulation
                -Nuclear Force
                    -postprocesseddata
    
acrretiondisk.py: Contains class AccetionDisk which manipulates data from the data folder to calculate different values assosiated with the neutrinos of the system for a specific nuclear force as well as create plots of the neutrino surface.
    
AvEnergyCompare.py: This file contains functions that will create visualizations to compare different methods of nucleonic interaction and see what their effect on different neutrino properties are.
    
structuremain.py: this file allows us to vissualize how different properties such as structure factors change according to location in the disk.
    
structfact.py: this file calculates structure factors based on Horowitz et al (2017)
    
datadictionaries.py: this contains a bunch of dictionaries that contain info about how data in a specific hydrodynamic simulation is ordered so that our AccetionDisk class can deal with it.
    
edense.f90: this file calculates electron and neutrino densities so that we can determine the mean free path from the cross section for annihilation and neutrino electron scattering events. Provides a test whether f2py is working on your system other than the more complex makefile in RPAStructure.
    
nointercross.py: this contains all the cross sections asides from those calculated using RPA.
    
RPAStructure : this folder contains code we have borrowed from Pastore that we have adapted for use with our systems. 
    
   -main.f90: this is the main interface we call from python
        
   -asym.f90: this file calculates our cross section.
        
   -asym_cstes.f90: this file sets a bunch of constants
        
   -asym_gauss.f90: this file contains a function and tools to implements gaussian integration
        
   -asym_inm.f90: this file provides subroutines to calculate many properties of nuclear matter(eg. effective mass)
        
   -asym_io.f90: input output functionality: Not used as currently configured as the input/output is done with python.
        
   -asym_param.f90: sets up parameters for use in our matrices based on input (readforce,type (eg. TENS, SORB ....)).
        
   -beta.f90: calculates various energy integrals that will make up the system of linear equations to be solved
        
   -wi.f90: sets up the matrices of the linear problem to be solved.
        
   -determinante.f90: solves linear equation using cramer's method based on matrices set up by wi.f90
        
   -----------------------------------------------------------------------------

   -forces.param/functional.param : contains parameter weights
        
   -pythonRPAtest.py: some code to examine static structure functions from RPA.
        
   -setup.py: This provides some functions that test that all is working well with our RPA code (energy squared dependence check, as well as check that non-interating and bare casses converge in certain limits), and that it is properly connected to python.
