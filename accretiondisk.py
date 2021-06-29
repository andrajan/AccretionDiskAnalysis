#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 11:13:25 2021

@author: rajananderson


This module contains the class AccretionDisk which can be used to preprocess data, as well as calculate
and plot many aspects of the accretion disk. A template parameter file is given at the bottom of this
module which can be edited to control the functionality of the object.

"""
import numpy as np
import nointercross as cr
import edense as ed
import scipy.integrate as integ
import scipy.special as sp
import structfact as st
from scipy.interpolate import Rbf
import os
import matplotlib.pyplot as plt
from RPAStructure.main import main as RPA
import time
from matplotlib import cm
import datadictionaries as dicts
import sys


print(os.getcwd())

class AccetionDisk():
    
   '''
   Class that takes hydrodynamical simulation data and turns processes it.  
   '''
   
   def __init__(self,param,datapath,datasetdict,nameofsave=None):
       '''
       

       Parameters
       ----------
       param : Tells class whether the data is processesd and how it should be processesd.
       
       datapath : where to find our data.
       
       datasetdict : The dictionary specific to this hydrodynamic simulation.
       
       nameofsave : What our processed data should be called if not already processed. Default is none,
           which results in the terminal prompting you for a savename after completion


       '''
       
       
       super(AccetionDisk, self).__init__()
       
       #Caution this dictionary may give wrong labels if using additional data
       self.dictforplot={0 : 'Non Interacting', 1: 'Interacting',
                         2 : "Electron Anti-Neutrino",3 : 'Electron Neutrino',
                         4 : 'Mu and Tau Neutrino'}
       
       self.label='' # to be used in order to better control plot labels
       
       self.spatialdim=len(datasetdict['addtional_spatial_dims'])+1
       self.pathohome="/".join(datapath.split('/')[:-1])+'/'
       self.interpolationcoefficient=10
       
       if param['preprocessed_crosssections']:
           self.fulldata=np.loadtxt(datapath)
           print(self.fulldata.shape)
       else:
           self.fulldata=generate_data(datapath,param['crossparam'],datasetdict)
           print(self.fulldata.shape)   
           #save_data for later use
           if nameofsave==None:
               print('enter name of dataset')
               directory=self.pathohome+str(input())   
               if not os.path.exists(directory):
                   os.mkdir(directory)
           else:
               directory=self.pathohome+nameofsave
               if not os.path.exists(directory):
                   os.mkdir(directory)

             
           filetosave=directory+'/crossoptdata.dat'
           print(os.getcwd())

           file=open(filetosave,'w')
           np.savetxt(file,self.fulldata)
           file.close()

            

           
           
   def surface(self,itypes=[12,15],interpolate=False,surfacedepth=2/3,plotyes=True):
      '''
       

       Parameters
       ----------
       itypes : TYPE, optional : [12,15]
           DESCRIPTION. Chooses what type of neutrinos we should plot the surface of
           
       interpolate : TYPE, optional : Default is False.
           DESCRIPTION. Caution Does not workwith all additional functions. Tells surface function to 
               interpolate in order to get a higher definition surface
           
       surfacedepth : TYPE, optional: The default is 2/3.
           DESCRIPTION. The optical depth at which we define the neutrino surface
           
       plotyes : TYPE, optional :  The default is True. A neutrino surface plot is created
           DESCRIPTION. The default is True.

       Returns
       -------
       surfindx : 
           DESCRIPTION. An index of where the in our self.fulldata the surface is located at

      '''
      if itypes==[10,13] or itypes==[13]:
          print(itypes,'whyt')
          eflavour='Electron'          
      elif itypes==[11,14] or itypes==[14] :
          eflavour='Anti Electron'
      elif itypes==[12,15] or itypes==[15]:
          eflavour='Mu and Tau'
      else:
          eflavour=str(itypes)
      if plotyes:
          print(eflavour,itypes)
          plt.title(eflavour +' Surfaces')
          plt.xlabel('radial distance from center of system (km)')
          plt.ylabel('Z distance from plane of rotation(km)')
          
      for j,i in enumerate(itypes):
           cross=self.fulldata[:,i]
           x=self.fulldata[:,-1]
           y=self.fulldata[:,0]

           if interpolate:               
               x=self.fulldata[:,-1]
               y=self.fulldata[:,0]
               starttime=time.time()
               crossinterpolate=Rbf(x,y,cross)
               xi=np.linspace(min(x),max(x),len(cross)*self.interpolationcoefficient)
               yi=np.linspace(min(y),max(y),len(cross)*self.interpolationcoefficient)
               cross=crossinterpolate(xi,yi)
               print('time_taken',time.time()-starttime)
               
               
           surfinder=np.diff(np.where(cross>surfacedepth,0,1),axis=0)
           sectionlen=len(np.split(self.fulldata, np.where(np.diff(self.fulldata[:,-self.spatialdim+1:],
                                                                   axis=0))[0]+1))*self.interpolationcoefficient
           surfinder=np.where((surfinder+1)%sectionlen==0,surfinder,0)
           surfindx=np.where(surfinder)[0]
           if plotyes:
               print(i)
               if (i<13):
                   label=self.dictforplot[0] + ' ' +self.dictforplot[i-8]
                   print(label)
               else:
                   label=self.dictforplot[1] + ' ' + self.dictforplot[i-11]
               x=x[surfindx]
               y=y[surfindx]
               plt.plot(x,y,label=label+' ('+self.label+')')
      
      if plotyes:
          plt.legend()
      else:
          return surfindx
          
   def luminocity(self,surfacedepth=2/3,itype=[10]):
       
       Lt=0
       surfindx=self.surface(surfacedepth=surfacedepth,itypes=itype,plotyes=False)
       
       for r,h,temp in zip(self.fulldata[surfindx,-1],self.fulldata[surfindx,0],self.fulldata[surfindx,1]):
           g_00=self.g00(r,h)
           L=g_00*abs(r)*2*cr.pi*7*temp**4*cr.pi**4/120
           Lt+=L
       return Lt
   
   def totflux(self,surfacedepth=2/3,itype=[10]):
       
       fluxt=0
       surfindx=self.surface(surfacedepth=surfacedepth,itypes=itype,plotyes=False)
       
       for r,h,temp in zip(self.fulldata[surfindx,-1],self.fulldata[surfindx,0],self.fulldata[surfindx,1]):
           g_00=self.g00(r,h)
           flux=g_00**0.5*abs(r)*2*np.pi*sp.zeta(3)*3/2*temp**3
           fluxt+=flux
       return fluxt
   
   def averageenergy(self,surfacedepth=2/3,itype=[10]):
       return self.luminocity(surfacedepth=surfacedepth,itype=itype)/self.totflux(surfacedepth=surfacedepth,itype=itype)
       
   def averageenergycompared(self,eflavour="Mu and Tau",surfdepthlimits=(0.01,2)):
       if eflavour=='Mu and Tau':
           itypes=[[12],[15]]
       elif eflavour=='Electron':
           itypes=[[10],[13]]
       elif eflavour=='Anti Electron':
           itypes=[[11],[14]]
           
       surfacedepths=np.linspace(*surfdepthlimits)
       compared=[]
       for surf in surfacedepths:
           Noninter=self.averageenergy(surfacedepth=surf,itype=itypes[0])
           Inter=self.averageenergy(surfacedepth=surf,itype=itypes[1])
           compared.append(Inter/Noninter)
       return compared,surfacedepths    
       
   def g00(self,rad,z,rs=0,a=1):
        '''
       

       Parameters
       ----------
       rad : Radius away from center of black hole.
       
       z : distance above plane of rotation.
       
       rs : Schwarschild radius
       
       a : black hole spin coefficient.

       Returns
       -------
       TYPE
           g00 GR modification.

        '''
        r=np.sqrt(rad**2+z**2)
        costheta=z/r
        epsi=r**2+a*costheta
        return 1-r*rs/epsi**2
    
   def effflux(self,itype,rs=0):
       surfindx=self.surface(itypes=itype,plotyes=False)
       r=self.fulldata[:,-1]
       temp=self.fulldata[:.1]
       #b=
       sflux=sp.zeta(3)*3/2*temp**3  # fermi dirac integral in the denominator 
       mod=(45*temp**5*sp.zeta(5))/2 #integral when the fermi dirac distribution is multiplied by E^2
       
       flux=0
       for r,h,temp in zip(self.fulldata[surfindx,-1],self.fulldata[surfindx,0],self.fulldata[surfindx,1]):
           flux+=1/(4*np.pi)
       return flux

def generate_data(datafile,param,datasetdict):
    
    '''

    Parameters
    ----------
    datafile : This is the unorganized data file we wish to organzie and calculate cross sections/optical
        depths for a given axis as defined in our datasetdictionary.
        
    param : This contains parameters that control how we account for nucleonic interaction
        
    datasetdict : This dictionary contains information about how our dataset is organized as well as
        information on which axis we calculate our neutrino emission from. Data 
        
        
        
    Returns
    -------
    dataout : This function returns data processed for use in our AccretionDisk Class. DataOut is organized
        as follows
        
        data[:,0]   : axis of neutrino emission
        data[:,1]   : temperature
        data[:,2]   : number density
        data[:,3]   : electron fraction
        data[:,4]   : non interacting electron neutrino cross sections
        data[:,5]   : non interacting electron anti neutrino cross sections
        data[:,6]   : non interacting mu and tau neutrino cross sections
        data[:,7]   : interacting electron neutrino cross sections
        data[:,8]   : interacting electron anti neutrino cross sections
        data[:,9]   : interacting mu and tau neutrino cross sections
        data[:,10]  : non interacting electron neutrino optical depth
        data[:,11]  : non interacting electron anti neutrino optical depth
        data[:,12]  : non interacting mu and tau neutrino optical depth
        data[:,13]  : interacting electron neutrino optical depth
        data[:,14]  : interacting electron anti neutrino optical depth
        data[:,15]  : interacting mu and tau neutrino optical depth
        data[:,16:] : other spatial dimensions
        


    '''
    
    datain=np.loadtxt(datafile)
    data=np.zeros((datain.shape[0],datain.shape[1]+6))
    data[:,0]=datain[:,datasetdict['axis_of_neutrino_emission']]
    data[:,1]=datain[:,datasetdict['temp_idx']]
    data[:,2]=datain[:,datasetdict['density_idx']]/1.6726219e15 #number density in natural units calculated from mass density in per grams cubed
    data[:,3]=datain[:,datasetdict['electronfrac_idx']]
    data[:,4:10]=gencrossdata(param,data[:,1],data[:,2],data[:,3])

    for i,idx in enumerate(datasetdict['addtional_spatial_dims']):
        data[:,-(i+1)]=datain[:,idx]
        
    dataout=genoptdata(data,len(datasetdict['addtional_spatial_dims'])+1)

    return dataout

def genoptdata(data,spatialdim):
    '''
    

    Parameters
    ----------
    data : function to calculate optical depth from the cross sections as organzied by generate_data
    
    spatialdim : number of spatial dimensions in simulation

    Returns
    -------
    dataout : fully organized dataset that will also be the output of generate_data.

    '''
        
    datas=(np.split(data, np.where(np.diff(data[:,-spatialdim+1:],axis=0))[0]+1))
    
    optdepths=[]
    newdats=[]
    for i,dat in enumerate(datas):
        dist=abs(data[1,0]-dat[0,0])
        dat=dat[dat[:,0].argsort(0)[::-1]]
        optdepths.append(np.cumsum(dat[:,4:10],0)*1e18*dist)
        newdats.append(dat)
    data=np.array(newdats).reshape(data.shape[0],-1)
    optdepths=np.array(optdepths).reshape(data.shape[0],-1)
    dataout=np.hstack((data[:,0:10],optdepths,data[:,-spatialdim+1:]))

    return dataout

    

def gencrossdata(param,temp,dens,electronfrac):
    '''
    

    Parameters
    ----------
    param : parameter controlling how the interacting cross section is calculated
        
    temp : Array of temperatures in simulation (MeV)
        
    dens : number density of nucleons (fm^{-3})
        
    electronfrac : Electron Fraction

    Returns
    -------
    crossdata : Ordered cross section data 

    '''
    rhop=dens*electronfrac
    rhon=dens-rhop
    
    #First calculate reactions we take to be non-interacting in this project

    #If we have an energy squared dependence, averages over the fermi dirac flux are analytic
    
    sflux=sp.zeta(3)*3/2*temp**3  # fermi dirac integral in the denominator 
    eavg=(7*temp**5*np.pi**4)/120 # integral when the fermi dirac distribution is multiplied by E
    mod=(45*temp**5*sp.zeta(5))/2 #integral when the fermi dirac distribution is multiplied by E^2
    avge=eavg/sflux #average energy for use in annihilation reactions
    
    #electron collisions hav more complex energy depedences 
    vhit=[]
    avhit=[]
    for t in temp:
        
        f=lambda e : cr.veve(e)*cr.fdflux(e,t)
        vhit.append(integ.quad(f,0,80)[0])
        
        f=lambda e : cr.aveve(e)*cr.fdflux(e,t)
        avhit.append(integ.quad(f,0,80)[0])

    vhit=np.array(vhit)
    avhit=np.array(avhit)

    spontreac=vhit*ed.edense.edens(rhop,temp,False)/sflux
    aspontreac=avhit*ed.edense.edens(rhop,temp,False)/sflux
    anile=ed.edense.edens(rhop,temp,True)*cr.veanil(1,avge)*eavg/sflux
    anilx=ed.edense.edens(rhop,temp,True)*cr.vanil(1,avge)*eavg/sflux
    vecharged=(cr.neutn(1)*rhon)*mod/sflux
    avecharged=cr.antineutp(1)*rhop*mod/sflux
    neutralreac=(cr.neutap(1)*rhop+cr.neutan(1)*rhon)*mod/sflux
    
    # Now lets evaluate the interacting neutral current reactions of interest
    if param['method']=='viri':
        
        '''
        Here we calculate virial crosssections with options to use the pure virial method or the fit for axial 
        We also allow for the turning on and off of the axial or vector structure factors.
        '''
        
        paramv=param['viriparam']
        SV=1
        SA=1
        
        
        if paramv['SVon']==True:
            SV=np.zeros(temp.shape)
            i=0
            for t,y,rho in zip(temp,electronfrac,dens):
                SV[i]=st.sv(t,y,rho)
                i+=1

        if paramv['SAfit']==True:
            SA=st.sa_f(temp,electronfrac,dens)
        elif paramv['SAlow']==True:
            SA=np.zeros(temp.shape)
            i=0
            for t,y,rho in zip(temp,electronfrac,dens):
                SA[i]=st.Sa_loworder(t,y,rho)
                i+=1
        neutralreac_inter=(cr.neutapi(1,SA,SV)*rhop+cr.neutani(1,SA,SV)*rhon)*mod/sflux
        
    elif param['method']=='RPA':
        
        paramRPA=param['paramRPA']
        
        functionalname=paramRPA['functional name'] # name of energy density functional to use in calculations
        ftype=paramRPA['type'] # this allows the turning on and off of certain types oforces (eg. Tensorial forces)
        readforce=paramRPA['readforce'] # are the parameters provided in terms of the residual force or energy density functional
        
        neutralreac_inter=np.zeros(temp.shape)
        i=0
        for t,y,rho in zip(temp,electronfrac,dens):
            os.chdir('RPAStructure')
            neutralreac_inter[i]=RPA.anmfp(rho,y,t,functionalname,ftype,readforce,1)**-1**2*mod[i]/sflux[i]/1e15
            os.chdir('..')
            i+=1
            print(i,temp.shape,neutralreac_inter[i])
                
    #Now we need to put the reactions together and find the total cross sections for different species
    
    vecross=neutralreac+spontreac+anile+vecharged # electron neutrinos can induce isospin flips and experience all other interactions
    avecross=neutralreac+aspontreac+anile+avecharged # anti-electron neutrinos can induce isospin and experience all other interactions
    vxcross=neutralreac+anilx # other flavors cannot induce isospin flips and do not interact with electrons
    #interacting cross sections of different species
    
    inter_vecross=neutralreac_inter+spontreac+anile+vecharged 
    inter_avecross=neutralreac_inter+spontreac+anile+avecharged
    inter_vxcross=neutralreac_inter+anilx 
    
    crossdata=np.transpose(np.array([vecross,avecross,vxcross,inter_vecross,inter_avecross,inter_vxcross]))
    return crossdata
    

    
    
    

        
    



        
param={
       'preprocessed_crosssections':True, #If you are running a new dataset or using a new method of interaction  set this to False, the other parameter control how interaction is accounted for
       #a standard preprocessed parameter set is can be accessed in the datadictionaries module as dicts.preparam
       
       
       'crossparam':
           {
               'method':'viri', #either RPA or viri
               'viriparam':
                   {
                       'SVon':True, #Turns the Virial Vector Structure factor on
                       'SAfit':False, #Turns the Virial RPA Fit Axial structure factor on
                       'SAlow':True, #Turns the true Virial Axial structure factor on
                       
                       },
                'paramRPA' :
                    {
                        'functional name' : 'BSK20', #see .param files in RPAStructure directory for possible inputs
                        'type' : 'NONE', 
                          # NONE -> complete functional
                          # SORB -> spin orbit a zero
                          # TENS -> tensor only at zero
                          # SOTE -> spin-orbit and tensor to zero
                        'readforce' : True # make sure whichever force or functional you select is in the appropriate .param file
                        }
                   
               
           }
       }
    

if __name__=='__main__': 
    print('time to process more data')
    acd=AccetionDisk(dicts.preparam,'data/NL3_t2.5ms_Y0/f0/crossoptdata.dat',dicts.NL3_datasetdict)
    acd.surface(itypes=[12,15],interpolate=True)
    
        
    

