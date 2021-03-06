'''
This module can be used to verify that the our RPA cross sections still demonstrate an energyu squared
dependence and that bare hartree fock response cross sections are close to the non-interacting elastic case.

'''


import main
import numpy as np
import scipy.stats as stat
import matplotlib.pyplot as plt
import scipy.special as sp
import time
import seaborn as sns
import pandas as pd
import math
import sys

mn=939.56563
hbarc=197.32697
cv=1/2+2*0.23
ca=0.5
me=0.511
ga=1.26
delt=1.293
gf=1.166379e-11 
pi=3.141592
phi0=4*(gf*me)**2*(hbarc)**2/pi

def neutap(e):  # neutral p + v --> p + v
    neutap=phi0*((cv-1)**2+3*ga**2*(ca-1)**2)*e**2/(4*me**2)
    return neutap

def neutan(e): # neutral n + v --> n + v
    neutan=phi0*(1+3*ga**2)*e**2/(16.*me**2)
    return neutan


def functional(rho,energy,temp,y,force='BARE',mode='None',forceorfunctional=True):
    '''
    Fortran function but with a bunch of default options

    '''
    return main.main.anmfp(rho,y,temp,force,mode,forceorfunctional,energy)
    

def testingnoninteractingsimilar(rho,e,t,y):
    '''
    

    Parameters
    ----------
    rho : density(fm^-3).
    e : Energy (mev) of incoming neutrino.
    t : Temperature (MeV).
    y : Electron Fraction.

    Returns
    -------
        Non interacting elastic over hartree fock cropss section.

    '''
    starttime=time.time()
    a=functional(rho,e,t,y)*1e15
    b=1/(neutan(e)*(1-y)*rho+neutap(e)*y*rho)
    endtime=time.time()
    print(a**-1*1e18,b**-1*1e18,'how similar :',a/b)
    print('{} seconds taken with current cores allocated'.format(endtime-starttime)' for one response function to be calculated')
    return a/b


def r2valscompared(temp,rho,y):
    '''
    Plots heatmap of r^2 values to determine suitability of assuming a energy squared dependence

    Parameters
    ----------
    temp : a list of temperatures(MeV) for our grid
    
    rho : a list of densities for our grid (fm^-3.)
    
    y : A single Electron Fraction.


    '''
    n=3 #number of sig figs in the label
    templabel=['{:g}'.format(float('{:.{p}g}'.format(i, p=n))) for i in temp][::-1]
    rholabel=['{:g}'.format(float('{:.{p}g}'.format(i, p=n))) for i in rho]
    
    rvalheat=np.zeros((temp.shape[0],rho.shape[0]))
    for it,t in enumerate(temp):
        sflux=sp.zeta(3)*3/2*t**3
        eavg=(7*t**4*pi**4)/120
        avge=eavg/sflux
        energy=np.array([0.5,1,1.5,2])*avge
        for irho,rh in enumerate(rho):
            test=[]
            fit_tester=[]
            for e in energy:
                fit_tester.append(functional(rh,e,t,y)*rh*e**3)
            r2val=stat.linregress(energy,fit_tester)[2]**2
            rvalheat[-(it+1),irho]=r2val
    
    rvalheat=pd.DataFrame(rvalheat,index=templabel,columns=rholabel)
    fig=plt.figure(dpi=600)
    fig=sns.heatmap(rvalheat,cmap="mako")
    plt.xlabel('Density ($fm^{-3}$)')
    plt.ylabel('Temperature (MeV)')
    plt.title('$R^{2}$ Values for $1/E^{2}$ Dependence'+' (y={})'.format(y))
    #plt.savefig('y{}_cmapSLY4'.format(y)) # unsupported apparently
    plt.tight_layout()

if __name__=='__main__':
    
    print(testingnoninteractingsimilar(1e-7,3,3,0.4))
    
    sys.exit()
    temp=np.linspace(0.1,15,num=2)
    rho=np.logspace(-10,-3,num=2)
    r2valscompared(temp,rho,0.2)
    
    endtime=time.time()

