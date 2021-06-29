#!/usr/bin/env python3 
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 14:02:49 2020

@author: rajananderson

This module calculates cross sections for the virial approaches and non-interacting approaches.
Other unaltered cross sections are also present here.

"""
import math as m
import structfact as st
import edense as ed
import numpy as np
import scipy.integrate as integ
import scipy.special as sp

print(ed.__doc__)

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


def neutn(e):    # charged (n + v -> p + e)
    wm=1+1.1*e/mn
    neutn=phi0/(4*me**2)*(1+3*ga**2)*(e+delt)**2*(1-(me/(e+delt))**2)**(-1/2)*wm
    return neutn

def antineutp(e):  # p + v --> n + v
    wm=1-7.1*e/mn
    if (wm<=0):
        antineutp=0
    elif((1-(me/(e-delt))**2)<0):
        antineutp=0
    else:
        antineutp=phi0/(4*me**2)*(1+3*ga**2)*(e-delt)**2*(1-(me/(e-delt))**2)**(-1/2)*wm
    return antineutp

def neutap(e):  # neutral p + v --> p + v
    neutap=phi0*((cv-1)**2+3*ga**2*(ca-1)**2)*e**2/(4*me**2)
    return neutap

def neutan(e): # neutral n + v --> n + v
    neutan=phi0*(1+3*ga**2)*e**2/(16.*me**2)
    return neutan

def neutapi(e,SA,SV): # neutral p + v --> p + v
    neutap=phi0*((cv-1)**2*SV+SA*3*ga**2*(ca-1)**2)*e**2/(4*me**2)
    return neutap

def neutani(e,SA,SV): # neutral n + v --> n + v
    neutan=phi0*(SV+3*SA*ga**2)*e**2/(16.*me**2)
    return neutan


def veve(e): # electron collision ve + e --> ve + e
    te=2*e**2/(me+2*e)
    veve=phi0/(8.*me)*((cv+ca)**2*te+(cv-ca)**2*(te**3/(3.*e**2)-te**2/e+te)-(cv**2-ca**2)*me*te**2/(2.*e**2))
    return veve

def aveve(e): 
    te=2*e**2/(me+2*e) # anti neutrino electron collision
    aveve=phi0/(8.*me)*((cv-ca)**2*te+(cv+ca)**2*(te**3/(3.*e**2)-te**2/e+te)-(cv**2-ca**2)*me*te**2/(2.*e**2))
    return aveve

def fdflux(e,temp): # fermi dirac flux
    fdflux=e**2/(np.exp(e/temp)+1)
    return fdflux

def veanil(e,avge):     # electron neutrino anilhilation  ve + bar(ve) --> e + bar(e)
    ke=(1+4*0.23**2.+8*0.23**2.)/(6*pi)
    vanil=4/3*ke*phi0*e*avge
    return vanil

def vanil(e,avge):   # other neutrino annhilation vx + bar(vx) --> x + bar(x)
    ke=(1-4*0.23**2.+8*0.23**2.)/(6*pi)
    vanil=4/3*ke*phi0*e*avge
    return vanil


def vecross(temp,dens,y,inter=False,SA=1,SV=1): 
    '''

    Parameters
    ----------
    inter : If inter is selected the axial fit structure factor and the vector virial structure factor will be turned on.

    SA : Provides option ato input directly the axial structure factor.
    SV : Provides option ato input directly the vwector structure factor.

    Returns
    -------
    cross : elecrtron neutrino cross section.

    '''

    rhop=dens*y
    rhon=dens-rhop
    
    if inter:
        SA=st.sa_f(temp,y,dens)
        SV=st.sv(temp, y, dens)
    
    
    sflux=sp.zeta(3)*3/2*temp**3
    eavg=(7*temp**4*pi**4)/120
    avge=eavg/sflux
    
    
    f=lambda e : neutapi(e,SA,SV)*fdflux(e,temp)
    ven=integ.quad(f,0,80)[0]
    f=lambda e : veve(e)*fdflux(e,temp)
    vhit=integ.quad(f,0,80)[0]

    mod=(45*temp**5*sp.zeta(5))/2

    charged=ven*rhon/sflux
    anil=ed.edense.edens(rhop,temp,True)*vanil(1,avge)*eavg/sflux
    neutralreac=(neutapi(1,SA,SV)*rhop+neutani(1,SA,SV)*rhon)*mod/sflux
    spontreac=vhit*ed.edense.edens(rhop,temp,False)/sflux
    cross=neutralreac+charged+spontreac+anil
    

    return cross

def avecross(temp,dens,y,inter=False,SA=1,SV=1):
    '''

    Parameters
    ----------
    inter : If inter is selected the axial fit structure factor and the vector virial structure factor will be turned on.

    SA : Provides option ato input directly the axial structure factor.
    SV : Provides option ato input directly the vwector structure factor.

    Returns
    -------
    cross : anti-elecrtron neutrino cross section.

    '''
    rhop=dens*y
    rhon=dens-rhop
    
    if inter:
        SA=st.sa_f(temp,y,dens)
        SV=st.sv(temp, y, dens)
        
    sflux=sp.zeta(3)*3/2*temp**3
    eavg=(7*temp**4*pi**4)/120
    avge=eavg/sflux
    
    
    f=lambda e : antineutp(e)*fdflux(e,temp)
    vep=integ.quad(f,0,80)[0]
    f=lambda e : aveve(e)*fdflux(e,temp)
    vhit=integ.quad(f,0,80)[0]

    mod=(45*temp**5*sp.zeta(5))/2


    charged=vep*rhon/sflux
    anil=ed.edense.edens(rhop,temp,True)*vanil(1,avge)*eavg/sflux

    neutralreac=(neutapi(1,SA,SV)*rhop+neutani(1,SA,SV)*rhon)*mod/sflux

    spontreac=vhit*ed.edense.edens(rhop,temp,False)/sflux
    cross=neutralreac+charged+spontreac+anil
    

    return cross

def neutralcross(temp,dens,y,inter,SA=1,SV=1,ST=1):
    '''

    Parameters
    ----------
    inter : If inter is selected the axial fit structure factor and the vector virial structure factor will be turned on.

    SA : Provides option ato input directly the axial structure factor.
    SV : Provides option ato input directly the vwector structure factor.
    stru: Provides option to manually input total structure factor

    Returns
    -------
    cross : Cross section of only neutral current reactions

    '''

    rhop=dens*y
    rhon=dens-rhop
    stru=1
    
    if inter:
        SA=st.sa_f(temp,y,dens)
        SV=st.sv(temp, y, dens)
    
    
    sflux=sp.zeta(3)*3/2*temp**3
    mod=(45*temp**5*sp.zeta(5))/2

    neutralreac=(neutapi(1,SA,SV)*rhop+neutani(1,SA,SV)*rhon)*mod/sflux*stru
    
    return neutralreac

def vcross(temp,dens,y,inter=True,SA=1,SV=1):
    '''

    Parameters
    ----------
    inter : If inter is selected the axial fit structure factor and the vector virial structure factor will be turned on.

    SA : Provides option ato input directly the axial structure factor.
    SV : Provides option ato input directly the vwector structure factor.
    stru: Provides option to manually input total structure factor

    Returns
    -------
    cross : Mu and Tau neutrino cross section

    '''

    rhop=dens*y
    rhon=dens-rhop
    
    if inter:
        SA=st.sa_f(temp,y,dens)
        SV=st.sv(temp, y, dens)
    
        
    sflux=sp.zeta(3)*3/2*temp**3
    eavg=(7*temp**4*pi**4)/120
    avge=eavg/sflux
    
    
    f=lambda e : aveve(e)*fdflux(e,temp)
    vhit=integ.quad(f,0,80)[0]

    mod=(45*temp**5*sp.zeta(5))/2


    anil=ed.edense.edens(rhop,temp,True)*vanil(1,avge)*eavg/sflux

    neutralreac=(neutapi(1,SA,SV)*rhop+neutani(1,SA,SV)*rhon)*mod/sflux

    spontreac=vhit*ed.edense.edens(rhop,temp,False)/sflux
    cross=neutralreac+spontreac+anil
    
    return cross


def crossfit(temp,dens,y):
    rhop=dens*y
    rhon=dens-rhop
    
    stru=1
    SA=1
    SV=1
        
    
    f=lambda e : neutapi(e,SA,SV)*fdflux(e,temp)
    vp=integ.quad(f,0,800)[0]
    f=lambda e : neutani(e,SA,SV)*fdflux(e,temp)
    vn=integ.quad(f,0,800)[0]
    f=lambda e : fdflux(e,temp)
    sflux1=integ.quad(f,0,800)[0]

    sflux=sp.zeta(3)*3/2*temp**3
    mod=(45*temp**5*sp.zeta(5))/2
    trial=(neutapi(1,SA,SV)*rhop+neutani(1,SA,SV)*rhon)*mod/sflux
    neutralreac=(vp*rhop+vn*rhon)/sflux1*stru
    cross=neutralreac

    return cross,trial,mod/sflux
    



if __name__=='__main__':
    temp=200
    dens=8.465704652079468e-6
    y=0.448302748
    print(ed.edense.edens(dens,y,temp))
    #11.1283928 8.465704652079468e-06 0.448302748
    #print(neutralcross(temp,dens,y,True)/neutralcross(temp,dens,y,False),neutralcross(temp,dens,y,False))


