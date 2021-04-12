#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 22:03:38 2020

@author: rajananderson
"""
import numpy as np
import sympy
from sympy import Symbol, solve, nsolve

def sa_f(temp,y,dens):
    a0,b0,c0,d0=920,3.05,6140,1.5e13
    a=a0*dens*(1.-y+y**2.)/temp**1.22
    b=b0/temp**0.75
    c=c0*dens*y*(1.-y)/temp**0.5+d0*dens**4./temp**6.
    #print(c,d0*dens**4./temp**6.,dens)
    sa=1/(1+a*(1+b*np.exp(-c)))
    return sa


def Sa_loworder(temp,y,dens):
    rhon,rhop=seperatedens(y,dens)
    lamb=(2.*np.pi/938.27208816/temp)**0.5*197.33

    
        
    bn=(0.3608245735172147 - 0.07166371507427416*temp - 0.001130419243529185*temp**2
                     + 9.249283273497513e-6*temp**3 + 0.07814148523187128*np.log(temp) + 
                     0.02564700493123892*temp*np.log(temp))

    bpn= (0.05866576478320349 + 2.1388457038062465*np.exp(2.22/temp) - 0.34531025628632056*temp - 
                       0.006483144053841453*temp**2 + 0.00006892075911127388*temp**3 + 
                       0.13340377588640795*np.log(temp) + 0.12602735446104169*temp*np.log(temp))
    
    
    fugn=fuga(bn,bpn,rhop,rhon,lamb)
    fugp=fuga(bn,bpn,rhon,rhop,lamb)
    fugmask=[zn>=0 and zp>=0 for zn,zp in zip(fugn,fugp)]
    fugn=fugn[fugmask]
    fugp=fugp[fugmask]

    npp=2/lamb**3*(fugn+fugn**2*bn+2*bpn*fugn*fugp)
    nn=2/lamb**3*(fugp+fugp**2*bn+2*bpn*fugn*fugp)
    choop=[rhop-npr for npr in npp ]
    choon=[rhon-nne for nne in nn]
    
    fugn=fugn[choon.index(min(choon))]
    fugp=fugp[choop.index(min(choop))]


    ba=-0.691863 +  0.0449494/temp+0.00888447*temp -0.000164002*temp**2
    bnp=-0.586369 + 0.801285*np.exp(2.13342/temp)+ 0.000419746*temp

    return (1+4/lamb**3*((fugn**2+fugp**2)*ba-2*fugn*fugp*bnp)/dens)


def seperatedens(y,dens):
    rhop=dens*y
    rhon=dens-rhop
    return rhon,rhop

def fuga(bn,bpn,rho1,rho2,lamb):
    a=lamb**6*rho1**2*bn/8/bpn**2
    b=lamb**3*rho1/4/bpn-lamb**3*rho1*bn/2/bpn**2
    c=-1/2/bpn+lamb**3*rho1/2-lamb**3*rho1*bn**2/bpn**2+bn/2/bpn**2-rho2*lamb**3/2
    d=-bn/bpn-1+2*bn**2/bpn**2
    e=-2*bn+2*bn**3/bpn**2
    coef=[e,d,c,b,a]
    fug=np.roots(coef)
    fug=fug[fug>=0]
    if len(fug)==0:
        print('no positive fugacity solution available')
    return fug


def sv(temp,y,dens):
    rhon,rhop=seperatedens(y,dens)
    cnv=0.5
    cpv=0.04
    
    lamb=(2.*np.pi/938.27208816/temp)**0.5*197.33
    
    bn=(0.3608245735172147 - 0.07166371507427416*temp - 0.001130419243529185*temp**2
                     + 9.249283273497513e-6*temp**3 + 0.07814148523187128*np.log(temp) + 
                     0.02564700493123892*temp*np.log(temp))

    bpn= (0.05866576478320349 + 2.1388457038062465*np.exp(2.22/temp) - 0.34531025628632056*temp - 
                       0.006483144053841453*temp**2 + 0.00006892075911127388*temp**3 + 
                       0.13340377588640795*np.log(temp) + 0.12602735446104169*temp*np.log(temp))
    
    
    fugp=fuga(bn,bpn,rhop,rhon,lamb)
    fugn=(lamb**3*rhop/2-fugp-2*fugp**2*bn)/(2*fugp*bpn)
    #fugmask=[zn>=0 and zp>=0 for zn,zp in zip(fugn,fugp)]
    #fugn=fugn[fugmask]
    #fugp=fugp[fugmask]

    #print(fugp,bn,bpn,fugn,temp,y,dens)
    npp=[2/lamb**3*(fn+fn**2*bn+2*bpn*fn*fp) for fn,fp in zip(fugn,fugp)]
    nn=[2/lamb**3*(fp+fp**2*bn+2*bpn*fn*fp) for fn,fp in zip(fugn,fugp)]
    choop=[rhop-npr for npr in npp]
    choon=[rhon-nne for nne in nn]
    fugn=fugn[choon.index(min(choon))]
    fugp=fugp[choop.index(min(choop))]
    
    print(np.log(fugn)*temp)
    svc=(1.+4./lamb**3.*(cnv**2.*fugn**2.*bn+2.*cpv*cnv*fugn*fugp*bpn+cpv**2.*fugp**2.*bn)
    /(cnv**2.*rhon+cpv**2.*rhop))
    return svc
'''
def fugax(ba,bapn,rhon,rhop,lamb):
    rhonp=0.5*rhon
    rhonm=rhon-rhonp
    rhopp=0.5*rhop
    rhopm=rhop-rhopp
    
    fugnp=Symbol('fugnp')
    fugnm=Symbol('fugnm')
    fugpp=Symbol('fugpp')
    fugpm=Symbol('fugpm')
    bap=Symbol('bap')
    bam=Symbol('bam')
    
    rhonp=1/lamb**3*(fugnp+2*fugnp**2*bap+fugnp*(bam*fugnm+))
    
'''

def stot_f(temp,y,dens):
    ga=1.93
    stot=(5*ga**2*sa_f(temp,y,dens)+(1-y)*sv(temp,y,dens))/(5*ga**2+1-y)
    return stot

if __name__=='__main__': 
    temp=[0.5,1,2,4,6,10,14]
    for t in temp:
        print(sv(t,0.0000000001, 1e-10))

