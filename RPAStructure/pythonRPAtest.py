#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 11 02:58:41 2021

@author: rajananderson
"""
import numpy as np
import scipy.integrate as integ
import matplotlib.pyplot as plt

def tp0(msn,omega,q,mun,T):
    Q=(msn*omega/q - q/2)/(2*msn)**0.5

    aux1=np.exp(-(Q**2-mun)/T)
    aux2=np.exp(-(Q**2-mun+omega)/T)

    beta0IM=T*(msn**2)/(4*np.pi*q)*np.log((1+aux1)/(1.+aux2)) 
    
    return beta0IM


def tp0R(msn,omega,q,mun,T,s):
    Q=(msn*omega/q - q/2)/(2*msn)**0.5

    aux1=np.exp(-((Q+s)**2-mun)/T)
    aux2=np.exp(-((Q-s)**2-mun)/T)

    beta0R=T*(msn**2)/(4*np.pi*q)*np.log((1+aux1)/(1.+aux2)) 
    
    
    return beta0R

def beta0R(msn,omega,q,mun,T):
    f=lambda s: tp0R(msn,-omega,q,mun,T,s)/s
    bR=integ.quad(f,0,np.inf)[0]
    return bR*T*(msn**2)/(4*np.pi**2*q)

def beta(msn,omega,q,mun,T):
    return complex(beta0R(msn,omega,q,mun,T), tp0(msn,omega,q,mun,T))

energy=6

msn=938.91818294668667
T=5
q=3*T
mu=[-35,-30,-20,-10,-5,0]

def densfrommu(p,mun,T,msn):
    return p**3/(2*np.pi)**3*1/(1+np.exp(1/T*(p**2/(2*msn)-mun)))/(197.3)**3

mun=0

q0=np.linspace(-7*energy,energy,num=100)

funcs=[]
funcshf=[]
'''
for omega in q0:
    #funcs.append((beta(msn,omega,q,mun,T)/(1+1.75e-5*beta(msn,omega,q,mun,T))/(1+np.exp(-omega/T))).imag)
    funcshf.append((beta(msn,omega,q,mun,T)/(1+4.5e-5*beta(msn,omega,q,mun,T))/(1+np.exp(-omega/T))).imag)

    
'''
#plt.plot(q0,funcs)
#plt.plot(q0,funcshf)
mun=-21


print('response', (beta(msn,2,q,mun,T)/(1+1.75e-5*beta(msn,2,q,mun,T))).imag)
truedens=[]
limits=energy
for mun in mu:
    densfrommun = lambda p : densfrommu(p,mun,T,msn)
    trudens=integ.quad(densfrommun,0,np.inf)[0]
    truedens.append(trudens)
    print(trudens,mun)
    
    
    f= lambda x : (beta(msn,x,q,mun,T)/(1+np.exp(-x/T))).imag/197.3**3
    integral=integ.quad(f,0,np.inf)
    integral=-integral[0]/np.pi
    hf=integral
    
    funcs.append(hf)
'''
    f= lambda x : (beta(msn,x,q,mun,T)/(1+1.75e-5*beta(msn,x,q,mun,T))/(1+np.exp(-x/T))).imag
    integral=integ.quad(f,0,np.inf,limit=10000)
    integral=-integral[0]/np.pi
    hf=integral
    
    funcshf.append(hf)
    '''

plt.plot(truedens,funcs)
#plt.plot(truedens,funcshf)
print(sum(funcs))
print(funcshf,funcs)
#   4.8007852268265054        938.91818294668667        5.5253536111378230        4.0000000000000000       -4.9168821862619101       0.30936677262698237        0.0000000000000000     


