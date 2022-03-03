import main
import numpy as np
import scipy.stats as stat
import matplotlib.pyplot as plt
import scipy.special as sp
import time
import seaborn as sns
import pandas as pd
import math

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




def functional(rho,energy,temp,y):
    return main.main.anmfp(rho,y,temp,'BARE','NONE',True,energy)
    

forcename=['SLY4','SLY4','BSK21']
temp=0.69940483570098877#np.linspace(0.1,15,num=10)
rho=188168.58617223453/1.6726219e15#np.logspace(-10,-3,num=2)
y=0.258980006


rho=3.90353271e-06
e=4
t=2
y=0.3

a=functional(rho,e,t,y)*1e15
b=1/(neutan(e)*(1-y)*rho+neutap(e)*y*rho)
print(a**-1*1e18,b**-1*1e18,a/b)

