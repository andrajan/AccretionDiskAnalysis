#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 21:30:56 2020

@author: rajananderson
"""


import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import structfact as st
    

data=np.loadtxt('M3a0.8t20.dat')
r=data[:,0]
z=data[:,1]
temp=data[:,2]
rho=data[:,3]
y=data[:,4]


surfdat=np.loadtxt('esurfdatan.dat')
rad=surfdat[:,0]
height=surfdat[:,1]-1

datas=np.split(data,len(set(r)))



for i,dat in enumerate(datas):
    strutfit_tot=[]
    strut_V=[]
    strutfit_A=[]
    strutlow_A=[]
    print(i)
    for j,data in enumerate(dat):
        
        strutfit_tot.append(1-st.stot_f(data[2], data[4], data[3]/1.6726219e15))
        strut_V.append(1-st.sv(data[2], data[4], data[3]/1.6726219e15))
        strutfit_A.append(1-st.sa_f(data[2], data[4], data[3]/1.6726219e15))
        strutlow_A.append(1-st.Sa_loworder(data[2], data[4], data[3]/1.6726219e15))
      
    if i==0:
        astrutfit_tot=np.array(strutfit_tot)
        astrutfit_A=np.array(strutfit_A)
        astrut_V=np.array(strut_V)
        astrutlow_A=np.array(strutlow_A)
    else:
        astrutfit_tot=np.vstack([astrutfit_tot,np.array(strutfit_tot)])
        astrutfit_A=np.vstack([astrutfit_A,np.array(strutfit_A)])
        astrut_V=np.vstack([astrut_V,np.array(strut_V)])
        astrutlow_A=np.vstack([strutlow_A,np.array(strutlow_A)])
    
file=open('strutfit_A.dat','w')
np.savetxt(file,astrutfit_A)
file.close()

file=open('STOTFIT.dat','w')
np.savetxt(file,astrutfit_tot)
file.close()

file=open('strutlow_A.dat','w')
np.savetxt(file,astrutlow_A)
file.close()

file=open('strutfit_V.dat','w')
np.savetxt(file,astrut_V)
file.close()


strutfit_tot=np.loadtxt('StructureFactors_Limited.dat')
plt.hist(strutfit_tot.tolist())
plt.show()
'''
strutfit_tot=np.loadtxt('StructureFactors_Limited.dat')
strutfit_tot=np.array(strutfit_tot)

strutfit_tot=[stru for stru in strutfit_tot if stru[0]>0.002]
strutfit_tot=np.array(strutfit_tot)
plt.scatter(strutfit_tot[:,2],strutfit_tot[:,1])
plt.title('Diminished Search Area after Eliminating Smallest Bin from Histogram')
plt.ylabel('Distance from Plane of Acretion Disk')
plt.xlabel('Distance from Center of Accretion Disk')
plt.show()
'''
ylab=list(range(0,251,25))
xlab=[18.5]
for xl in range(50,550,50):
    xlab.append(xl)
    
xlab=[str(x) for x in xlab]
ylab=[str(y) for y in ylab]



for yl in range((int(50-18.5)*2),int((500-18.5)*2)):
    ylab.append(yl)

plt.title('Axial Component from Virial')
ax=sns.heatmap(astrutlow_A,xticklabels=xlab,yticklabels=ylab)
ax.set_yticks([0,50,100,150,200,250,300,350,400,450,500])
ax.set_xticks([0,63,163,263,363,463,563,663,763,863,963])
plt.ylabel('Distance from Plane of Acretion Disk')
plt.xlabel('Distance from Center of Accretion Disk')
plt.show()