#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 21:30:56 2020

@author: rajananderson
"""

import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import structfact as st
import pandas as pd
import sys
    
file='NL3_t2.5ms_Y0'
file='time3500Y0'
#file='M3a0.8t60'
file='data/'+file+'/'+file+'.dat'
data=np.loadtxt(file)
r=data[:,0]
z=data[:,2]

datas=np.split(data,len(set(r)))
xlab=sorted(list(set(data[:,2])))
ylab=sorted(list(set(data[:,0])))
tempidx=3
density_dx=5
yidx=4
heatmapmatrix=np.zeros((3,len(xlab),len(ylab)))
for i,data in enumerate(datas):
    for j,dat in enumerate(data):
        x=j
        y=i
            
        heatmapmatrix[0,x,y]=1-st.sv(dat[tempidx], dat[yidx], dat[density_dx]/1.6726219e15)
        #heatmapmatrix[1,-x,y]=1-st.sa_f(dat[tempidx], dat[yidx], dat[density_dx]/1.6726219e15)
        #heatmapmatrix[2,-x,y]=1-st.Sa_loworder(dat[tempidx], dat[yidx], dat[density_dx]/1.6726219e15)
xx=0
vertical=np.zeros([201,201*3])
horiz=np.hstack([np.zeros([201,201]),heatmapmatrix[xx,:,:],np.zeros([201,201])])
newdata=np.vstack([vertical,horiz,vertical])
xlab=np.arange(-42.5,42.5,85/603)
ylab=np.arange(-42.5,42.5,85/603)
#newdata=heatmapmatrix[xx,:,:]
n=3 #number of sig figs in the label
xlab=['{:g}'.format(float('{:.{p}g}'.format(i, p=n))) for i in xlab][::-1]
ylab=['{:g}'.format(float('{:.{p}g}'.format(i, p=n))) for i in ylab]
plt.figure(dpi=700)
plt.title('Vector Virial')
ax=sns.heatmap(pd.DataFrame(newdata,index=xlab,columns=ylab))
#ax.set_yticks([0,50,100,150,200,250,300,350,400,450,500])
#ax.set_xticks([0,63,163,263,363,463,563,663,763,863,963])
#plt.ylabel('Distance from Plane of Acretion Disk')
plt.xlabel('Distance from Center of Accretion Disk')
plt.show()
