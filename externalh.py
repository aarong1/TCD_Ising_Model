#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: aarongorman
"""

import numpy as np
import matplotlib.pyplot as plt
import functions as IS  #module containing our packaged functions

N=16
h=0  #vary the external field
Tc=2.269 #set around the literature theoretical value of critical temp

#initialising thermodynamic statistics to be collected
magnetisation  = []
energy = []

energyh=list(np.zeros(7))
magnetisationh=list(np.zeros(7))
exh=[]
t=[]

for index,h in enumerate((1,2,5,10,15)):
    print('The external field for this run is '+str(h))
    print('The index is '+str(index))
    exh.append(h)
    
    for m in range(0,60):
        print(m)
        T = m / 30.0+1
        print('Temp, T is '+str(T))
        t.append(T)

        mag,magsq,eng,engsq,S= IS.Metropolis(N,T)

        energy.append(eng/(800.0*N*N))
        magnetisation.append(np.abs(mag)/(800.0*N*N))
    energyh[index]=energy
    magnetisationh[index]=magnetisation

#plotting magnetisation for a range of applied external fields
f=plt.figure(2)
plt.xlabel("Temperature ($k_b$ K)")
plt.ylabel("Magnetization per spin for applied h")
plt.plot(t, magnetisation[1], '.', color="orange",label='h=1')
plt.plot(t, magnetisation[2], '.', color="blue",label='h=2')
plt.plot(t, magnetisation[5], '.', color="green",label='h=5')

f=plt.figure(3)
plt.plot(exh,energyh[90])

f=plt.figure(4)
plt.plot(exh,magnetisationh[90])
