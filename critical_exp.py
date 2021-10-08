#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is used to find the critical exponents by the finite scaling method
each is found by varying the log of the thermodynamic quantity in line 32
i.e magnetisation -> specificH or susceptibility
the result prints to screen
the error in fitting is quite high => take the average of multiple simulations
@author: aarongorman
"""

import numpy as np
import matplotlib.pyplot as plt
import functions as IS

T=2.5
h= 0  #no external magnetic field

#init data structures
magnetisation  = []
susceptibility = []
energy =[]
specificH=[]
n=[]

#for different sized (square)lattices
for N in (20,30,50,70):
    print(N)
    n.append(N)
    mag,magsq,eng,engsq,S = IS.Metropolis(N,T)
    energy.append(eng/(800.0*N*N))
    magnetisation.append(np.abs(mag)/(800.0*N*N))
    susceptibility.append((magsq/(800.0*N*N) - mag*mag/(800.0*800*N*N))/T)
    specificH.append((engsq/(800.0*N*N) - eng*eng/(800.0*800*N*N))/T**2)

#fitting best line
x,y=np.log(n),np.log(magnetisation)
fig = plt.figure()
p=np.polyfit(x,y,1)
slope,intercept=p

#plotting best fit line
best_fit=slope*x+intercept
plt.plot(x,best_fit, label='best fit line')
plt.plot(x,y,'o')
print(p)
#ax = fig.add_subplot(1, 1, 1)

fig.legend(loc='upper right')
#recalculating, by hard coding to find the fitting errors
no=len(x)
s_x=sum(x)
s_y=sum(y)
s_xx=sum(x**2)
s_xy=sum(x*y)
denom=no*s_xx-s_x**2
c=(s_xx*s_y-s_x*s_xy)/denom
m=(no*s_xy-s_x*s_y)/denom

sigma = np.sqrt(sum(( y - ( c + m*x ) )**2 )/(no -2))
sigma_c = np.sqrt( sigma**2 * s_xx / denom )
sigma_m = np.sqrt( sigma**2 * no / denom )

print('slope m = ',m,' with error= ',sigma_m)
print('intercept c = ',c,' with error= ',sigma_c)

IS.plot(S)
#plot final matrix

plt.show()
