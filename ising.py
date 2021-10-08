#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file stores the code to simulate the basic thermodynamics of the ising model
@author: aarongorman
"""

import numpy as np
import matplotlib.pyplot as plt
import functions as IS

#define parameters
N=16
h=0

#initialise holding data structures
magnetisation  = []
susceptibility = []
energy =[]
specificH=[]

t=[] #list

#monte carlo walk for a range of temperatures
#1 to reach equilibrium for a temperature, t
#2 to collect statistics on an ensemble of microstates at that temperature,t

for m in range(30,90):

    T=m/30.0
    t.append(T)# over a temp range 1->3

    mag,magsq,eng,engsq,S = IS.Metropolis(N,T)
    
    energy.append(eng/(800.0*N*N))
    magnetisation.append(np.abs(mag)/(800.0*N*N))
    print(magnetisation)
    susceptibility.append((magsq/(800.0*N*N) - mag*mag/(800.0*800*N*N))/T)
    print(susceptibility)
    specificH.append((engsq/(800.0*N*N) - eng*eng/(800.0*800*N*N))/T**2)

#plotting variables derived
f = plt.figure(1)

plt.plot(t, energy, '.', color="yellow")
plt.xlabel("Temperature (($k_b$ K))")
plt.ylabel("Energy ")

f=plt.figure(2)
plt.xlabel("Temperature ($k_b$ K)")
plt.ylabel("Magnetization per spin ")
plt.plot(t, magnetisation, '.', color="orange")

f=plt.figure(3)
plt.plot(t, specificH, '.', color="blue")
plt.xlabel("Temperature (($k_b$ K))")
plt.ylabel("Specific Heat")

f=plt.figure(4)
plt.plot(t, susceptibility, '.', color="purple")
plt.xlabel("Temperature ($k_b$ K)")
plt.ylabel("Susceptibility")

#show map of spins
IS.plot(S)

plt.show()
