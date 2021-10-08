#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file stores all the funnctions needed.
@author: aarongorman
"""

h=0
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import random as rand

#function to initialise a random matrix with collinear spins randomly distributed
def init(N):   
    state = np.random.randint(2, size=(N,N))
    for i in range(N):
       for j in range(N):
           # converting to a matrix of and one and minus one.
           if state[i,j]==0:
               state[i,j]=-1
       #applying periodic boundary conditions
       state[0,j]=state[N-1,j]  #mapping top to bottom 
       state[j,0]=state[0,N-1]  #mapping left to right
    return state

#the main metropolis algorithm function coding for the monte carlo simulation 
def randomwalk(S, T, N):

    for i in range(N):
        for j in range(N):
                a = np.random.randint(N)
                b = np.random.randint(N)
                s =  S[a, b]
                #remainder a smarter way of applying periodic boundary conditions.
                Ham = ham(a,b,S,N)
                deltaE = -2*Ham
                if deltaE < 0:
                    s *= -1
                elif rand() < np.exp(-deltaE/T):
                    s *= -1
                S[a, b] = s
    return S

#define the Hamiltionian at the (i,j) site, of Spin matrix of size N^2
def ham(i,j,S,N):
    s = S[i,j]
    H = S[(i+1)%N, j] + S[i,(j+1)%N] + S[(i-1)%N, j] + S[i,(j-1)%N]
    hinner=-h*Mag(S)
    return -s*H+hinner

#calculating the energy of the spins in the whole microstate
def Eng(S,N):

    E = 0
    for i in range(len(S)):
        for j in range(len(S)):
            E+=ham(i,j,S,N)
    return E/4.0
# returning normalised energy

def Mag(S): #defining magnetisation as the sum of the spins in the system
    x = np.sum(S)
    return x

#monte carlo walk for a range of temperatures
#1 to reach equilibrium for a temperature, t
#2 to collect statistics on an ensemble of microstates at that temperature,t

def Metropolis(N,T):

    mag = 0
    magsq = 0
    eng = 0
    engsq = 0
    
    #init spin matrix w/ periodic boudary conditions
    S = init(N)

#burn in 
    for i in range(500):         # equilibrate steps
        S = randomwalk(S, T, N)           # Monte Carlo moves

#stats collecting steps
    for i in range(250):        #statistic collecting steps
        S = randomwalk(S, T, N)           
        E = Eng(S,N)      # calculate the energy
        M = Mag(S)       # calculate the magnetisation

        mag += M
        magsq += M**2.
        eng += E
        engsq += E**2.
    return mag, magsq, eng, engsq, S

#a function to show the spin map -all good
def plot(A):

    plt.figure()
    plt.imshow(A, cmap=plt.cm.gray) 
    plt.show()
    return