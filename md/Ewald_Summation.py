###################################
############ PACKAGES  ############
###################################

import numpy as np
import math
import matplotlib.pyplot as plt

###################################
############ FUNCTIONS  ###########
###################################

### EWALD LONG RANGE ###

# DUMMY 1D Long Range Ewald 
#THIS FUNCTION DO NOT WORK !!!!
def ewald_long_range_1D(trajParticles,Kmax,N,Lx,sigma,epsilon0):
  #calculates the structural factor S(k)
  #with k=2pi/Lx * n with |n|=0,1,2,...Kmax
  #s. lecture notes p. 39
  S=np.zeros(2*Kmax)
  for n in range(-Kmax,Kmax):
    for a in range(1,N):
      S[n+Kmax]=S[n+Kmax]+trajParticles{q[a]}*np.exp(-1j*(2*np.pi*n/Lx)*trajParticles{r[a]})
  
  #s. lecture notes p.40
  fi=np.zeros(N)
  for i in range(1,N):
    for n in range(-Kmax,Kmax):
      fi[i]=fi[i]+S[n+Kmax]*np.exp(-sigma**2*(2*np.pi*n/Lx)**2*0.5)*np.exp(1j*(2*np.pi*n/Lx)*trajParticles{r[i]})*(2*np.pi*n/Lx)*trajParticles{q[i]}/(epsilon0*Lx)
      
  return fi
