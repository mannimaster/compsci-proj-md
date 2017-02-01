
# coding: utf-8

# In[ ]:

import numpy as np
from boxvectors import directions as directions
import Initial_Parameters as ip
from md import System
from md import md
from distribution import maxwellboltzmann
import matplotlib.pyplot as plt
#get_ipython().magic(u'matplotlib inline')


# In[ ]:

Symbols = ip.Symbols
Coefficients = ip.Coefficients
Charges = ip.Charges
N = ip.N*np.sum(Coefficients)
L = ip.L
T = ip.T
dt = ip.dt
p_rea = ip.p_rea
p = ip.p
std = ip.std
k_cut = ip.k_cut
k_max = ip.k_max_long_range
n_boxes_LJ = ip.n_boxes_LJ
n_boxes_coulomb = ip.n_boxes_short_range
Sys= System(Symbols, Coefficients, Charges, N/2)
Labels = Sys.get_Labels()
Sigma, Epsilon = Sys.get_LJ_parameter()

m = Labels[:,0]


# In[ ]:

switch_parameter = np.array([1,-1,0,0])
r_switch = 1.5


# In[ ]:

def get_random_starting_Positions(N,L):
    Positions = np.zeros((N,3))
    Positions[:,0] = np.linspace(0.1,L[0],N, endpoint = False)
    Positions[:,1] = np.linspace(0.1,L[1],N, endpoint = False)
    Positions[:,2] = np.linspace(0.1,L[2],N, endpoint = False)
    np.random.shuffle(Positions[:,0])
    np.random.shuffle(Positions[:,1])
    np.random.shuffle(Positions[:,2])
    return Positions
Positions = get_random_starting_Positions(N,L)
Velocities = maxwellboltzmann().sample_distribution(N,m,T)
Forces = np.zeros((N,3))
R = np.linalg.norm(Positions,axis=1)


# In[ ]:

MD = md(
    Positions, 
    R, 
    Labels, 
    Velocities,
    Forces, 
    L, 
    T, 
    std, 
    Sigma, 
    Epsilon, 
    switch_parameter,
    r_switch, 
    n_boxes_coulomb, 
    k_max, 
    dt, 
    p_rea,
    k_cut)
MD.forces = MD.get_forces()

