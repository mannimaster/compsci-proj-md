
# coding: utf-8

# In[2]:

import numpy as np
from boxvectors import directions as directions
import Initial_Parameters as ip
from md import System
from md import md
from distribution import maxwellboltzmann
import matplotlib.pyplot as plt
from scipy.special import erf
from scipy.special import erfc
from scipy.constants import epsilon_0
get_ipython().magic(u'matplotlib inline')


# In[3]:

Symbols = ip.Symbols
Coefficients = ip.Coefficients
Charges = ip.Charges
N = ip.N*np.sum(Coefficients)
L = ip.L
T = ip.T
dt = ip.dt
p_rea = ip.p_rea
std = ip.std
n_boxes_short_range = ip.n_boxes_short_range
p_error = ip.p_error
Sys= System(Symbols, Coefficients, Charges, N/2)
Labels = Sys.get_Labels()
Sigma, Epsilon = Sys.get_LJ_parameter()
r_cut_coulomb = ip.r_cut_coulomb
r_cut_LJ = ip.r_cut_LJ
r_switch = ip.r_switch
switch_parameter = ip.switch_parameter
m = Labels[:,0]


# In[4]:

def get_random_starting_Positions(N,L):
    Positions = np.zeros((N,3))
    Positions[:,0] = np.linspace(0.1/N,L[0],N, endpoint = False)
    Positions[:,1] = np.linspace(0.1/N,L[1],N, endpoint = False)
    Positions[:,2] = np.linspace(0.1/N,L[2],N, endpoint = False)
    np.random.shuffle(Positions[:,0])
    np.random.shuffle(Positions[:,1])
    np.random.shuffle(Positions[:,2])
    return Positions
Positions = get_random_starting_Positions(N,L)
Velocities = maxwellboltzmann().sample_distribution(N,m,T)
Forces = np.zeros((N,3))
R = np.linalg.norm(Positions,axis=1)


# In[5]:

MD = md(
    Positions, 
    R, 
    Labels, 
    Velocities,
    Forces, 
    L, 
    T, 
    Sigma, 
    Epsilon, 
    switch_parameter,
    r_switch,
    r_cut_LJ,
    dt, 
    std,
    n_boxes_short_range,
    p_rea,
    p_error)


# In[6]:

MD.forces = MD.get_forces()


# In[7]:

print MD.get_energy()
print MD.get_potential()
print MD.forces


# In[9]:

Temperature = np.zeros(10)
for i in np.arange(10):
    Positions_New, Velocities_New, Forces_New = MD.propagte_system()
    MD.positions = Positions_New
    MD.velocities = Velocities_New
    MD.forces = Forces_New
    MD.neighbours_LJ  = MD.get_neighbourlist_LJ()[0]
    Temperature[i] = MD.get_Temperature()
plt.plot(Temperature)


# In[12]:

MD.neighbours = 0
MD.neighbours


# In[ ]:



