
import numpy as np
# TESTFILE

#Number of Particles
N = 3

#Boxsize
L_x=3.
L_y=3.
L_z=3.
L = np.array([L_x, L_y, L_z])

#LJ Cutoff Radius
r_cut_LJ = 1.5

#Short-Range Potential Cutoff Radius
r_cut_coulomb = 1.5

#Accuracy Factor, the cutofferror is given by exp(-p)
p = 10.0

#Temperature 
T = 100 # Kelvin 

#Timestep
dt = 1e-3 # 1 ns

#Characetristic coupling time for Thermostat, must be larger than dt
tau = 1e-1
assert tau>dt, "tau must be larger than dt"

#labels : first column the masses, second the charge
labels = np.array(([1, 1, 0], [1, -2, 1], [1, 1, 0]))
neighbours = {0: [1, 2], 1: [0, 2], 2: [0, 1]}
distances = {0: [np.sqrt(3.), np.sqrt(3.)], 1: [np.sqrt(3.), np.sqrt(3.)], 2: [np.sqrt(3.), np.sqrt(3.)]}
positions = np.array(([1., 1., 1.], [2., 2., 2.], [3., 3., 3.]))
sigma     = np.array([1.,1.,1.])
epsilon   = np.array([1.,1.,1.])

#switch-parameter for Lennard-Jones-Forces
switch_parameter = np.array([1,-1,0,0])

#distance where the switch-function kicks in (Lennard-Jones-Forces)
r_switch = r_cut_LJ*0.8
###############################################
# !!!  DO NOT CHANGE THE FOLLOWING LINES  !!! #
###############################################

# Summarizing Dimension in one array
L = np.array([L_x, L_y, L_z])

#Reassignment Probability
p_rea = dt/tau

#Coulomb interaction sigma
std = r_cut_coulomb/np.sqrt(2*p)

#K_cut
k_cut = 2*p/r_cut_coulomb

#number of Boxes to consider for LJ-Potential
n_boxes_LJ = np.ceil(r_cut_LJ/np.max(L)).astype(int) 

##number of Boxes to consider for short ranged Potential
n_boxes_short_range = ( np.ceil(r_cut_coulomb/np.max(L)) ).astype(int)

#largest values of k to consider for long range Potential
k_max_long_range = int(np.floor((k_cut*L[0])/(2*np.pi)))


