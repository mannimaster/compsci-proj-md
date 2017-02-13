
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
r_cut_LJ = 1.0

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

# Calculate Switch Parameter by solving the following System of linear equations
#
# A*switch_parameters = [1, 0, 1, 1]
# <=> switch_parameters = A^-1 * [1, 0, 1, 1]
#
# s(x) = a + bx + cx^2 + d*x^3
# s'(x) = b + 2cx + 3dx^2
# s''(x) = 2c + 6dx
# (I):    s(r_switch) = 1
# (II):   s(r_cutoff) = 0
# (III): s'(r_switch) = 1
# (IV):  s''(r_switch = 1
A= np.array([
        [1, r_switch, r_switch**2, r_switch**3],
        [1, r_cut_LJ, r_cut_LJ**2, r_cut_LJ**3],
        [0, 1, 2*r_switch, 3*r_switch**2],
        [0, 0, 2, 6*r_switch]
    ])
switch_parameter = np.dot(np.linalg.inv(A),np.array([1,0,1,1]))
