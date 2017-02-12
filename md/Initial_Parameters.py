#############################################################################################
########################################## PACKAGES #########################################
#############################################################################################

import numpy as np

#############################################################################################
##################################### INITIAL PARAMETERS ####################################
#############################################################################################
# Fill these before you start the Simulation

# Substance to be simulated
# e.g NaCl
Symbols = np.array(['Na','Cl'])

# Stochiometric Coefficients
# e.g. Na_1_Cl_1
Coefficients = np.array([1,1])

# Atomic Charges in e
# e.g. Na_1+_Cl_1-
Charges = np.array([1.0 ,-1.0])

# Number of Particles
N = 32

# Boxsize in Angström
L_x=22.56
L_y=22.56
L_z=22.56
L = np.array([L_x, L_y, L_z])

# LJ Cutoff Radius in Angström
r_cut_LJ = 0.4*L_x

# Accuracy Factor
# the cutoff-error is given by exp(-p)
p_error = 10.0

# Temperature in Kelvin
T = 100 # Kelvin 

# Timestep in femtoseconds
dt = 0.1

#Characetristic coupling time for Thermostat in femtoseconds
# must be larger than dt
tau = 100

#Switch Radius in Angström
r_switch = r_cut_LJ*0.9

#############################################################################################
#############################################################################################
#############################################################################################

#############################################################################################
############################### DO NOT CHANGE THESE LINES !!! ###############################
#############################################################################################

assert tau>dt, "tau must be larger than dt"
assert r_switch < r_cut_LJ, "switch radius must be smaller than LJ cutoff Radius"

# Summarizing Dimension in one array
L = np.array([L_x, L_y, L_z])

#Reassignment Probability
p_rea = dt/tau

#number of Boxes to consider for LJ-Potential
n_boxes_LJ = np.ceil(r_cut_LJ/np.max(L)).astype(int) 

##number of Boxes to consider for short ranged Potential
n_boxes_short_range = ( np.ceil((L_x / (float)(2))/np.max(L)) ).astype(int)