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

# Boxsize in Angstroem
L_x=22.56
L_y=22.56
L_z=22.56


# LJ Cutoff Radius in Angstroem
r_cut_LJ = 0.4*L_x
assert r_cut_LJ < L_x/2, "LJ cutoff radius must be smaller than half a box length"

# Accuracy Factor
# the cutoff-error is given by exp(-p)
p_error = 10.0

# Temperature in Kelvin
T = 100 # Kelvin 

# Timestep in seconds
timestep = 0.1e-15

#Characetristic coupling time for Thermostat in femtoseconds
# must be larger than dt
tau = 100

#Switch Radius in Angstroem
r_switch = r_cut_LJ*0.9

#Number of itereations
N_steps=20

#stopping condition for the loop
threshold=1e-11

#Every how many iterations the Energy should be saved
Energy_save=1

#Every how many iterations the Positions should be saved
Frame_save=1

#Every how many itereations the Temperature of the system should be saved
Temperature_save=1


#############################################################################################
#############################################################################################
#############################################################################################

#############################################################################################
############################### DO NOT CHANGE THESE LINES !!! ###############################
#############################################################################################

assert tau>timestep, "tau must be larger than dt"
assert r_switch < r_cut_LJ, "switch radius must be smaller than LJ cutoff Radius"

# Summarizing Dimension in one array
L = np.array([L_x, L_y, L_z])

#Reassignment Probability
p_rea = timestep/tau

##number of Boxes to consider for short ranged Potential
n_boxes_short_range = ( np.ceil((L_x / (float)(2))/np.max(L)) ).astype(int)

dt = timestep/48.8882*1e15 #correcting Unit s --> dt = 48.8882 fs = 48.8882e-15 s


#test -p