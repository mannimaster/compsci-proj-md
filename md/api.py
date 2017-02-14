#   md - Molecular Dynamics Applied to ionic solids.
#   Copyright (C) 2017 Nils Harmening, Marco Manni,
#   Darian Steven Viezzer, Stefanie Kieninger, Henrik Narvaez
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
import numpy as np
import sys
from shutil import copyfile
import os.path



def get_random_starting_Positions(N, L):
    Positions = np.zeros((N, 3))
    Positions[:, 0] = np.linspace(0, L[0], N, endpoint=False)
    Positions[:, 1] = np.linspace(0, L[1], N, endpoint=False)
    Positions[:, 2] = np.linspace(0, L[2], N, endpoint=False)
    np.random.shuffle(Positions[:, 0])
    np.random.shuffle(Positions[:, 1])
    np.random.shuffle(Positions[:, 2])
    return Positions


def api(N_steps, threshold, Energy_save, Frame_save,Temperature_save):
    import Initial_Parameters as ip
    from md import System
    from md import md
    from distribution import maxwellboltzmann
    from scipy.special import erf
    from scipy.constants import epsilon_0
    import os
    import time
    cwd = os.getcwd()


    Symbols = ip.Symbols
    Coefficients = ip.Coefficients
    Charges = ip.Charges
    N = ip.N * np.sum(Coefficients)
    L = ip.L
    T = ip.T
    dt = ip.dt
    p_rea = ip.p_rea
    n_boxes_short_range = ip.n_boxes_short_range
    p_error = ip.p_error
    Sys = System(Symbols, Coefficients, Charges, N / 2)
    Labels = Sys.get_Labels()
    Sigma, Epsilon = Sys.get_LJ_parameter()
    r_cut_LJ = ip.r_cut_LJ
    r_switch = ip.r_switch
    m = Labels[:, 0]
    max_displacement = ip.max_displacement

    Positions = get_random_starting_Positions(N, L)
    Velocities = maxwellboltzmann().sample_distribution(N, m, T)
    Forces = np.zeros((N, 3))


    R = np.linalg.norm(Positions, axis=1)


    MD = md(
        Positions,
        Labels,
        Velocities,
        Forces,
        L,
        T,
        Sigma,
        Epsilon,
        r_switch,
        r_cut_LJ,
        n_boxes_short_range,
        dt,
        p_rea,
        p_error,
        Symbols)

    MD.forces = MD.get_forces()
    MD.minmimize_Energy(N_steps=N_steps, threshold=threshold, Energy_save=Energy_save, Frame_save=Frame_save, path=cwd,
                        max_displacement=max_displacement)

    MD.get_traj(N_steps, Energy_save, Temperature_save, Frame_save, path=cwd)



    return


flaggs = np.zeros(6)
if len(sys.argv)==1:
    import Initial_Parameters as ip

    # Starting api function
    api(N_steps=ip.N_steps, threshold=ip.threshold, Energy_save=ip.Energy_save, Frame_save=ip.Frame_save,
        Temperature_save=ip.Temperature_save)

elif len(sys.argv)>=2:

    #Help command
    if sys.argv[1] == "-h" or sys.argv[1] == "help":
        print "\npython api.py"
        print "Starts the api.py script with the parameters from the Initial_Parameters.py file that is located in the same directory as api.py. \n"
        print "\n###.Optional Arguments.###\n"
        print "-p <file_path>"
        print "<file_path> has to point on your modified Initial_Parameters.py (e.g. ./Initial_Parameters.py). Copies the target <file_path> to the directory of api.py and runs the script with these new parameters. \n"
        print "-N <integer>"
        print "Number of iterations the simulation should run.\n"
        print "-thr <float>"
        print "threshold: The stopping condition. (e.g. 1e-11)\n"
        print "-EnS <integer>"
        print "Every how many steps the Energy of the system sould be saved.\n"
        print "-FrS <integer>"
        print "Every how many steps the Positions of every particle should be saved.\n"
        print "TeS <integer>"
        print "Every how many steps the temperatur of every particle should be saved.\n"

    #checking for flaggs
    else:
        for i in range(len(sys.argv)):
            if sys.argv[i]=="-N":
                vN_steps     = int(sys.argv[i+1])
                flaggs[0] = 1
            if sys.argv[i]=="-thr":
                vthreshold   = float(sys.argv[i+1])
                flaggs[1] = 1
            if sys.argv[i]=="-EnS":
                vEnergy_save = int(sys.argv[i+1])
                flaggs[2] = 1
            if sys.argv[i]=="-FrS":
                vFrame_save  = int(sys.argv[i+1])
                flaggs[3] = 1
            if sys.argv[i]=="TeS":
                vTemperature_save= int(sys.argv[i+1])
                flaggs[4] = 1
            if sys.argv[i]=="-p":
                if os.path.isfile(sys.argv[i+1]):
                    flaggs[5] = 1
                    copyfile(sys.argv[i+1],'./Initial_Parameters.py')
                else:
                    print "Error: -p target file does not exist."

        import Initial_Parameters as ip
        for i in range(6):
            if flaggs[0] == 0:
                vN_steps = ip.N_steps
            if flaggs[1] == 0:
                vthreshold = ip.threshold
            if flaggs[2] == 0:
                vEnergy_save = ip.Energy_save
            if flaggs[3] == 0:
                vFrame_save = ip.Frame_save
            if flaggs[4] == 0:
                vTemperature_save = ip.Temperature_save


        #Starting api function
        api(N_steps=vN_steps,threshold=vthreshold,Energy_save=vEnergy_save,Frame_save=vFrame_save,Temperature_save=vTemperature_save)