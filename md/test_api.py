#   md - Molecular Dynamics Applied to ionic solids.
#   Copyright (C) 2017 Nils Harmening, Marco Manni,
#   Darian Steven Viezzer, Steffanie Kieninger, Henrik Narvaez
#
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

from .api import md
from particle_interaction import coulomb
import Initial_Parameters as ip
from md import System

Symbols = ip.Symbols
Coefficients = ip.Coefficients
Charges = ip.Charges
N = ip.N*np.sum(Coefficients)
std = ip.std
k_cut = ip.k_cut
k_max = ip.k_max_long_range
n_boxes_short_range = ip.n_boxes_short_range
n_boxes_LJ = ip.n_boxes_LJ
r_cut_LJ = ip.r_cut_LJ
r_switch = r_cut_LJ*0.9
Sys= System(Symbols, Coefficients, Charges, N/2)
Sigma, Epsilon = Sys.get_LJ_parameter()
switch_parameter = np.array([0,0,0,0])

Test_Positions = np.array([[1,0,0],
                           [3,0,0]])
Test_R = np.linalg.norm(Test_Positions)
Test_L = np.array([5,5,5])
Test_Labels = np.array([[1,+1.0,0],
                        [1,-1.0,1]])

    
def test_get_dircetions():
    from boxvectors import directions
    #Create Test Array
    d = directions(3)
    K = d.get_directions()
    #Shift all entries up, to avoid negative entries
    n_boxes = 3
    K += n_boxes
    #define Base
    base = (2*n_boxes+1)

    #Make a test Array
    K_test = np.zeros(base**3)
    K_test = K[:,0]*base**2 +K[:,1]*base +K[:,2]
    assert np.size(np.unique(K_test)) == base**3, "get_dircetions is broken"
    return "Passed"



def test_neighborlist():
    N = 100
    R=np.random.rand(N,3)
    box_length=1.0
    r_cutoff=0.1
   
    naiveneighbors = {}
    dx = np.empty(3)
    for i in range(N):
        naiveneighbors[i] = []
        for j in range(N):
            d = 0.0
            for x in range(3):
                dx[x] = R[i][x]-R[j][x]
                if (dx[x] < -box_length/2):
                    dx[x] += box_length
                elif (dx[x] > box_length/2):
                    dx[x] -= box_length

                d += dx[x]**2
            
            d = np.sqrt(d)
            if (d <= r_cutoff):
                if (i>j):
                    naiveneighbors[i].append(j)
                    naiveneighbors[j].append(i)

   
    from neighbourlist import neighbourlist as nbl
    n2 = naiveneighbors
    n_inst = nbl()
    n1, dist = n_inst.compute_neighbourlist(R, box_length, r_cutoff)
    for i in range(N):
      n1[i].sort()
      n2[i].sort()

    assert n1 == n2



def test_SymmetriesPotC():
    #tests coulomb potential function with equidistant charges where the middle one has twice the negativ charge
    from particle_interaction import coulomb

    potential        = coulomb(ip.std, ip.n_boxes_short_range,ip.L, ip.k_max_long_range, ip.k_cut)
    result           = potential.compute_potential(positions=ip.positions, labels=ip.labels, neighbours=ip.neighbours, distances=ip.distances)
    assert ( abs(result[0]/result[2])<1+10**(-8) ) , "Potential does not have the symmetrie that the should follow from the particle position and charge. P1 and P3 should be the same."
    assert ( abs(result[0]/result[1])<0.5+10**(-8) ) , "Potential does not have the symmetrie that the should follow from the particle position and charge. P2 should be P1*2."
    assert ( abs(result[2]/result[1])<0.5+10**(-8) ) , "Potential does not have the symmetrie that the should follow from the particle position and charge. P2 should be P3*2."
    return


def test_SymmetriesPotLJ():
    #tests LJ Potential for  with equidistant identical charges
    from particle_interaction import lennard_jones

    potential = lennard_jones()
    result    = potential.compute_potential(sigma=ip.sigma, epsilon=ip.epsilon, labels=ip.labels, distances=ip.distances, neighbours=ip.neighbours)
    assert ( abs(result[0]/result[2])<1+10**(-8) ) , "Potential does not have the symmetrie that the should follow from the particle position and charge. P1 and P3 should be the same."
    assert ( abs(result[0]/result[1])<1+10**(-8) ) , "Potential does not have the symmetrie that the should follow from the particle position and charge. P1 and P2 should be the same."
    assert ( abs(result[2]/result[1])<1+10**(-8) ) , "Potential does not have the symmetrie that the should follow from the particle position and charge. P3 and P2 should be the same."
    return


def test_SymmetriesPotLJ2():
    from particle_interaction import lennard_jones

    potential = lennard_jones()
    result    = potential.compute_potential(sigma=ip.sigma, epsilon=ip.epsilon, labels=ip.labels, distances={0: [np.sqrt(12), np.sqrt(3)], 1: [np.sqrt(12), np.sqrt(12)], 2: [np.sqrt(3), np.sqrt(12)]}, neighbours=ip.neighbours)
    assert ( abs(result[0]/result[2])<1+10**(-8) ) , "Potential does not have the symmetrie that the should follow from the particle position and charge. P1 and P3 should be the same."
    assert ( result[0]!=result[1] )                , "Potential does not have the symmetrie that the should follow from the particle position and charge. P1 and P2 should not be the same."
    assert ( result[2]!=result[1] )                , "Potential does not have the symmetrie that the should follow from the particle position and charge. P3 and P2 should not be the same."
    return

    
def test_coulomb_forces():
    c = coulomb(std,
                n_boxes_short_range,
                Test_L,
                k_max,
                k_cut)
    Force = c.compute_forces(Test_Positions,
                             Test_Labels,
                             Test_L)
    assert np.all(Force[0,:] == -Force[1,:]), "coulomb force is broken"
    
def test_LJ_range_forces():
    LJ = lennard_jones()
    Force = LJ.compute_forces(Test_Positions,
                              Sigma,
                              Epsilon,
                              Test_Labels,
                              Test_L,
                              switch_parameter,
                              r_switch,
                              r_cut_LJ)
    assert np.all(Force[0,:] == -Force[1,:]), "lennard Jones force is broken"

