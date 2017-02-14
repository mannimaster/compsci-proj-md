#   md - Molecular Dynamics Applied to ionic solids.
#   Copyright (C) 2017 Nils Harmening, Marco Manni,
#   Darian Steven Viezzer, Stefanie Kieninger, Henrik Narvaez
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

#from .api import md
from particle_interaction import coulomb
from particle_interaction import lennard_jones
from neighbourlist import neighbourlist 
import Initial_Test_Parameters as ip
from md import System

Symbols = ["Na" , "Cl"]
Coefficients = [1, 1]
Charges = [+1.0, -1.0]
N = 1*np.sum(Coefficients)
#k_cut = ip.k_cut
#k_max = ip.k_max_long_range
n_boxes_short_range = ip.n_boxes_short_range
n_boxes_LJ = ip.n_boxes_LJ
r_cut_LJ = ip.r_cut_LJ
r_switch = r_cut_LJ*0.9
Sys= System(Symbols, Coefficients, Charges, N/2)
Sigma, Epsilon = Sys.get_LJ_parameter()
p_error = ip.p
switch_parameter = ip.switch_parameter
Test_Positions = np.array([[1,1,1],
                           [2,2,2]])
Test_R = np.linalg.norm(Test_Positions)
Test_L = ip.L
Test_Labels = np.array([[1,+1.0,0],
                        [1,-1.0,1]])

Test_d_Pos = np.zeros((N,N,3))
Test_d_Pos[:,:,0] = np.subtract.outer(Test_Positions[:,0], Test_Positions[:,0])
Test_d_Pos[:,:,1] = np.subtract.outer(Test_Positions[:,1], Test_Positions[:,1])
Test_d_Pos[:,:,2] = np.subtract.outer(Test_Positions[:,2], Test_Positions[:,2])

neighbours = neighbourlist().compute_neighbourlist(Test_Positions, Test_L[0], r_cut_LJ)[0]
  
  
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

def test_coulomb_forces():
    c = coulomb(n_boxes_short_range,
                Test_L,
                p_error)
    c.compute_optimal_cutoff(p_error, Test_L, Test_Labels, neighbours[0], neighbours[1], r_switch, 0.49 * Test_L[0], Test_Positions)
    Force = c.compute_forces(Test_d_Pos,
                             Test_Labels,
                             Test_L)
    assert np.all(Force[0,:] == -Force[1,:]), "coulomb force is broken"
    
def test_LJ_forces():
    LJ = lennard_jones()
    Force = LJ.compute_forces(Test_Positions,
                              Sigma,
                              Epsilon,
                              Test_Labels,
                              Test_L,
                              switch_parameter,
                              r_switch,
                              neighbours)
    assert np.all(Force[0,:] == -Force[1,:]), "lennard Jones force is broken"

def test_neighborlist():
    N = 100
    R=np.random.rand(N,3)
    box_length=1.0
    r_cutoff=0.11
   
    naiveneighbors = {}
    dx = np.empty(3)
    for i in range(N):
        naiveneighbors[i] = []
        for j in range(N):
            d = 0.0
            for x in range(3):
                dx[x] = np.abs(R[i][x]-R[j][x])
                if (dx[x] < -box_length/2):
                    dx[x] += box_length
                elif (dx[x] > box_length/2):
                    dx[x] -= box_length

                d += dx[x]**2
            
            d = np.sqrt(d)
            if (d <= r_cutoff):
                if i>j:
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
    potential        = coulomb(ip.n_boxes_short_range,ip.L, ip.p)
    potential.compute_optimal_cutoff(p_error, Test_L, Test_Labels, neighbours[0], neighbours[1], r_switch, 0.49 * Test_L[0], Test_Positions)
    result           = potential.compute_potential(positions=ip.positions, labels=ip.labels, neighbours=ip.neighbours,
                                                   distances=ip.distances, r_s=np.sqrt(3*2**2), r_c=np.sqrt(3*2**2)+0.1)
    assert ( abs(result[0]/result[2])<1+10**(-8) ) , "Potential does not have the expected symmetrie. P1 and P3 should be the same."
    assert ( abs(result[0]/result[1])<0.5+10**(-8) ) , "Potential does not have the expected symmetrie. P2 should be P1*2."
    assert ( abs(result[2]/result[1])<0.5+10**(-8) ) , "Potential does not have the expected symmetrie. P2 should be P3*2."
    return


def test_SymmetriesPotLJ():
    #tests LJ Potential for  with equidistant identical charges

    potential = lennard_jones()
    result    = potential.compute_potential(sigma=ip.sigma, epsilon=ip.epsilon, labels=ip.labels,
                                            distances=ip.distances, neighbours=ip.neighbours, r_s=np.sqrt(3*2**2), r_c=np.sqrt(3*2**2)+0.1)
    assert ( abs(result[0]/result[2])<1+10**(-8) ) , "Potential does not have the expected symmetrie. P1 and P3 should be the same."
    assert ( abs(result[0]/result[1])<1+10**(-8) ) , "Potential does not have the expected symmetrie. P1 and P2 should be the same."
    assert ( abs(result[2]/result[1])<1+10**(-8) ) , "Potential does not have the expected symmetrie. P3 and P2 should be the same."
    return


def test_SymmetriesPotLJ2():

    potential = lennard_jones()
    result    = potential.compute_potential(sigma=ip.sigma, epsilon=ip.epsilon, labels=ip.labels,
                                            distances={0: [np.sqrt(12), np.sqrt(3)], 1: [np.sqrt(12), np.sqrt(12)], 2: [np.sqrt(3), np.sqrt(12)]},
                                            neighbours=ip.neighbours, r_s=np.sqrt(3*2**2), r_c=np.sqrt(3*2**2)+0.1)
    assert ( abs(result[0]/result[2])<1+10**(-8) ) , "Potential does not have the expected symmetrie. P1 and P3 should be the same."
    assert ( result[0]!=result[1] )                , "Potential does not have the expected symmetrie. P1 and P2 should not be the same."
    assert ( result[2]!=result[1] )                , "Potential does not have the expected symmetrie. P3 and P2 should not be the same."
    return

def test_create_k_list():
    Coulomb = coulomb(ip.n_boxes_short_range, ip.L, ip.p_error)
    list = Coulomb._coulomb__create_k_list(2, 5)
    result = np.where(np.linalg.norm(list, axis=1) >=2)[0]
    assert result.shape == np.array([]).shape, "Should not exist."

def test_switchfunction():
    Coulomb = coulomb(ip.n_boxes_short_range, ip.L, ip.p_error)
    assert Coulomb._coulomb__switchFunction(1.,1.,2.)==1.,"The switchfunction should be one here."
    assert Coulomb._coulomb__switchFunction(2., 1., 2.) == 0., "The switchfunction should be zero here."

    LJpotential = lennard_jones()
    assert LJpotential._lennard_jones__switchFunction(1.,1.,2.)==1.,"The switchfunction should be one here."
    assert LJpotential._lennard_jones__switchFunction(2., 1., 2.) == 0., "The switchfunction should be zero here."
