#   md - Molecular Dynamics Applied to ionic solids.
#   Copyright (C) 2017 Nils Harmening, Marco Manni,
#   Darian Steven Viezzer, Steffi, Hendrik
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

"""
def test_md_callable():
    assert False
"""


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
