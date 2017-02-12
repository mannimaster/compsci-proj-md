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
cimport numpy as np


cdef extern from 'src_fast_neighbourlist.h':
  double _fast_neighbourlist(double **R, int N, double box_length, double r_cutoff)


def fast_neighbourlist(
    np.ndarray[double, ndim=2, mode='c'] R not None,
    double box_length, double r_cutoff):

    cdef:
        np.ndarray[np.double_t, ndim=2, mode='c'] dist
        int i
        int j
        int k
    
    dist = _fast_neighbourlist(
        <double**> np.PyArray_DATA(R), R.shape[0], box_length, r_cutoff)

    neighbors = {}
    distances = {}
    for i in range(R.shape[0]):
        #empty list of neighbors
        neighbors[i] = []
        distances[i] = []
        j = 0
        k = 0
        for j in range(R.shape[0]):
            if i<j and dist[i][j]!=0.0:
                neighbors[i].append(j)
                neighbors[j].append(i)
                distances[i].append(dist[i][j])
                distances[j].append(dist[j][i])

    return neighbors, distances
