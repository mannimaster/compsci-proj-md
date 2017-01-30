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
from boxvectors import directions

    
def test_get_dircetions():

        #Create Test Array
        d = directions(3)
        K = d.get_directions()
        #Shift all entries up, to avoid negative entries
        K += self.n_boxes
        #define Base
        base = (2*self.n_boxes+1)

        #Make a test Array
        K_test = np.zeros(base**3)
        K_test = K[:,0]*base**2 +K[:,1]*base +K[:,2]
        assert np.size(np.unique(K_test)) == base**3, "get_dircetions is broken"
        return "Passed"

