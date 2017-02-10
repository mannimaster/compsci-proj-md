/*  md - Molecular Dynamics Applied to ionic solids.
*   Copyright (C) 2017 Nils Harmening, Marco Manni,
*   Darian Steven Viezzer, Steffanie Kieninger, Henrik Narvaez
*
*   This program is free software: you can redistribute it and/or modify
*   it under the terms of the GNU General Public License as published by
*   the Free Software Foundation, either version 3 of the License, or
*   (at your option) any later version.
*
*   This program is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*   You should have received a copy of the GNU General Public License
*   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

double _fast_neighbourlist(double *R, double box_length, double r_cutoff)
{
  return 5;
}


/*

        N, dim = np.shape(R)
        # assume same size in all N dimensions
        n_cells = np.int(box_length / r_cutoff)
        # divide simulation box into small cells of equal size r_c >= r_cutoff
        r_c = box_length / n_cells 

        #define head and list 
        head = [-1] * n_cells**3
        cllist = [-1] * N
        for i in range(0, N):
            #empty list of neighbors
            neighbors[i] = []
            distances[i] = []
            # cell index of particle by its position
            x = np.int(R[i][0] / r_c)
            y = np.int(R[i][1] / r_c)
            z = np.int(R[i][2] / r_c)
            ind_vec = np.int(box_length / r_c)
            cell_index = x*ind_vec*ind_vec + y*ind_vec + z
            cllist[i] = head[cell_index]
            # The last one goes to the head
            head[cell_index] = i


        # For all cells: Look for neighbors within neighboring cells
        for cell in range(n_cells**3):
            x = cell/(ind_vec*ind_vec)
            y = (cell/ind_vec) % ind_vec
            z = cell % ind_vec

            nb = np.empty(dim)
            r_shift = np.empty(dim)
            # Scan the neighboring cells (including itself)
            nb[0] = x-1
            nb[1] = y-1
            nb[2] = z-1
            for nbcell_ind in range(1,3**3+1):
                if nbcell_ind != 1:
                    nb[2] += 1
                if nbcell_ind % 3 == 0:
                    nb[2] = z-1
                    nb[1] += 1
                if nbcell_ind % 9 == 0:
                    nb[1] = y-1
                    nb[0] += 1

                # Shift image position of simulation box?
                for d in range(dim):
                    if (nb[d] < 0):
                        r_shift[d] = -box_length
                    elif (nb[d]>=ind_vec):
                        r_shift[d] = box_length
                    else:
                        r_shift[d] = 0.0
                  
                # Calculate cell index nbcell of neighbor cell
                nbcell = np.int(((nb[0]+ind_vec)%ind_vec)* ind_vec*ind_vec
                          + ((nb[1]+ind_vec)%ind_vec) * ind_vec 
                          + ((nb[2]+ind_vec)%ind_vec))
                # where % pulls index back into appr. range

                # Scan particle i in cell 
                i = head[cell]
                while(i != -1):
                # Scan particle j in cell nbcell
                    j = head[nbcell]
                    while(j != -1):
                    # Avoid double counting of pair (i, j)
                        if (i<j):
                            # dist of i, j smaller than cutoff?
                            dist = np.linalg.norm(R[i]-(R[j]+r_shift))
                            if (dist <= r_cutoff):
                                neighbors[i].append(j)
                                distances[i].append(dist)
                                neighbors[j].append(i)
                                distances[j].append(dist)


                        j = cllist[j]
                    i = cllist[i]
                                    
        return neighbors, distances
*/
