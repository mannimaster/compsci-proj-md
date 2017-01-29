class neighbourlist(object):

    def __init__(self):
        return

    def compute_neighbourlist(self, R, box_length, r_cutoff):
        # computeneighbourlist
        #return None

        """
        compute_neighbourlist(self, R, box_length, r_cutoff):
        returns list of neighbors within the cutoff radius for all particles
        Parameters
        ----------
        R : 2dim np.array (partictle, dim).
            Distance to origin of every particle in all 3 dim
        box_length : scalar number, positiv.
            Length of simulation box
        r_cutoff : scalar number, positiv.
            Cutoff radius, above that the interaction of two particles
            are neglegted.


        Returns
        -------
        neighbors : dictionary of list
            For each particle the list of its neighbors
            (i.e. distance < r_cutoff) are returned.
        """

        import numpy as np
        from collections import defaultdict 

        neighbors = defaultdict(list)
        N, dim = np.shape(R)
        # assume same size in all N dimensions
        n_cells = np.int(box_length / r_cutoff)
        # divide simulation box into small cells of equal size r_c >= r_cutoff
        r_c = box_length / n_cells 

        #define head and list 
        head = [-1] * n_cells**3
        cllist = [-1] * N
        for i in range(0, N):
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
        for x in range(ind_vec):
          for y in range(ind_vec):
            for z in range(ind_vec):
              cell = x*ind_vec*ind_vec + y*ind_vec + z

              nb = np.empty(dim)
              r_shift = np.empty(dim)
              # Scan the neighboring cells (including itself)
              for nb[0] in [x-1, x, x+1]:
                for nb[1] in [y-1, y, y+1]:
                  for nb[2] in [z-1, z, z+1]:
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
                          if (np.linalg.norm(R[i]-(R[j]+r_shift)) <= r_cutoff):
                            neighbors[i].append(j)
                            neighbors[j].append(i)
                                    

                        j = cllist[j]
                      i = cllist[i]
                                    
        return neighbors
