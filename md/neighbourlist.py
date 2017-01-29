class neighbourlist(object):

    def __init__(self):
        raise NotImplementedError('You cannot create a neighbourlist object, just use the classmethods')
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

        from collections import defaultdict 
        neighbors = defaultdict(list)
        N, dim = np.shape(R)
        # assume same size in all 3 dimensions (i.e. L = L_x, L_y, L_z)
        # divide simulation box into small cells of equal size r_c >= r_cutoff
        n_boxes = np.int(box_length / r_cutoff)
        r_c = box_length / n_boxes 

        head = [-1] * n_boxes**3
        list = [-1] * N
        for i in range(0, N):
            # box index of particle by its position
            x = np.int(R[i][0] / r_c)
            y = np.int(R[i][1] / r_c)
            z = np.int(R[i][2] / r_c)
            ind_vec = np.int(box_length / r_c)
            box_index = x*ind_vec*ind_vec + y*ind_vec + z
            list[i] = head[box_index]
            # The last one goes to the head
            head[box_index] = i



        for x in range(ind_vec):
         for y in range(ind_vec):
          for z in range(ind_vec):
            box_index = x*ind_vec*ind_vec + y*ind_vec + z
            # Scan the neighbor boxes (including itself) of box c 
            nb = np.empty(dim)
            r_shift = np.empty(dim)
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
                        
                        
                        # Calculate box index of neighbor box # % pulls index back into appr. range
                        c1 = np.int(((nb[0]+ind_vec)%ind_vec) * ind_vec*ind_vec + ((nb[1]+ind_vec)%ind_vec) * ind_vec + ((nb[2]+ind_vec)%ind_vec))


                        # Scan particle i in cell box_index
                        i = head[box_index]
                        while(i != -1):
                            # Scan particle j in cell c1
                            j = head[c1]
                            while(j != -1):
                                # Avoid double counting of pair (i, j)
                                if (i<j):
                                    if (np.linalg.norm(R[i]-R[j]+r_shift) <= r_cutoff):
                                        neighbors[i].append(j)
                                        neighbors[j].append(i)
                                        

                                j = list[j]
                            i = list[i]
                                    
            

        return neighbors
