import numpy as np

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

        #from collections import defaultdict 
        #neighbors = defaultdict(list)
        neighbors = {}
        distances = {}
        N, dim = np.shape(R)
        # assume same size in all N dimensions
        n_cells = np.ceil(box_length / r_cutoff).astype(int)
        # divide simulation box into small cells of equal size r_c >= r_cutoff
        r_c = box_length / n_cells 
        ind_vec = n_cells#np.int(box_length / r_c)

        #define head and list 
        head = [-1] * (n_cells+1)**3
        cllist = [-1] * N
        for i in range(0, N):
            #empty list of neighbors
            neighbors[i] = []
            distances[i] = []
            # cell index of particle by its position
            x = np.int(R[i][0] / r_c)
            y = np.int(R[i][1] / r_c)
            z = np.int(R[i][2] / r_c)
            cell_index = x*ind_vec*ind_vec + y*ind_vec + z
            cllist[i] = head[cell_index]
            # The last one goes to the head
            head[cell_index] = i


        # For all cells: Look for neighbors within neighboring cells
        for cell in range(n_cells**3):
            x = np.int(cell/(ind_vec*ind_vec)) % ind_vec
            y = np.int(cell/ind_vec) % ind_vec
            z = cell % ind_vec

            nb = np.empty(dim)
            r_shift = np.empty(dim)
            # Scan the neighboring cells (including itself)
            nb[0] = x-1
            nb[1] = y-1
            nb[2] = z-1
            for nbcell_ind in range(3**3):
                if nbcell_ind != 0:
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

      
      

    ####################################################
#THIS IS FOR TESTING AND DOES NOT WORK (WELL) 
#####################################################

    def naive_Dists(R, r_cut=0.8):
        '''Creates Neighbouring List using a naive approach

        Paramters
        ----------------

        R: Nx1 Array
            Array with N entries. Contains each Particles Distance to the coordinate origin.

        r_cut: float
            desired cutoff radius



        Returns
        -----------------
        Neighbors: List
        Neighbors[i] returns an array that contains the indices of all particles within the cutoff radius, with respect to particle i. 

        '''
        Neighbors = {}
        N = np.size(R)
        for i in np.arange(N):
            Neighbors[i]=[]
            for j in np.arange(N):
                if i != j:
                    dist = np.abs(R[i]-R[j])
                    if dist < r_cut:                     
                        Neighbors[i].append(j)
        return Neighbors


#    def compute_neighbourlist(self):
        # computeneighbourlist
#        return None

    def minDistance(self, position1, position2, box_length, half_box_length):
        """
        minDistance(position1,position2, box_length, half_box_length)

        minimal distancebetween two particles

        Parameters
        ----------
        position1, position2 : a skalar or np.array
            ..representing a position

        box_length : scalar number, positiv
            side length of the box otherwise denoted as "L"

        half_box_length : scalar number
            half the box_length

        Returns
        -------
        d : scalar number
            minimal distance between the inputed position considering periodic boundary conditions
        dr : vector distance
            minimal distance between the inputed position considering periodic boundary conditions
        """
        #d = np.linalg.norm((position1 - position2 + half_box_length) % box_length - half_box_length)
        dr = ((position1 - position2 + half_box_length) % box_length - half_box_length)
        d = np.linalg.norm(dr)
        return d, dr

    def neighborListInclDist(self, n_particles, positions, r_cutoff, r_s, box_length):
        """
        neighborList(n_particles, positions, r_cutoff, r_s, box_length)
        
        creates a neighbor list
        algorithm taken from https://github.com/markovmodel/compsci-2016/blob/master/lecture-notes/neighbor_lists.pdf p.10
        
        Parameters
        ----------
        n_particles : scalar number
                number of particles
                
        positions : np.array or single numbers
                array of the positionvectors with each line a position
                
        r_cutoff : scalar number
                cutoff radius
                
        r_s : scalar number
                skin radius of the maximal possible movement in one time step
                
        box_length : scalar number, positiv
                side length of the box otherwise denoted as "L"
                
        Returns
        -------
        neighbor_list : dictionary of lists
                all neighbors within given cutoff radius + skin radius
                
        distances : same as neighbor_list
                all distances as scalars
        """
        neighbor_list = {}
        distances = {}
        half_box_length = box_length/2
        
        for i in range(0,n_particles):
                if i not in neighbor_list:
                        #empty list of neighbors
                        neighbor_list[i] = []
                        distances[i] = []
                for j in range(i+1, n_particles):
                        d,dr = self.minDistance(positions[i],positions[j], box_length, half_box_length)
                        if  d < r_cutoff + r_s:
                                #neighborship i => j
                                neighbor_list[i].append(j)
                                distances[i].append(d)
                                #neighborship j => i
                                if j not in neighbor_list:
                                        neighbor_list[j] = [i]
                                        distances[j] = [d]
                                else:
                                        neighbor_list[j].append(i)
                                        distances[j].append(d)
        return neighbor_list, distances
