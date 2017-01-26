import numpy as np

class neighbourlist(object):

    def __init__(self):
        #raise NotImplementedError('You cannot create a neighbourlist object, just use the classmethods')
        return
    
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