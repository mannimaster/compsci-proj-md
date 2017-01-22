import numpy as np
class directions(object):
    
    def __init__(self, n_boxes):
        self.n_boxes = n_boxes
        return
    
    def get_directions(self):
        """ Returns
        -------------------
        K (2*n_boxes+1)**3 x 3 Array
            Each Row of the Array denotes a 1x3 Array that points to a unique dircetion"""
        
        assert type(self.n_boxes) is int, "n_boxes must be an integer"
        N_total_boxes = (2*self.n_boxes+1)**3
        decimals = np.arange(N_total_boxes)
        base = (2*self.n_boxes+1)

        K = np.zeros((N_total_boxes,3))

        K[:,0] = np.floor(decimals/base**2)
        K[:,1] = np.floor((decimals - K[:,0]*base**2) / base)
        K[:,2] = decimals - K[:,0]*base**2 - K[:,1]*base

        K -= self.n_boxes             # make sure that highest possible coefficient is nnumber of boxes in one direction
        return K
    
    def test_get_dircetions(self):

        #Create Test Array
        K = self.get_directions()
        #Shift all entries up, to avoid negative entries
        K += self.n_boxes
        #define Base
        base = (2*self.n_boxes+1)

        #Make a test Array
        K_test = np.zeros(base**3)
        K_test = K[:,0]*base**2 +K[:,1]*base +K[:,2]
        assert np.size(np.unique(K_test)) == base**3, "get_dircetions is broken"
        return "Passed"



