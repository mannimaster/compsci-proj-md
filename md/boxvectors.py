import numpy as np
class directions(object):
    
    def __init__(self, n_boxes):
        self.n_boxes = int(n_boxes)
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
    
