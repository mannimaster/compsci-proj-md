'''
Class for the Molecular Dynamics Simulation
'''
class md(object):

    '''
    Initializes the md object.

    Parameters
    ----------
        particleposition : array with n rows and m columns
                            The number of rows represents the number of particles
                            The number of columns represents the number of dimensions (typical up to 3)
                            Each row represents the position of a particle
        particleproperties : array with n rows and 4 columns
                            The number of rows represents the number of prticles
                            The columns represent the properties [q, m, LJsigma, LJepsilon]
                            q = particle charge
                            m = particle mass
                            LJsigma = Sigma-Parameter for the Lennard Jones (Zero cross distance)
                            LJepsilon = Epsilon-Parameter for the Lennard Jones (valley depth)
        box : array with 1 row and m coloumns
                The array contains the length of the box in each dimension
                A 1 dimensional Problem has to be given as a 1x1 array

        The number of rows in particleposition and particleproperties has to be the same
        The number of coloumns in particleposition and box has to be the same
    Returns
    -------
        nothing
    '''
    def __init__(self,particleposition,particleproperties,box):
        #check input parameters
        self.position=particleposition
        self.labels=particleproperties
        self.L=box
        #return
        pass