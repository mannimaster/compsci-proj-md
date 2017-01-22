'''
Class for the Molecular Dynamics Simulation
'''
import PSE
class System(object):
    """ This Class defines the Chemical Systems
    input must look as follows:
    e.g. NaCl
    Symbols = ['Na', Cl']  #Chemical Symbols of the types of Atmos
    Coeffcients = [1,1] #Stochiometric coefficients in the sum formula
    Charges = [1, -1] #Charges of the individual chemical species
    n = int #Number of Particles in the System
    """
    def __init__(self,Symbols,Coefficients, Charges,n):
        self.Symbols = Symbols
        self.Coefficients = Coefficients
        self.Charges = Charges
        self.n = n
        return
    
    def get_Labels(self):
        """ Creates an n x 3 Array Containing the Labels of the System.
        Each Row belongs to 1 Particle and contains [masses, charges, chemical_label]"""
        s = sum(Coefficients)
        m = np.zeros(n*s)
        q = np.zeros(n*s)
        chemical_labels  = np.zeros(n*s)
        size = np.size(Symbols)

        m[:n*Coefficients[0]] = PSE[ Symbols[0] ][1].astype('float64')
        for j in np.arange((np.size(Coefficients)-1)):
            m[n*np.cumsum(Coefficients)[j]:n*np.cumsum(Coefficients)[j+1]] = PSE[ Symbols[j+1] ][1].astype('float64')

        q[:n*Coefficients[0]] = Charges[0]
        for j in np.arange((np.size(Coefficients)-1)):
            q[n*np.cumsum(Coefficients)[j]:n*np.cumsum(Coefficients)[j+1]] = Charges[j+1]


        index = np.zeros(np.size(Symbols))
        for i in np.arange(np.size(index)):
            index[i] = 2**i-1

        chemical_labels[:n*Coefficients[0]] = index[0]
        for j in np.arange((np.size(Coefficients)-1)):
            chemical_labels[n*np.cumsum(Coefficients)[j]:n*np.cumsum(Coefficients)[j+1]] = index[j+1]

        Labels  = np.zeros((n*s,3))
        Labels[:,0] = m
        Labels[:,1] = q
        Labels[:,2] =chemical_labels
        return Labels
    
    def get_LJ_parameter(self):
        index = np.zeros(np.size(Symbols))
        for i in np.arange(np.size(index)):
            index[i] = 2**i-1

        index = index.astype(int)
        Sigma = np.zeros((2*max(index)+1).astype(int) )
        Epsilon = np.zeros((2*max(index)+1).astype(int) )
        
        for i in index:
            l2_i = np.log2(i+1).astype(int)
            for j in index:
                l2_j = np.log2(j+1).astype(int)
                Sigma[i+j] = (PSE[Symbols[l2_i]][2].astype(float)+PSE[Symbols[l2_j]][2].astype(float)) /2.0

        for i in index:
            l2_i = np.log2(i+1).astype(int)
            for j in index:
                l2_j = np.log2(j+1).astype(int)
                Epsilon[i+j] = np.sqrt(PSE[Symbols[l2_i]][3].astype(float)*PSE[Symbols[l2_j]][3].astype(float))
        return Sigma, Epsilon
    
class md(System):

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
            
        Temperature: int
            Temperature of the System in Kelvin.
                        

        The number of rows in particleposition and particleproperties has to be the same
        The number of coloumns in particleposition and box has to be the same
    Returns
    -------
        nothing
    '''
    def __init__(self,particleposition,particleproperties,box,Temperatre):
        #check input parameters
        self.position=particleposition
        self.labels=particleproperties
        self.L=box
        self.T= Temperature
        #return
        pass
