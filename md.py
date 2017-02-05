'''
Class for the Molecular Dynamics Simulation
'''
import numpy as np
import PSE
from neighbourlist import neighbourlist
from particle_interaction import coulomb
from particle_interaction import lennard_jones
from dynamics import dynamics
PSE = PSE.PSE

class System(object):
    """ This Class defines the Chemical Systems
    input must look as follows:
    
    e.g. NaCl
    
    Symbols = ['Na', Cl']  #Chemical Symbols of the types of Atmos
    
    Coeffcients = [1,1] #Stochiometric coefficients in the sum formula
    
    Charges = [1, -1] #Charges of the individual chemical species
    
    n = int #Number of Particles in the System
    
    Returns:
    ----------
    Nothing
    """
    def __init__(self,Symbols,Coefficients, Charges,n):
        self.Symbols = Symbols
        self.Coefficients = Coefficients
        self.Charges = Charges
        self.n = int(n)
        return
    
    def get_Labels(self):
        """Create the Nx3 Array Labels containing each Particles Mass, Charge and chemical Label """
        s = np.sum(self.Coefficients)
        m = np.zeros(self.n*s)
        q = np.zeros(self.n*s)
        chemical_labels  = np.zeros(self.n*s)
        size = np.size(self.Symbols)

        m[:self.n*self.Coefficients[0]] = PSE[ self.Symbols[0] ][1].astype('float64')
        for j in np.arange((np.size(self.Coefficients)-1)):
            m[self.n*np.cumsum(self.Coefficients)[j]:self.n*np.cumsum(self.Coefficients)[j+1]] = PSE[ self.Symbols[j+1] ][1].astype('float64')
        m *= 1.660539040e-27  # Correcting Unit, amu --> kg
        
        q[:self.n*self.Coefficients[0]] = self.Charges[0]
        for j in np.arange((np.size(self.Coefficients)-1)):
            q[self.n*np.cumsum(self.Coefficients)[j]:self.n*np.cumsum(self.Coefficients)[j+1]] = self.Charges[j+1]
        q *= 1.6021766208e-19 #Correcting Unit, 1 --> C


        index = np.zeros(np.size(self.Symbols))
        index = 2**np.arange(np.size(index))-1

        chemical_labels[:self.n*self.Coefficients[0]] = index[0]
        for j in np.arange((np.size(self.Coefficients)-1)):
            chemical_labels[self.n*np.cumsum(self.Coefficients)[j]:self.n*np.cumsum(self.Coefficients)[j+1]] = index[j+1]

        Labels  = np.zeros((self.n*s,3))
        Labels[:,0] = m
        Labels[:,1] = q
        Labels[:,2] =chemical_labels
        return Labels
    
    def get_LJ_parameter(self):
        
        index = np.zeros(np.size(self.Symbols))
        index = 2**np.arange(np.size(index))-1

        index = index.astype(int)
        Sigma = np.zeros((2*max(index)+1).astype(int) )
        Epsilon = np.zeros((2*max(index)+1).astype(int) )
        
        for i in index:
            l2_i = np.log2(i+1).astype(int)
            for j in index:
                l2_j = np.log2(j+1).astype(int)
                Sigma[i+j] = (PSE[self.Symbols[l2_i]][2].astype(float)+PSE[self.Symbols[l2_j]][2].astype(float)) /2.0

        for i in index:
            l2_i = np.log2(i+1).astype(int)
            for j in index:
                l2_j = np.log2(j+1).astype(int)
                Epsilon[i+j] = np.sqrt(PSE[self.Symbols[l2_i]][3].astype(float)*PSE[self.Symbols[l2_j]][3].astype(float))
        return Sigma, Epsilon
    
    
class md(object):

    '''
    Initializes the md object.

    Parameters
        ----------
        
        positions: Nx3 Array 
            Array with N rows and 3 columns. Contains the Positions of each Particle component wise
            
        R : Nx1 Array 
            Array with N rows and 1 column. Contains the euclidic distances of each particle to the coordinate origin. 
            
        properties : Nx3 Array
            The number of rows represents the number of prticles
            The columns represent the properties [q, m, label]
            q = particle charge
            m = particle mass
            label = number to identify the type of atom
            
        velocities : Nx3 Array
            Array with N rows and 3 columns. Contains the Velocites of each Particle component wise
            
        box : 3xArray
           Array containing the lengths of the edges of the box. 
            
        Temperature : float
            Temperature of the System in Kelvin.
            
        std : int 
            Standart deviation of the Gauss Distribution, used in the Ewald Summation
            
        Sigma_LJ : Array
            Contains the Lennard Jones Parameter sigma for each interaction pair
        
        Epsilon_LJ : Array
             Contains the Lennard Jones Parameter epsilon for each interaction pair
             
        switch_parameter : 4x1 Array
            Contains the values of the parameters used in the switch polynomial
        
        r_switch: float
            Radius where the switch funtion is applied
        
        n_boxes_short_range: int
            Number of Boxes to consider for the short ranged interactions
        
        k_max_long_range: int
            Upper limit for any entry of any vektor k of the reciprocal space
        
        dt: int
            Timestep of the Simulation
            
        k_cut: float
            Cutoff-radius in reciprocal space
        
        p_rea: float
             Reassingment probability. Denotes the coupling strength to the thermostat. It is the probability with which a particle will undergo a velocity reassignment. 

    Returns
    -------
        nothing
    '''
    def __init__(self, 
                 positions,
                 R,
                 properties, 
                 velocities,
                 forces,
                 box,
                 Temperature,
                 Sigma_LJ,
                 Epsilon_LJ, 
                 switch_parameter, 
                 r_switch,
                 r_cut_LJ,
                 std,
                 n_boxes_short_range,
                 dt,
                 p_rea,
                 p_error):
        #check input parameters
        self.positions=positions
        self.R=R
        self.labels=properties
        self.velocities = velocities
        self.forces = forces
        self.L=box
        self.T= Temperature
        self.std =  std
        self.lennard_jones = lennard_jones()
        self.Sigma_LJ = Sigma_LJ
        self.Epsilon_LJ = Epsilon_LJ
        self.switch_parameter = switch_parameter
        self.r_switch = r_switch
        self.r_cut_LJ = r_cut_LJ
        self.dt = dt
        self.p_rea = p_rea
        self.k_cut = 2 * p_error * 2 / (float)(box[0])
        self.n_boxes_short_range= n_boxes_short_range
        self.coulomb = coulomb(self.std, n_boxes_short_range, box, self.k_cut)
        self.r_cut_coulomb, self.k_cut = self.coulomb.compute_optimal_cutoff(positions,properties,box,p_error)
        self.neighbours_LJ, self.distances_LJ= neighbourlist().compute_neighbourlist(positions, box[0], self.r_cut_LJ)
        self.neighbours_coulomb, self.distances_coulomb= neighbourlist().compute_neighbourlist(positions, box[0], self.r_cut_coulomb)
        return
    
    @property
    def positions(self):
        return self._positions
    
    @positions.setter
    def positions(self,xyz):
        self._positions = xyz
     
    @property
    def R(self):
        return self._R
    
    @R.setter
    def R(self,new_R):
        self._R = new_R
        
    @property
    def velocities(self):
        return self._velocities
    
    @velocities.setter
    def velocities(self,xyz):
        self._velocities = xyz
    
    @property
    def forces(self):
        return self._forces    
    @forces.setter
    def forces(self,xyz):
        self._forces = xyz 
        
    @property
    def potential(self):
        return self._potential  
    
    @potential.setter
    def potential(self,xyz):
        self._potential = xyz
    
    
    def get_neighbourlist(self):
        """Compute the neighbourlist according to r_cut_coulomb and given configuration 
       
        
        Returns
        ..........
        
        neighbours: cell-linked list
            list at entry i contains all neighbours of particle i within cutoff-radius
        distances: cell-linked list
            list at entry i contains the distances to all neighbours of particle i within cutoff-radius
        """  
        neighbours, distances = neighbourlist().compute_neighbourlist(self.positions, self.L[0], self.r_cut_coulomb)
        return neighbours, distances
    
    
    
    
    # work in progress    
    def get_potential(self):
        """Compute the potential for the current configuration of the System
        
        Potential_total = Potential_Coulomb + Potential_Lennard_Jones
        
        Returns
        ..........
        
        Forces : float
            Number containg the Potential of the system with N particels 
        """        
        Potential = self.lennard_jones.compute_potential(sigma = self.Sigma_LJ,
                                                         epsilon = self.Epsilon_LJ, 
                                                         labels = self.labels,
                                                         neighbours = self.neighbours_LJ,
                                                         distances = self.distances_LJ)+(
        self.coulomb.compute_potential(labels = self.labels,
                                       positions = self.positions,
                                       neighbours = self.neighbours_coulomb,
                                       distances = self.distances_coulomb))
        
        return Potential
    
    
    def get_energy(self):
        
        Energy = self.lennard_jones.compute_energy(sigma = self.Sigma_LJ,
                                                   epsilon = self.Epsilon_LJ,
                                                   labels = self.labels,
                                                   neighbours = self.neighbours_LJ,
                                                   distances = self.distances_LJ)+(
        self.coulomb.compute_energy(labels = self.labels,
                                    positions = self.positions, 
                                    neighbours = self.neighbours_coulomb, 
                                    distances = self.distances_coulomb))
        return Energy
    
    
    
    
    
    
    def get_forces(self):
        """Compute the forces for the current configuration of the System
        
        F_total = F_Coulomb + F_Lennard_Jones
        
        Returns
        ..........
        
        Forces : N x 3 Array
            Array containg the Forces that act upon each particle component wise. 
        """
        Forces = self.lennard_jones.compute_forces(Positions =self.positions, 
                                                   Sigma =self.Sigma_LJ,
                                                   Epsilon = self.Epsilon_LJ, 
                                                   Labels =self.labels, 
                                                   L =self.L, 
                                                   switch_parameter = self.switch_parameter, 
                                                   r_switch = self.r_switch,
                                                   neighbours = self.neighbours_LJ)+(
        self.coulomb.compute_forces(Positions =self.positions,
                                      Labels = self.labels,
                                      L = self.L) )
        return Forces
    
    def propagte_system(self):
        """ Propagates the system by one timestep of length dt. Uses the Velocity Verlet Integrator and Andersen Thermostat.
        
        Returns 
        ----------
        Positions_new : N x 3 Array
            Array with N rows and 3 columns. Contains each particles positons component wise.
        
        Velocities_new : N x 3 Array
            Array with N rows and 3 columns. Contains each particles velocity component wise.
        
        Forces_new : N x 3 Array
            Array with N rows and 3 columns. Contains each the force acting upon each particle component wise.
        
        """
        Positions, Velocities, Forces = dynamics().compute_dynamics(self.positions,
                                                                    self.velocities, 
                                                                    self.forces, 
                                                                    self.labels,
                                                                    self.Sigma_LJ, 
                                                                    self.Epsilon_LJ,
                                                                    self.dt,
                                                                    self.L, 
                                                                    self.T,
                                                                    self.switch_parameter, 
                                                                    self.r_switch,
                                                                    self.neighbours_LJ,
                                                                    self.p_rea,
                                                                    self.coulomb,
                                                                    self.lennard_jones)
        
        return Positions, Velocities, Forces
    
    def get_Temperature(self):
        """ Calculate the instantaneous Temperature of the given Configuration.
        
        Returns
        -----------
        Temperature : float
            The instantaneous Temperature.
        """
       
        T = dynamics().Thermometer(self.labels, self.velocities)
        return T
        
                         
    
    
    
