'''
Class for the Molecular Dynamics Simulation
'''
import numpy as np
import PSE
from neighbourlist import neighbourlist
from particle_interaction import coulomb
from particle_interaction import lennard_jones
from dynamics import dynamics
from scipy.constants import epsilon_0
import sys
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
        m *= 1# Correcting Unit, g/Mol --> g/Mol


        q[:self.n*self.Coefficients[0]] = self.Charges[0]
        for j in np.arange((np.size(self.Coefficients)-1)):
            q[self.n*np.cumsum(self.Coefficients)[j]:self.n*np.cumsum(self.Coefficients)[j+1]] = self.Charges[j+1]
        q *= 1 #Correcting Unit, e --> e  [1.6021766208e-19C]


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
                 properties, 
                 velocities,
                 forces,
                 box,
                 Temperature,
                 Sigma_LJ,
                 Epsilon_LJ, 
                 r_switch,
                 r_cut_LJ,
                 n_boxes_short_range,
                 dt,
                 p_rea,
                 p_error,
                 Symbols):
        #check input parameters
        self.positions=positions
        self.R=np.linalg.norm(self.positions, axis=1)
        self.N = np.size(self.positions[:,0])
        first_d_Pos = np.zeros((self.N,self.N,3))
        first_d_Pos[:,:,0] = np.subtract.outer(self.positions[:,0], self.positions[:,0])
        first_d_Pos[:,:,1] = np.subtract.outer(self.positions[:,1], self.positions[:,1])
        first_d_Pos[:,:,2] = np.subtract.outer(self.positions[:,2], self.positions[:,2])
        self.d_Pos = first_d_Pos
        self.labels=properties
        self.velocities = velocities
        self.forces = forces
        self.L=box
        self.T= Temperature
        
        self.lennard_jones = lennard_jones()
        self.Sigma_LJ = Sigma_LJ
        self.Epsilon_LJ = Epsilon_LJ
        
        
        self.dt = dt
        self.p_rea = p_rea
        self.n_boxes_short_range= n_boxes_short_range

        self.r_switch = r_switch
        self.r_cut_LJ = r_cut_LJ

        # epsilon0 = (8.854 * 10^-12) / (36.938 * 10^-9) -> see Dimension Analysis
        self.coulomb = coulomb(n_boxes_short_range, box, p_error, epsilon0 = epsilon_0 / (36.938 * 10**-9))

        self.neighbours_coulomb, self.distances_coulomb = neighbourlist().compute_neighbourlist(positions, box[0],
                                                                                                0.49 * box[0])
        self.r_cut_coulomb, self.k_cut, self.std = self.coulomb.compute_optimal_cutoff(p_error, box, properties, self.neighbours_coulomb, self.distances_coulomb, r_switch, 0.49 * box[0], positions)
        self.neighbours_coulomb, self.distances_coulomb = neighbourlist().compute_neighbourlist(positions, box[0],
                                                                                                self.r_cut_coulomb)

        self.coulomb.n_boxes_short_range = np.ceil( self.r_cut_coulomb/self.L[0] ).astype(int)
        self.switch_parameter = self.__get_switch_parameter()
        self.neighbours_LJ, self.distances_LJ= neighbourlist().compute_neighbourlist(positions, box[0], self.r_cut_LJ)
        self.Symbols = Symbols

        return
    
    @property
    def positions(self):
        return self._positions
    
    @positions.setter
    def positions(self,xyz):
        self._positions = xyz
        self.R = np.linalg.norm(xyz, axis=1)
        new_d_Pos = np.zeros( (xyz.shape[0], xyz.shape[0], 3) )
        new_d_Pos[:,:,0] = np.subtract.outer(xyz[:,0], xyz[:,0])
        new_d_Pos[:,:,1] = np.subtract.outer(xyz[:,1], xyz[:,1])
        new_d_Pos[:,:,2] = np.subtract.outer(xyz[:,2], xyz[:,2])
        self.d_Pos = new_d_Pos
        return
     
    @property
    def R(self):
        return self._R
    
    @R.setter
    def R(self,new_R):
        self._R = new_R
        return
    
    @property
    def velocities(self):
        return self._velocities
    
    @velocities.setter
    def velocities(self,xyz):
        self._velocities = xyz
        return
    
    @property
    def forces(self):
        return self._forces    
    
    @forces.setter
    def forces(self,xyz):
        self._forces = xyz 
        return
    
    @property
    def potential(self):
        return self._potential  
    
    @potential.setter
    def potential(self,xyz):
        self._potential = xyz
        return
    
    def __get_switch_parameter(self):
        A = np.array([ 
            [1, self.r_switch, self.r_switch**2, self.r_switch**3, self.r_switch**4, self.r_switch**5], 
            [1, self.r_cut_LJ, self.r_cut_LJ**2, self.r_cut_LJ**3, self.r_cut_LJ**4, self.r_cut_LJ**5],
            [0, 1, 2*self.r_switch, 3*self.r_switch**2, 4*self.r_switch**3, 5*self.r_switch**4], 
            [0, 1, 2*self.r_cut_LJ, 3*self.r_cut_LJ**2, 4*self.r_cut_LJ**3, 5*self.r_cut_LJ**4],
            [0, 0, 2, 6*self.r_switch, 12*self.r_switch**2, 20*self.r_switch**3],
            [0, 0, 2, 6*self.r_cut_LJ, 12*self.r_cut_LJ**2, 20*self.r_cut_LJ**3]])
        switch_parameter = np.dot(np.linalg.inv(A), np.array([1,0,0,0,0,0]))
        return switch_parameter

    
    
    
    def get_neighbourlist_coulomb(self):
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
    
    def get_neighbourlist_LJ(self):
        """Compute the neighbourlist according to r_cut_coulomb and given configuration 
       
        
        Returns
        ..........
        
        neighbours: cell-linked list
            list at entry i contains all neighbours of particle i within cutoff-radius
        distances: cell-linked list
            list at entry i contains the distances to all neighbours of particle i within cutoff-radius
        """  
        neighbours, distances = neighbourlist().compute_neighbourlist(self.positions, self.L[0], self.r_cut_LJ)
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
                                                         distances = self.distances_LJ,
                                                         r_s=self.r_switch,
                                                         r_c=self.r_cut_coulomb)+(
        self.coulomb.compute_potential(labels = self.labels,
                                       positions = self.positions,
                                       neighbours = self.neighbours_coulomb,
                                       distances = self.distances_coulomb,
                                       r_s=self.r_switch,
                                       r_c=self.r_cut_coulomb,))
        
        return Potential
    
    
    def get_energy(self):
        
        Energy = self.lennard_jones.compute_energy(sigma = self.Sigma_LJ,
                                                   epsilon = self.Epsilon_LJ,
                                                   labels = self.labels,
                                                   neighbours = self.neighbours_LJ,
                                                   distances = self.distances_LJ,
                                                   r_c=self.r_cut_LJ,
                                                   r_s=self.r_switch)+(
        self.coulomb.compute_energy(labels = self.labels,
                                    positions = self.positions, 
                                    neighbours = self.neighbours_coulomb, 
                                    distances = self.distances_coulomb,
                                    r_s=self.r_switch,
                                    r_c=self.r_cut_coulomb))
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
        self.coulomb.compute_forces(d_Pos = self.d_Pos,
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
        Positions, Velocities, Forces = dynamics().velocity_verlet_integrator(self.positions,
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
                                                                               self.lennard_jones,
                                                                               self.d_Pos,
                                                                               thermostat = True)
        
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
        
                         
    
    def get_traj(self, N_steps, Energy_save, Temperature_save, Frame_save, path):
        """Propagates the System unitil convergence is reached or the maximum Number of Steps is reached

        Parameters
        ----------------
        N_steps: int
            The number of steps the simulation will maximally run for.

        Energy_save: int
            The intervall at which energies are measured and saved.

        Temperature_save: int
            The intervall at which temperatures are measured and saved.

        Frame_save: int
            The intervall at which frames are created and saved. 

        Path: string
            Location where results will be saved.

        Returns
        -----------------
        "Simulation Completed"
        Results will be saved in the specified path. Files will be named as follows:
        Trajectory : traj.xyz
        Energies: Energies
        Temperature: Temperature

        """
        traj_file = ''.join([path,"\\traj.xyz"])
        Energy_file = ''.join([path,"\\Energies"])
        Temperature_file = ''.join([path,"\\Temperature"])
        rdf_file = ''.join([path,"\\rdf"])
        string1 = (''.join([str(self.N), "\n", "\n"]))

        #write header
        myfile = open(traj_file,'w')
        myfile.write(str(self.N)+"\n"+"\n")
        myfile.close()

        frame = np.zeros((self.N), dtype=[('var1', 'U16'), ('var2',float), ('var3', float), ('var4',float)])
        frame['var1'] = self.Symbols[self.labels[:,2].astype(int)]

        #Positions in Angstroem
        frame['var2'] = self.positions[:, 0]
        frame['var3'] = self.positions[:, 1]
        frame['var4'] = self.positions[:, 2]

        myfile = open(traj_file,'ab')
        np.savetxt(myfile,frame, fmt = "%s %f8 %f8 %f8", )
        myfile.close()
        myfile = open(traj_file,'a')
        myfile.write(string1)
        myfile.close()

        Energy = np.zeros(np.ceil(N_steps/Energy_save).astype(int))
        Temperature = np.zeros(np.ceil(N_steps/Temperature_save).astype(int))

        counter_Energy = 0
        counter_Temperature = 0
        counter_Frame = 0
        E_index = -1

        #######################################################################
        ### write radial distribution function of current frame into a file ###
        #######################################################################
        number_part_A = np.unique(self.labels[:, 2], return_counts=True)[1][0]
        number_part_B = np.unique(self.labels[:, 2], return_counts=True)[1][1]

        totdismatrix = np.linalg.norm(self.d_Pos, axis=2)
        distancematrix_particle_A = totdismatrix[number_part_B:self.N, 0:number_part_A]
        distancematrix_particle_B = totdismatrix[0:number_part_A, number_part_B:self.N]

        # Volumenelement hat Dicke 1 Angstroem
        # maximal bis L/2 wird die Distribution genommen
        histvector_particle_A = self.radial_distribution(distancematrix_particle_A, 1, self.L[0]/2, self.N/2, self.N / self.L[0]**3)
        histvector_particle_B = self.radial_distribution(distancematrix_particle_B, 1, self.L[0]/2, self.N/2, self.N / self.L[0]**3)

        myfile = open(rdf_file, 'w')
        myfile.write(str(histvector_particle_A) + "\n" + str(histvector_particle_B) + "\n")
        myfile.close()
        #######################################################################
        #######################################################################
        #######################################################################

        for i in np.arange(N_steps):

            Positions_New, Velocities_New, Forces_New = self.propagte_system()
            self.positions = Positions_New
            self.velocities = Velocities_New
            self.forces = Forces_New
            self.neighbours_LJ  = self.get_neighbourlist_LJ()[0]

            counter_Energy += 1
            counter_Temperature += 1
            counter_Frame += 1

            #Save Energy
            if counter_Energy > Energy_save-1:
                E_index += 1
                Energy[E_index] = self.get_energy()
                counter_Energy = 0

            #Save Temperature
            if counter_Temperature > Temperature_save-1:
                Temperature[np.floor(i/Temperature_save).astype(int)] = self.get_Temperature()
                counter_Temperature = 0

            # Save Frame
            if counter_Frame > Frame_save-1 :
                #create frame
                frame = np.zeros((self.N), dtype=[('var1', 'U16'), ('var2',float), ('var3', float), ('var4',float)])
                frame['var1'] = self.Symbols[self.labels[:,2].astype(int)]

                frame['var2'] = self.positions[:,0]
                frame['var3'] = self.positions[:,1]
                frame['var4'] = self.positions[:,2]

                #save frame
                myfile = open(traj_file,'ab')
                np.savetxt(myfile,frame, fmt = "%s %f8 %f8 %f8", )
                myfile.close()
                myfile = open(traj_file,'a')
                myfile.write(string1)
                myfile.close()

                #######################################################################
                ### write radial distribution function of current frame into a file ###
                #######################################################################
                number_part_A = np.unique(self.labels[:, 2], return_counts=True)[1][0]
                number_part_B = np.unique(self.labels[:, 2], return_counts=True)[1][1]

                totdismatrix = np.linalg.norm(self.d_Pos, axis=2)
                distancematrix_particle_A = totdismatrix[number_part_B:self.N, 0:number_part_A]
                distancematrix_particle_B = totdismatrix[0:number_part_A, number_part_B:self.N]

                histvector_particle_A = self.radial_distribution(distancematrix_particle_A, 1, self.L[0] / 2,
                                                                 self.N / 2, self.N / self.L[0] ** 3)
                histvector_particle_B = self.radial_distribution(distancematrix_particle_B, 1, self.L[0] / 2,
                                                                 self.N / 2, self.N / self.L[0] ** 3)

                myfile = open(rdf_file, 'a')
                myfile.write(str(histvector_particle_A) + "\n" + str(histvector_particle_B) + "\n")
                myfile.close()
                #######################################################################
                #######################################################################
                #######################################################################

                counter_Frame = 0

            sys.stdout.write("\r")
            sys.stdout.write( ''.join([str(float(i+1)/N_steps*100), "% of steps completed"]))
            sys.stdout.flush()

        # save Energy       
        np.savetxt(Energy_file, Energy)
        # save Temperature
        np.savetxt(Temperature_file, Temperature)

        print("Simulation Completed")
        return

    def import_traj(self, file_path):
        '''
        converts a xyz-file back to a trajectory

        Parameters
        ----------------
        file_path: string
            full path to the xyz-file containing the trajectories

        Returns
        ----------------
        particle_type_A: string
            type of particle A

        particle_num_A: int
            number of particle A in the trajectory
            the first particle_num_A entries of the position_array in each frame belongs to particle A

        particle_type_B: string
            type of particle B

        particle_num_B:
            number of particle B in the trajectory
            the last particle_num_B entries in each frame belongs to particle B

        position_array: M x N x 3 array
            array containing the frames with the positions of all particles
            M -> frame
            N -> particle
            3 -> 3 dimensional (x, y, z)

        '''

        #read in the lines and convert the list into a numpy array of strings
        myfile = open(file_path, 'r')
        lines = myfile.readlines()
        position_array = np.array(lines).astype(str)

        #reads the first line to determine the total number of particles
        particle_num = int(np.fromfile(file_path, int, 1, "\n"))
        #reads the first table entry to determine the particle type A
        particle_type_A = position_array[2][0:position_array[2].find(' ')]
        #reads the last table entry to determine the particle type A
        particle_type_B = position_array[particle_num + 1][0:position_array[particle_num + 1].find(' ')]

        #deletes all blank lines and all lines containing the number of particles
        position_array = np.delete(position_array, np.where(position_array == str(particle_num) + '\n'))
        position_array = np.delete(position_array, np.where(position_array == '\n'))

        #counts the number of particle A and B
        particle_num_A = np.array([l.split(' ') for l in position_array], dtype=None)[0:particle_num, 0].tolist().count(
            'A')
        particle_num_B = np.array([l.split(' ') for l in position_array], dtype=None)[0:particle_num, 0].tolist().count(
            'B')

        #deletes the particle type in all rows and the \n
        position_array = np.char.replace(position_array, particle_type_A + ' ', '')
        position_array = np.char.replace(position_array, particle_type_B + ' ', '')
        position_array = np.char.replace(position_array, '\n', '')
        #converts the blankspace seperated numbers in the string to float values
        position_array = np.array([l.split(' ') for l in position_array], dtype=np.float32)

        #calculates the number of trajectories in the array
        number_of_traj = int(position_array.shape[0] / particle_num)

        #reshapes the 1D array to the M x N x 3 array
        position_array = position_array.reshape(number_of_traj, particle_num, 3)

        return (particle_type_A, particle_num_A, particle_type_B, particle_num_B, position_array)


    def minmimize_Energy(self,N_steps, threshold, Energy_save, Frame_save, max_displacement, path):
        """Minimizes the Energy by steepest descent

        Parameters
        ----------------
        N_steps: int
            The number of steps the simulation will maximally run for.
            
        threshold: float
            if the line sum norm of the forces falls below the threshold, the minimization is completed

        Energy_save: int
            The intervall at which energies are measured and saved.

        Frame_save: int
            The intervall at which frames are created and saved. 
            
        max_displacement: float
            Fraction of the boxlenth that particles will move with one Simulation step.

        Path: string
            Location where results will be saved.

        Returns
        -----------------
        Means by which the minimization was completed.
        either "Energy converged" 
        or "Maximum Number of steps reached
        
        Results will be saved in the specified path. Files will be named as follows:
        Trajectory : traj.xyz
        Energies: Energies"""
        
        traj_file = ''.join([path,"\\traj_minimization.xyz"])
        Energy_file = ''.join([path,"\\Energies_minimization"])
        string1 = (''.join([str(self.N), "\n", "\n"]))

        #write header
        myfile = open(traj_file,'w')
        myfile.write(str(self.N)+"\n"+"\n")
        myfile.close()

        frame = np.zeros((self.N), dtype=[('var1', 'U16'), ('var2',float), ('var3', float), ('var4',float)])
        frame['var1'] = self.Symbols[self.labels[:,2].astype(int)]

        #Positions in Angstroem
        frame['var2'] = self.positions[:,0]
        frame['var3'] = self.positions[:,1]
        frame['var4'] = self.positions[:,2]

        myfile = open(traj_file,'ab')
        np.savetxt(myfile,frame, fmt = "%s %f8 %f8 %f8", )
        myfile.close()
        myfile = open(traj_file,'a')
        myfile.write(string1)
        myfile.close()

        Energy = np.zeros(np.ceil(N_steps/Energy_save).astype(int))

        counter_Energy = 0
        counter_Frame = 0
        E_index = -1
        
        for i in np.arange(N_steps):
            
            
            constant = max_displacement/np.mean(np.linalg.norm(self.forces, axis=1))*self.L[0]
        
            #Update Positions
            Positions_new = dynamics().steepest_descent(self.positions,self.forces,self.L, constant)

            #Update Self
            self.positions = Positions_new
            self.neighbours_LJ = self.get_neighbourlist_LJ()[0]
            self.forces = self.get_forces()
            
            counter_Energy += 1
            counter_Frame += 1
            #Save Energy
            if counter_Energy > Energy_save-1:
                E_index += 1
                Energy[E_index] = self.get_energy()
                counter_Energy = 0

            # Save Frame
            if counter_Frame > Frame_save-1 :
                #create frame
                frame = np.zeros((self.N), dtype=[('var1', 'U16'), ('var2',float), ('var3', float), ('var4',float)])
                frame['var1'] = self.Symbols[self.labels[:,2].astype(int)]

                #Positions in Angstroem
                #frame['var2'] = self.positions[:,0]*1e10
                #frame['var3'] = self.positions[:,1]*1e10
                #frame['var4'] = self.positions[:,2]*1e10
                frame['var2'] = self.positions[:, 0]
                frame['var3'] = self.positions[:, 1]
                frame['var4'] = self.positions[:, 2]

                #save frame
                myfile = open(traj_file,'ab')
                np.savetxt(myfile,frame, fmt = "%s %f8 %f8 %f8", )
                myfile.close()
                myfile = open(traj_file,'a')
                myfile.write(string1)
                myfile.close()

                counter_Frame = 0
            
            #show Progress
            sys.stdout.write("\r")
            sys.stdout.write( ''.join([str(float(i+1)/N_steps*100), "% of steps completed"]))
            sys.stdout.flush()
            
            # Calculate norm
            norm_Forces = np.max(np.sum(np.abs(self.forces),1))
            
            if norm_Forces < threshold:
                #create frame
                frame = np.zeros((self.N), dtype=[('var1', 'U16'), ('var2',float), ('var3', float), ('var4',float)])
                frame['var1'] = self.Symbols[self.labels[:,2].astype(int)]

                #Positions in Angstroem
                #frame['var2'] = self.positions[:,0]*1e10
                #frame['var3'] = self.positions[:,1]*1e10
                #frame['var4'] = self.positions[:,2]*1e10
                frame['var2'] = self.positions[:, 0]
                frame['var3'] = self.positions[:, 1]
                frame['var4'] = self.positions[:, 2]

                #save frame
                myfile = open(traj_file,'ab')
                np.savetxt(myfile,frame, fmt = "%s %f8 %f8 %f8", )
                myfile.close()
                myfile = open(traj_file,'a')
                myfile.write(string1)
                myfile.close()
                
                # save Energy       
                np.savetxt(Energy_file, Energy)

                print("Energy Converged")
                return 
            
            
        # save Energy       
        np.savetxt(Energy_file, Energy)

        print("Maximum Number of Steps reached")
        return
      
      
    def radial_distribution(self, distancematrix, dr, rmax, numb_of_probes, rho0):

        """Computes the radial distribution

        rdf(R) = { H(R) / [num * V(R)] } / rho0
        for more information see:
        https://de.wikipedia.org/wiki/Radiale_Verteilungsfunktion

        Parameters
        ----------------
        distancematrix: numpy.array
            Array with the distance of each particle of type A to each particle of type B

        dr: float
            thickness of spherical shell

        rmax: float
            maximum radius that will be considered

        numb_of_probes: int
            number of particles of type A, that are considered within the distancematrix

        rho0: float
            ideal gas density (= number of particles / box volume)
        Returns
        -----------------

        histlist: 1 x m numpy.array
            vector containing the radial distribution
        """

        histlist=np.zeros(int(np.ceil(rmax/dr)))
        volume_factor = 4 * np.pi / 3

        for iteration in range (0,len(histlist)):
            histlist[iteration] = len(distancematrix[distancematrix < iteration*dr + dr]) \
                                  - len(distancematrix[distancematrix <= iteration*dr])

            binvolume = volume_factor * ((iteration * dr + dr) ** 3 - (iteration * dr) ** 3)
            histlist[iteration] /= binvolume

        histlist /= (numb_of_probes * rho0)

        return histlist
