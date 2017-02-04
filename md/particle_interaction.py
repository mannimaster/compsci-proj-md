from abc import ABCMeta, abstractmethod, abstractproperty
from boxvectors import directions as directions
import numpy as np
from scipy.special import erfc
from scipy.constants import epsilon_0

class __particle_interaction(object):

    '''
    Abstract class for implementing particle interaction
    '''

    __metaclass__ = ABCMeta

    def __init__(self):
        raise NotImplementedError('You cannot create a particle_interaction object')

    # The parameter input is missing
    @abstractmethod
    def compute_potential(self):
        raise NotImplementedError('The method is not implemented for the abstract class')


    @abstractmethod
    def compute_forces(self):
        raise NotImplementedError('The method is not implemented for the abstract class')

    @abstractmethod
    def compute_energy(self):
        raise NotImplementedError('The method is not implemented for the abstract class')

    # @ MARCO
    # yes this works in the child classes, but it is not necessary like this; you just can call it directly, as you do
    # not any calculations, just call another function.

    def LinComb(self,n):
        return directions(n).get_directions()
    


class  coulomb(__particle_interaction):
    '''
    class to compute coulomb interaction
    the class offers methods to calculate
        potentials
        total energy
        forces
    '''

    def __init__(self,std, n_boxes_short_range,L, k_max_long_range,k_cut):
        '''
        creates a coloumb object with properties that do not change over time

        description is missing

        '''

        self.std = std
        self.n_boxes_short_range = n_boxes_short_range
        self.constant = 1 / (8 * np.pi * epsilon_0)        # prefactor for the short range potential/forces
        self.volume = L[0] ** 3                            # Volume of box

        #in dependency of k_max_long_range and the box length L all linear combination are computed
        lincomb_k = directions(k_max_long_range).get_directions() * (2 * np.pi / L)
        #the row k(k1,k2,k3) = [0,0,0] is deleted, as this one is not needed
        lincomb_k = np.delete(lincomb_k, (len(lincomb_k) - 1) / 2, axis=0)
        #all k-vectors (rows) that are longer than k_cut are deleted
        self.k_list = np.delete(lincomb_k, (np.where(np.linalg.norm(lincomb_k, axis=1) > k_cut)), axis=0)
        
        return 


    def compute_potential(self,labels, positions, neighbours, distances):
        """
        short range potential
        compute_potential(self,labels,positions, neighbours, distances)

        Calculates the short range potential

        Parameters
        ----------
        labels: Nx? Array
            Array with N rows and ? Columns. The third Column should contain labels, that specify the chemical species of the Particles.
            Particle A should have the label 1 and Particle B should have the label 0. The first column contains the masses, the second the charge.

        neighbours : dictionary of lists
            all neighbours within given cutoff radius + skin radius

        distances : like neighbours
            All distances of all neighbours within given cutoff radius + skin radius. The position of the distance in the list is the same as its index in neighbours.

        positions: N x 3 Array
            Array with N rows, where each row represents the position of the particle with the same index

        Returns
        -------
        shortPotential : 1D np.array
            ..the short range potential at the position of each particle. The potentials within the array are in the same order as the positions in the positions array.
        """
        return self.__short_range_potential(labels, neighbours, distances) + self.__long_range_potential(labels[:,1],positions)


    def __short_range_potential(self, labels, neighbours, distances):#distances should have the same format/order as the neighborlist
        """
        short range Coulomb potential
        __short_range_potential(self, positions,  labels, n_particles, neighbours, distances,)

        Calculates the short range coulomb potential using a neighbor list.

        Parameters
        ----------
        labels: Nx? Array
            Array with N rows and ? Columns. The third Column should contain labels, that specify the chemical species of the Particles.
            Particle A should have the label 1 and Particle B should have the label 0. The first column contains the masses, the second the charge.

        neighbours : dictionary of lists
            all neighbours within given cutoff radius + skin radius

        distances : like neighbours
            All distances of all neighbours within given cutoff radius + skin radius. The position of the distance in the list is the same as its index in neighbours.

        Returns
        -------
        shortPotential : 1D np.array
            ..the short range potential at the position of each particle. The potentials within the array are in the same order as the positions in the positions array.
        """
        n_particles = np.shape(labels)[0]
        shortPotential = np.zeros(n_particles,dtype=float)

        for i in neighbours:
            for j,absDistance in zip(neighbours[i],distances[i]):#TIMEPROBLEM http://stackoverflow.com/questions/1663807/how-can-i-iterate-through-two-lists-in-parallel-in-python

                shortPotential[i] += labels[j,1]/absDistance*erfc(absDistance/(np.sqrt(2)*self.std))               #calculating and summing the short range coulomb potential
        shortPotential *= self.constant
        return shortPotential

    def __long_range_potential(self,charges,positions):

        ''' Calculate the Long Range Potential

                Parameters
                ---------------
                charges: N x 1 Array
                    Array with the charge of each particle. The index+1 indicates the particle number. However the order
                    is not important.

                positions: N x 3 Array
                    Array with N rows, where each row represents the position of the particle with the same index

                Returns
                --------------
                vector : N x 1 Array
                    Array with the potential at each particle position
        '''
        
        # initializes all necessary arrays, which will be reused => COMMENT (@daviezz): I Would delete them here. They get over written anyway.
        #matrix = np.zeros((self.k_list.shape[1],len(positions)))
        #vector = np.zeros(self.k_list.shape[1])
        #return_vector = np.zeros(positions.shape[1])

        # compute 1 x K vector of |k|^2 =
        #
        # [
        # k11^2 + k21^2 + k31^2 ,
        # k12^2 + k22^2 + k32^2 ,
        # k13^2 + k23^2 + k33^2 ,
        # ...
        # ]
        # with ki = (k1i k2i k3i)
        # each row is like going through k

        vector = np.sum(np.square(self.k_list), axis = 1)

        # compute 1 x K vector of the Gaussian part exp(-sigma^2 * |k|^2) / |k|^2
        #
        # [
        # exp(-sigma^2 * (k11^2 + k21^2 + k31^2)) / (k11^2 + k21^2 + k31^2) ,
        # exp(-sigma^2 * (k12^2 + k22^2 + k32^2)) / (k12^2 + k22^2 + k32^2) ,
        # exp(-sigma^2 * (k13^2 + k23^2 + k33^2)) / (k13^2 + k23^2 + k33^2) ,
        # ...
        # ]
        # with ki = (k1i k2i k3i)
        # each row is like going through k

        vector = np.divide(np.exp(-np.multiply(vector, self.std ** 2 / (float)(2))), vector)

        # compute K x R matrix of the scalar <k,r> =
        #
        # [
        # [k1r1, k1r2, k1r3, ...] ,
        # [k2r1, k2r2, k2r3, ...] ,
        # [k3r1, k3r2, k3r3, ...] ,
        # [...]
        # ]
        # with ri = (rxi, ryi, rzi) and ki = (k1i k2i k3i)
        # going through the rows is like going through the k values
        # going through the columns is like going through the positions r

        matrix = np.dot(self.k_list, np.transpose(positions))

        # compute 1 x K vector of the elementwise product of the Gaussian part and the structural factor {elementwise product}
        #
        # [
        # k1r1q1 + k1r2q2 + k1r3q3 + ... ,
        # k2r1q1 + k2r2q2 + k2r3q3 + ... ,
        # k3r1q1 + k3r2q2 + k3r3q3 + ... ,
        # ...
        # ]
        #
        # *
        #
        # [
        # (exp(-sigma^2 * (k11^2 + k21^2 + k31^2)) / (k11^2 + k21^2 + k31^2)) * (k1r1q1 + k1r2q2 + k1r3q3 + ...) ,
        # (exp(-sigma^2 * (k12^2 + k22^2 + k32^2)) / (k12^2 + k22^2 + k32^2)) * (k2r1q1 + k2r2q2 + k2r3q3 + ...) ,
        # (exp(-sigma^2 * (k13^2 + k23^2 + k33^2)) / (k13^2 + k23^2 + k33^2)) * (k3r1q1 + k3r2q2 + k3r3q3 + ...) ,
        # ...
        # ]
        # with ki = (k1i k2i k3i) and qi = [q1, q2, q3, ...]
        # each row is like going through k

        vector = np.multiply(vector, np.dot(np.exp(-1j * matrix), charges))

        # compute 1 x R vector of the potential phi = transpose[matrix] * elementwise product {inner product}
        #
        # [
        # [k1r1, k2r1, k3r1, ...] ,
        # [k1r2, k2r2, k3r2, ...] ,
        # [k1r3, k2r3, k3r3, ...] ,
        # [...]
        # ]
        #
        # *
        #
        # [
        # (exp(-sigma^2 * (k11^2 + k21^2 + k31^2)) / (k11^2 + k21^2 + k31^2)) * (k1r1q1 + k1r2q2 + k1r3q3 + ...) ,
        # (exp(-sigma^2 * (k12^2 + k22^2 + k32^2)) / (k12^2 + k22^2 + k32^2)) * (k2r1q1 + k2r2q2 + k2r3q3 + ...) ,
        # (exp(-sigma^2 * (k13^2 + k23^2 + k33^2)) / (k13^2 + k23^2 + k33^2)) * (k3r1q1 + k3r2q2 + k3r3q3 + ...) ,
        # ...
        # ]
        #
        # =
        #
        # [
        # matrix[0,0] * vector[0] + matrix[0,1] * vector[1] + matrix[0,2] * vector[2] + ...
        # matrix[1,0] * vector[0] + matrix[1,1] * vector[1] + matrix[1,2] * vector[2] + ...
        # matrix[2,0] * vector[0] + matrix[2,1] * vector[1] + matrix[2,2] * vector[2] + ...
        # ]

        return_vector = np.dot(np.exp(1j*np.transpose(matrix)), vector)
        return_vector = np.multiply(return_vector, 1 / (float)(self.volume * epsilon_0))

        return return_vector.real

    def compute_energy(self,labels,positions, neighbours, distances):

        ''' Calculate the total coloumb energy

                Parameters
                ---------------
                labels: Nx? Array
                    Array with N rows and ? Columns. The third Column should contain labels, that specify the chemical species of the Particles.
                    Particle A should have the label 1 and Particle B should have the label 0. The first column contains the masses, the second the charge.

                charges: N x 1 Vector
                    Vector with the charge of each particle. The index indicates the particle number.

                positions: N x 3 Array
                    Array with N rows, where each row represents the position of the particle with the same index

                distances : like neighbours
                    All distances of all neighbours within given cutoff radius + skin radius. The position of the distance in the list is the same as its index in neighbours.

                Returns
                --------------
                total energy : float
                    float value of the total coloumb energy
        '''

        return self.__short_range_energy(labels, neighbours, distances) + self.__long_range_energy(labels,positions)

    def __short_range_energy(self, labels, neighbours, distances):
        
        ''' Calculate the total short range energy of the system

                Parameters
                ---------------
                labels: Nx? Array
                    Array with N rows and ? Columns. The third Column should contain labels, that specify the chemical species of the Particles.
                    Particle A should have the label 1 and Particle B should have the label 0. The first column contains the masses, the second the charge.

                charges: N x 1 Vector
                    Vector with the charge of each particle. The index indicates the particle number.

                positions: N x 3 Array
                    Array with N rows, where each row represents the position of the particle with the same index

                Returns
                --------------
                coulomb_short_energy : float
                    float value of the total long range energy
        '''

        #calculates the short range potential
        short_range_potential = self.__short_range_potential(labels, neighbours, distances)

        return 0.5 * np.sum(np.multiply(short_range_potential,labels[:,1]))
    
    

    def __long_range_energy(self,labels,positions):

        ''' Calculate the total long range energy of the system

                Parameters
                ---------------
                labels: Nx? Array
                    Array with N rows and ? Columns. The third Column should contain labels, that specify the chemical species of the Particles.
                    Particle A should have the label 1 and Particle B should have the label 0. The first column contains the masses, the second the charge.

                positions: N x 3 Array
                    Array with N rows, where each row represents the position of the particle with the same index

                Returns
                --------------
                coulomb_long_energy : float
                    float value of the total long range energy
        '''
        
        #calculates the long range potential
        long_range_potential = self.__long_range_potential(labels[:,1],positions)

        # calculates the self-interaction potential
        self_energy = np.sum(np.array(labels[:,1]) ** 2) / (float)(2 * epsilon_0 * self.std * np.power(2 * np.pi, 1.5))

        return 0.5 * np.sum(np.multiply(long_range_potential,labels[:,1])) - self_energy

    def compute_forces(self,Positions,R, Labels,L):
        '''

        please add here description

        '''
        
        Coulumb_forces = self.__short_range_forces(Positions,R, Labels,L) + self.__long_range_forces(Positions, Labels)

        return Coulumb_forces

    
    def __short_range_forces(self,Positions,R, Labels,L):#COMMENT (@mannimaster): Do you intend to use "R" again? else delete it here.
        ''' Calculate the Force resulting from the short range coulomb interaction between the Particles

        Parameters
        ---------------

        Positions: Nx3 Array
            Array with N rows and 3 Columns, each row i contains the x,y and z coordinates of particle i.

        R: Nx1 Array
            Array with N entries. Contains each Particles Distance to the coordinate origin.

        Labels: Nx? Array
            Array with N rows and ? Columns. The third Column should contain labels, that specify the chemical species of the Particles.
            Particle A should have the label 1 and Particle B should have the label 0. The first column contains the masses, the second the charge.

        L:3x1 Array
            Array containg the Dimensions of the Simulation box


        Returns
        --------------

        Force_short_range: Nx3 Array
            Array with N rows and 3 Columns. Each Row i contains the short-range-Force acting upon Particle i componentwise. 

        '''
        K = self.LinComb(self.n_boxes_short_range)
        K[:,0] *=L[0]
        K[:,1] *=L[1]
        K[:,2] *=L[2]

        N = np.size(Positions[:,0])
        Force_short_range = np.zeros((N,3))
        k = np.size(K[:,0])
        charges_pre_delete = np.zeros((N,3))
        charges_pre_delete[:,0]=Labels[:,1]  
        charges_pre_delete[:,1]=Labels[:,1] 
        charges_pre_delete[:,2]=Labels[:,1] 
        for i in np.arange(N):

            dists_single_cell_pre_delete = Positions[i,:]-Positions
            dists_single_cell = np.delete(dists_single_cell_pre_delete,i,0)
            charges = np.delete(charges_pre_delete,i,0)
            #This Skript paralellizes the sum over j and executes the sum over vectors n within the loop
            for j in np.arange(k):
                dists = dists_single_cell+K[j]

                norm_dists = np.linalg.norm(dists)

                Force_short_range[i,:] += np.sum(charges*dists/norm_dists**2 
                *( erfc( norm_dists/np.sqrt(2)/self.std )/norm_dists 
                + np.sqrt(2.0/np.pi)/self.std*np.exp(-norm_dists**2/(2*self.std**2) )) ,0)

        #Getting the Pre-factor right        
        Force_short_range = Force_short_range* charges_pre_delete /(8*np.pi*epsilon_0)
        return Force_short_range


    def __long_range_forces(self,Positions,Labels):
        ''' Calculate the Force resulting from the long range Coulomb interaction between the Particles

        Parameters
        ---------------

        Positions: Nx3 Array
            Array with N rows and 3 Columns, each row i contains the x,y and z coordinates of particle i.

        Labels: Nx? Array
            Array with N rows and ? Columns. The third Column should contain labels, that specify the chemical species of the Particles.
            Particle A should have the label 1 and Particle B should have the label 0. The first column contains the masses, the second the charge.

            

        Returns
        --------------

        Force_long_range: Nx3 Array
            Array with N rows and 3 Columns. Each Row i contains the long-range-Force acting upon Particle i componentwise. 

        '''
        
        # two k-vectors i,j have property k_i = -k_j respectively, delete one of them
        k = np.delete(self.k_list, (np.arange((np.shape(self.k_list)[0])/2)), axis=0)
        
        # setup data needed for calculation
        k_betqua = np.sum(np.square(k), axis = 1)   # |k|^2
        num_k = np.shape(k_betqua)[0]               # number of k-vectors
        num_Par = np.shape(Positions)[0]            # number of particles
        f_1 = np.zeros((num_k,3))                   # placeholder for Forces during summation over k outside of for second loop
        f_2 = np.zeros((num_k))                     # placeholder for Forces during summation over k in second for-loop
        Force_long_range = np.zeros((num_Par,3))    # placeholder for long ranged forces
        charges = np.zeros((num_Par,3))             # create num_Par x 3 matrix with charges
        charges[:,0]=Labels[:,1]                    # in each column L[:,1] repeated
        charges[:,1]=Labels[:,1] 
        charges[:,2]=Labels[:,1]

        for h in np.arange(num_Par):    # loop over all particles

            for j in np.arange(num_k):  # loop over all k-vectors
                
                # "structure factor" sum (right sum in equation)
                f_2[j] = sum(Labels[:,1]*np.sin(np.dot(k[j,:],(Positions[h,:]*np.ones((num_Par,3)) - Positions).transpose())))
            
            # complete sum over all k, (left part of equation)
            f_1 = (((np.exp(-self.std**2/2*k_betqua)/k_betqua)*f_2)*np.ones((3,num_k))).transpose()*k  
            
            # actually sum over all k
            Force_long_range[h,:] = sum(f_1)   
        
        # get prefactor right                                      
        Force_long_range *= charges/(self.volume*epsilon_0)*2      # multiply by 2 because of symmetry properties of sine
        return Force_long_range


class lennard_jones(__particle_interaction):

    def __init__(self):
        self.constant = 1 / (8 * np.pi * epsilon_0)  #prefactor for the short range potential/forces
        return


    def compute_potential(self, sigma, epsilon, labels, neighbours, distances):
        """
        Lennard Jones
        potentialShort(self, sigma, epsilon, labels, neighbours, distances)

        Calculates the short range Lennard Jones potential using a neighbor list.

        Parameters
        ----------
        labels: Nx? Array
            Array with N rows and ? Columns. The third Column should contain labels, that specify the chemical species of the Particles.
            Particle A should have the label 1 and Particle B should have the label 0. The first column contains the masses, the second the charge.

        neighbours : dictionary of lists
            all neighbours within given cutoff radius + skin radius

        distances : like neighbours
            All distances of all neighbours within given cutoff radius + skin radius. The position of the distance in the list is the same as its index in neighbours.

        sigma : 1 dim np.array, float
            Contains the sigmas for the Lennard Jones potential.

        epsilon : 1 dim np.array, float
            Contains the epsilons of the Lennard Jones potential

            Returns
        -------
        shortPotentialC : 1D np.array
                ..the short range potential at the position of each particle. The potentials within the array are in the same order as the positions in the positions array.
        """
        n_particles = np.shape(labels)[0]
        shortPotentialL = np.zeros(n_particles,dtype=float)

        for i in neighbours:
            for j,absDistance in zip(neighbours[i],distances[i]):#TIMEPROBLEM http://stackoverflow.com/questions/1663807/how-can-i-iterate-through-two-lists-in-parallel-in-python
                #print i,j,absDistance,labels[i,2]+labels[j,2]#, (sigma[labels[i,2]+labels[j,2]]/absDistance)**6
                sigmaPoSix = (sigma[int(labels[i,2]+labels[j,2])]/absDistance)**6                                       #precalulating the sixth power
                shortPotentialL[i] += 4*epsilon[int(labels[i,2]+labels[j,2])]*(sigmaPoSix**2-sigmaPoSix)                #calculating and summing the Lennard Jones potential
        shortPotentialL *= self.constant
        return shortPotentialL

    def compute_energy(self, sigma, epsilon, labels, neighbours, distances):
        
        ''' Calculate the total Lennard-Jones energy

        Parameters
        ---------------
        labels: Nx? Array
            Array with N rows and ? Columns. The third Column should contain labels, that specify the chemical species of the Particles.
            Particle A should have the label 1 and Particle B should have the label 0. The first column contains the masses, the second the charge.

        neighbours : dictionary of lists
            all neighbours within given cutoff radius + skin radius

        distances : like neighbours
            All distances of all neighbours within given cutoff radius + skin radius. The position of the distance in the list is the same as its index in neighbours.

        sigma : 1 dim np.array, float
            Contains the sigmas for the Lennard Jones potential.

        epsilon : 1 dim np.array, float
            Contains the epsilons of the Lennard Jones potential

        Returns
        --------------
        total energy : float
            float value of the total Lennard-Jones energy
        '''
        energy_LJ = self.compute_potential(sigma, epsilon, labels, neighbours, distances)

        return np.sum(np.multiply(energy_LJ,labels[:,0]))
    
    

    def compute_forces(self,Positions,R,Sigma, Epsilon, Labels,L, switch_parameter, r_switch, neighbours):
        ''' Calculate the Force resulting from the lennard Jones Interaction between the Particles

        Parameters
        ---------------

        Positions: Nx3 Array
            Array with N rows and 3 Columns, each row i contains the x,y and z coordinates of particle i.

        R: Nx1 Array
            Array with N entries. Contains each Particles Distance to the coordinate origin.

        Sigma: 3x1 Array
            Array with 3 entries. Regarding 2 Particles the entries are: Sigma = [sigma_AA, simga_AB, sigma_BB].

        Epsilon: 3x1 Array
            Array with 3 entries. Regarding 2 Particles the entries are: Epsilon = [epsilon_AA, epsilon_AB, epsilon_BB].

        Labels: Nx? Array
            Array with N rows and ? Columns. The third Column should contain labels, that specify the chemical species of the Particles.
            Particle A should have the label 1 and Particle B should have the label 0. The first column contains the masses, the second the charge.

        switch_parameter: 4x1 Array
            Array that contains the coefficient values for the switch polynomial

        r_switch: float
            Distance where the switch function kicks in. 

        neighbours: dictionary of lists
            all neighbours within given cutoff radius + skin radius
        Returns
        --------------

        Force_LJ: Nx3 Array
            Array with N rows and 3 Columns. Each Row i contains the LJ-Force acting upon Particle i componentwise. 

        '''
        #Get Neighbour-List, naive_Dists should be replaced in the future by Nils' Work
        

        N = np.size(R)
        Force_LJ = np.zeros((N,3))

        for i in np.arange(N):
            #Find the Indices that should be used for sigma and eps (by adding the corresponding labels)
            index_LJ = (Labels[i,2]*np.ones(N) +Labels[:,2]).astype('int')

            #Create Arrays that contain the approriate Values for sigma and eps, depending on the interaction pair.
            sig6 =(Sigma*np.ones( (N,3) ) )[0,index_LJ]**6
            eps = (Epsilon*np.ones( (N,3) ) )[0,index_LJ]

            Positions_Difference = Positions[i,:] - Positions
            
            Positions_Difference[:,0] = Positions_Difference[:,0]%(L[0]/2)
            Positions_Difference[:,1] = Positions_Difference[:,1]%(L[1]/2)
            Positions_Difference[:,2] = Positions_Difference[:,2]%(L[2]/2)
            for j in neighbours[i]:
                #Calculate the Distances 
                dist = np.linalg.norm(Positions_Difference[j,:])

                #Calculate the Value of the switch Function
                Switch = (switch_parameter[0]+switch_parameter[1]*dist+switch_parameter[2]*dist**2+switch_parameter[3]*dist**3)
                #Calculate the derivate of the switch Function
                d_Switch = (switch_parameter[1]/(2*np.sqrt(dist)) +switch_parameter[2] + 1.5*switch_parameter[3]*np.sqrt(dist))

                #sigma to distance ratio
                sig6_dist_ratio = sig6[j] /dist**6

                if dist < r_switch:
                    #This is the Analytical Expression for the LJ force
                    Force_LJ[i,:] += 48.0 *eps[j] *sig6[j] *Positions_Difference[j,:] /dist**8 *(-sig6_dist_ratio +0.5)
                else:
                    #This is the derivative of the switch-function, used to cut the LJ-Potential
                    Force_LJ[i,:] += 8.0 *eps[j] *Positions_Difference[i,:] *sig6_dist_ratio*(
                    (6 /dist**2*(-sig6_dist_ratio +0.5)*Switch ) +(d_Switch)*( sig6_dist_ratio-1)) 

        return Force_LJ
     
