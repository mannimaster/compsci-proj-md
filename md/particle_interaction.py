from abc import ABCMeta, abstractmethod, abstractproperty
from boxvectors import directions as directions
import numpy as np
import time
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


class  coulomb(__particle_interaction):
    '''
    class to compute coulomb interaction
    the class offers methods to calculate
        potentials
        total energy
        forces
    '''

    def __init__(self, n_boxes_short_range, L, p_error, **kwargs):
        '''
        creates a coloumb object with properties that do not change over time

        description is missing

        '''

        self.epsilon0 = kwargs.pop('epsilon0', epsilon_0)
        k_cut = kwargs.pop('k_cut', 4 * p_error / (float)(L[0]))

        self.std = np.sqrt(2 * p_error) / (float)(k_cut)
        self.n_boxes_short_range = n_boxes_short_range
        self.constant = 1 / (8 * np.pi * self.epsilon0)        # prefactor for the short range potential/forces
        self.volume = L[0] ** 3                            # Volume of box
        self.k_list = self.__create_k_list(k_cut,L[0])
        return

    def compute_optimal_cutoff(self, Positions, d_Pos, Labels, L, p_error):

        R_cut = (L[0] / (float)(2))
        K_cut = 2 * p_error / (float)(R_cut)

        start_time = time.time()
        self.__short_range_forces(d_Pos, Labels, L)
        T_r = time.time() - start_time

        start_time = time.time()
        self.__long_range_forces(d_Pos, Labels)
        T_k = time.time() - start_time

        factor = 8 * np.pi * Positions.shape[1] ** 2 * R_cut ** 3 / (float)(self.volume ** 2 * K_cut ** 3)

        R_opt_cut = np.sqrt(p_error / (float)(np.pi)) * (factor * T_k / (float)(T_r)) ** (1 / 6.0) * (L[0] / (Positions.shape[1] ** (1 / 6.0)))
        K_opt_cut = 2 * p_error / R_opt_cut

        self.k_list = self.__create_k_list(K_opt_cut,L[0])
        self.std = np.sqrt(2 * p_error) / (float)(K_opt_cut)

        return (R_opt_cut, K_opt_cut, self.std)

    def __create_k_list(self, k_cutoff, boxlength):

        # maximum factor for k (k_max = 2*pi/L *n_max)
        n_max = np.ceil(k_cutoff * boxlength / (float)(2 * np.pi))
        # in dependency of n_max and the box length all linear combination are computed
        lincomb_k = directions(n_max).get_directions() * (2 * np.pi / boxlength)
        # the row k(k1,k2,k3) = [0,0,0] is deleted, as this one is not needed
        k_list = np.delete(lincomb_k, (np.where(np.linalg.norm(lincomb_k, axis=1) == 0)), axis=0)
        # all k-vectors (rows) that are longer than k_cut are deleted
        k_list = np.delete(k_list, (np.where(np.linalg.norm(k_list, axis=1) > k_cutoff)), axis=0)

        if any(np.equal(k_list, [0, 0, 0]).all(1)):
            raise ValueError('The list for the k_values contains the prohibited row [0,0,0]')

        return k_list

    def compute_potential(self,labels, positions, neighbours, distances, r_s, r_c):
        """
        short range potential
        compute_potential(self,labels,positions, neighbours, distances)

        Calculates the short range potential

        Parameters
        ----------
        labels: Nx3 Array
            Array with N rows and ? Columns. The third Column should contain labels, that specify the chemical species of the Particles.
            Particle A should have the label 0 and Particle B should have the label 1. The first column contains the masses, the second the charge.

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
        return self.__short_range_potential(labels, neighbours, distances, r_s, r_c) + self.__long_range_potential(labels[:,1],positions)


    def __switchFunction(self, x, r_0, r_c):
        """
        switchFunction(self,x,r_0,r_c):
        Calculates a switch function that is for calculating the Potential.

        x1 = (x-r_0)/(r_c-r_0)
        return 2x1**2 -3x**2 +1

        Parameters
        ----------
        x : scalar
            absolute distance between two particles

        r_0 : scalar
            switch radius

        r_c : scalar
            cutoff radius

        return : scalar
            value between 1 and 0
        """
        rStar = (x - r_0) / (r_c - r_0)
        rStarSquare = rStar**2
        rStarSquareSquare = rStarSquare**2
        return -6 * rStar * rStarSquareSquare + 15 * rStarSquareSquare - 10*rStarSquareSquare*rStarSquare + 1


    def __short_range_potential(self, labels, neighbours, distances, r_s, r_c):#distances should have the same format/order as the neighborlist
        """
        short range Coulomb potential
        __short_range_potential(self, positions,  labels, n_particles, neighbours, distances,)

        Calculates the short range coulomb potential using a neighbor list.

        Parameters
        ----------
        labels: Nx3 Array
            Array with N rows and ? Columns. The third Column should contain labels, that specify the chemical species of the Particles.
            Particle A should have the label 0 and Particle B should have the label 1. The first column contains the masses, the second the charge.

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

                if absDistance < r_s:
                    shortPotential[i] += labels[j, 1] / absDistance * erfc(absDistance / (np.sqrt(2) * self.std))       #calculating and summing the short range coulomb potential
                else:
                    shortPotential[i] += self.__switchFunction(absDistance, r_s, r_c) * (labels[j, 1] / absDistance * erfc(absDistance / (np.sqrt(2) * self.std)))

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
        #matrix = np.zeros((self.k_list.shape[1],positions.shape[1]))
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
        return_vector = np.multiply(return_vector, 1 / (float)(self.volume * self.epsilon0))

        return return_vector.real

    def compute_energy(self,labels,positions, neighbours, distances, r_s, r_c):

        ''' Calculate the total coloumb energy

                Parameters
                ---------------
                labels: Nx3 Array
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

        return self.__short_range_energy(labels, neighbours, distances, r_s, r_c) + self.__long_range_energy(labels,positions)

    def __short_range_energy(self, labels, neighbours, distances, r_s, r_c):
        
        ''' Calculate the total short range energy of the system

                Parameters
                ---------------
                labels: Nx3 Array
                    Array with N rows and ? Columns. The third Column should contain labels, that specify the chemical species of the Particles.
                    Particle A should have the label 0 and Particle B should have the label 1. The first column contains the masses, the second the charge.

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
        short_range_potential = self.__short_range_potential(labels, neighbours, distances, r_s, r_c)

        return 0.5 * np.sum(np.multiply(short_range_potential,labels[:,1]))
    
    

    def __long_range_energy(self,labels,positions):

        ''' Calculate the total long range energy of the system

                Parameters
                ---------------
                labels: Nx? Array
                    Array with N rows and ? Columns. The third Column should contain labels, that specify the chemical species of the Particles.
                    Particle A should have the label 0 and Particle B should have the label 1. The first column contains the masses, the second the charge.

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
        self_energy = np.sum(np.array(labels[:,1]) ** 2) / (float)(2 * self.epsilon0 * self.std * np.power(2 * np.pi, 1.5))

        return 0.5 * np.sum(np.multiply(long_range_potential,labels[:,1])) - self_energy

    def compute_forces(self,Positions,d_Pos, Labels,L):
        '''

        please add here description

        '''
        
        Coulumb_forces = self.__short_range_forces(Positions, Labels,L) + self.__long_range_forces(d_Pos, Labels)

        return Coulumb_forces
      

    def __short_range_forces(self,d_Pos, Labels,L):

        ''' Calculate the Force resulting from the short range coulomb interaction between the Particles

        Parameters
        ---------------

        d_Pos: NxNx3 Array
            Array with N rows and 3 Columns, each row i contains the x,y and z coordinates of particle i.

        R: Nx1 Array
            Array with N entries. Contains each Particles Distance to the coordinate origin.

        Labels: Nx3 Array
            Array with N rows and ? Columns. The third Column should contain labels, that specify the chemical species of the Particles.
            Particle A should have the label 0 and Particle B should have the label 1. The first column contains the masses, the second the charge.

        L:3x1 Array
            Array containg the Dimensions of the Simulation box


        Returns
        --------------

        Force_short_range: Nx3 Array
            Array with N rows and 3 Columns. Each Row i contains the short-range-Force acting upon Particle i componentwise. 

        '''
        K = directions(self.n_boxes_short_range).get_directions()
        K[:,0] *=L[0]
        K[:,1] *=L[1]
        K[:,2] *=L[2]

        N = d_Pos.shape[0]
        Force_short_range = np.zeros((N,3))
        
        
 
        d_Pos_fix = np.zeros(((N-1 )*N,3))
        d_Pos_fix[:,0] = np.delete(d_Pos[:,:,0].reshape(-1),np.where(np.linalg.norm(d_Pos,axis=2).reshape(-1)==0)).reshape((N-1)*N)
        d_Pos_fix[:,1] = np.delete(d_Pos[:,:,1].reshape(-1),np.where(np.linalg.norm(d_Pos,axis=2).reshape(-1)==0)).reshape((N-1)*N)
        d_Pos_fix[:,2] = np.delete(d_Pos[:,:,2].reshape(-1),np.where(np.linalg.norm(d_Pos,axis=2).reshape(-1)==0)).reshape((N-1)*N)

        Pos_Comb = np.zeros((K.shape[0], (N-1)*N,3))
        Pos_Comb[:,:,0]=np.add.outer(K[:,0],d_Pos_fix[:,0]) # x
        Pos_Comb[:,:,1]=np.add.outer(K[:,1],d_Pos_fix[:,1]) # y
        Pos_Comb[:,:,2]=np.add.outer(K[:,2],d_Pos_fix[:,2]) # z


        ind = np.where(np.linalg.norm(d_Pos,axis=2).reshape(-1)==0)[0]
        ind_1 =(np.cumsum(np.ones(ind.shape))*N-N).astype(int)

        Pos_Comb_norm_1 = np.linalg.norm(Pos_Comb, axis=2)
        Pos_Comb_norm_2 = np.power(Pos_Comb_norm_1,2)

        charges = np.multiply.outer(Labels[:,1],np.ones(N)).reshape(-1)
        charges = np.delete(charges, np.where(np.linalg.norm(d_Pos,axis=2).reshape(-1)==0))

        F_short = np.zeros((MD.N,3))
        #x
        F_short[:,0] =np.sum(np.insert(np.sum(charges*Pos_Comb[:,:,0]/Pos_Comb_norm_2
                                              *(erfc( Pos_Comb_norm_1/np.sqrt(2.0) /self.std) /Pos_Comb_norm_1 
                                                +np.sqrt(2.0/np.pi) /self.std *np.exp( -Pos_Comb_norm_2 /2.0 /self.std**2 ) )
                                              ,axis=0)
                                       ,ind_1,0).reshape(MD.N,MD.N) 
                             ,axis=0)
        #y
        F_short[:,1] =np.sum(np.insert(np.sum(charges*Pos_Comb[:,:,1]/Pos_Comb_norm_2
                                              *(erfc( Pos_Comb_norm_1/np.sqrt(2.0) /self.std) /Pos_Comb_norm_1 
                                                +np.sqrt(2.0/np.pi) /self.std *np.exp( -Pos_Comb_norm_2 /2.0 /self.std**2 ) )
                                              ,axis=0)
                                       ,ind_1,0).reshape(MD.N,MD.N) 
                             ,axis=0)
        #z
        F_short[:,2] =np.sum(np.insert(np.sum(charges*Pos_Comb[:,:,2]/Pos_Comb_norm_2
                                              *(erfc( Pos_Comb_norm_1/np.sqrt(2.0) /self.std) /Pos_Comb_norm_1 
                                                +np.sqrt(2.0/np.pi) /self.std *np.exp( -Pos_Comb_norm_2 /2.0 /self.std**2 ) )
                                              ,axis=0)
                                       ,ind_1,0).reshape(N,N) 
                             ,axis=0)
        F_short*=np.outer(Labels[::-1,1],np.ones(3))/(8*np.pi*self.epsilon0)

        F_short
        
        
        return F_short


    def __long_range_forces(self,d_Pos,Labels):
        ''' Calculate the Force resulting from the long range Coulomb interaction between the Particles

        Parameters
        ---------------

        Positions: Nx3 Array
            Array with N rows and 3 Columns, each row i contains the x,y and z coordinates of particle i.

        Labels: Nx3 Array
            Array with N rows and ? Columns. The third Column should contain labels, that specify the chemical species of the Particles.
            Particle A should have the label 0 and Particle B should have the label 1. The first column contains the masses, the second the charge.

            

        Returns
        --------------

        Force_long_range: Nx3 Array
            Array with N rows and 3 Columns. Each Row i contains the long-range-Force acting upon Particle i componentwise. 

        '''
        
        k1 = np.delete(self.k_list, (np.arange((np.shape(self.k_list)[0])/2)), axis=0)
        SPM = np.tensordot(k1,d_Pos,(1,2))

        sum1= (np.tensordot(Labels[:,1],np.sin(SPM), axes =(0,2)))

        k_betqua = np.linalg.norm(k1,axis=1)**2
        k_betqua_right_size = np.outer(k_betqua, np.ones(3))

        k = np.exp(-self.std**2*k_betqua_right_size/2)*k1/k_betqua_right_size
        Force_long_range = np.tensordot(sum1, k, (0,0))*(np.outer(Labels[:,1],np.ones(3))/self.volume/self.epsilon0)*2
        return Force_long_range


class lennard_jones(__particle_interaction):

    def __init__(self):
        return

    def __switchFunction(self,x,r_0,r_c):
        """
        switchFunction(self,x,r_0,r_c):
        Calculates a switch function that is for calculating the Potential.

        x1 = (x-r_0)/(r_c-r_0)
        return 2x1**2 -3x**2 +1

        Parameters
        ----------
        x : scalar
            absolute distance between two particles

        r_0 : scalar
            switch radius

        r_c : scalar
            cutoff radius

        return : scalar
            value between 1 and 0
        """
        rStar       = (x-r_0)/(r_c-r_0)
        rStarSquare = rStar**2
        rStarSquareSquare = rStarSquare**2
        return -6 * rStar * rStarSquareSquare + 15 * rStarSquareSquare - 10*rStarSquareSquare*rStarSquare + 1

    def compute_potential(self, sigma, epsilon, labels, neighbours, distances, r_s, r_c):
        """
        Lennard Jones
        potentialShort(self, sigma, epsilon, labels, neighbours, distances)

        Calculates the short range Lennard Jones potential using a neighbor list.

        Parameters
        ----------
        labels: N3? Array
            Array with N rows and ? Columns. The third Column should contain labels, that specify the chemical species of the Particles.
            Particle A should have the label 0 and Particle B should have the label 1. The first column contains the masses, the second the charge.

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
                if absDistance<r_s:
                    shortPotentialL[i] += 4*epsilon[int(labels[i,2]+labels[j,2])]*(sigmaPoSix**2-sigmaPoSix)                #calculating and summing the Lennard Jones potential
                else:
                    shortPotentialL[i] += self.__switchFunction(absDistance,r_s,r_c)*(4 * epsilon[int(labels[i, 2] + labels[j, 2])] * (sigmaPoSix ** 2 - sigmaPoSix))
        return shortPotentialL

    def compute_energy(self, sigma, epsilon, labels, neighbours, distances, r_s, r_c):
        
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
        energy_LJ = self.compute_potential(sigma, epsilon, labels, neighbours, distances, r_s=r_s, r_c=r_c,)

        return np.sum(energy_LJ)
    
    

    def compute_forces(self,Positions,Sigma, Epsilon, Labels,L, switch_parameter, r_switch, neighbours):
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

        Labels: Nx3 Array
            Array with N rows and ? Columns. The third Column should contain labels, that specify the chemical species of the Particles.
            Particle A should have the label 0 and Particle B should have the label 1. The first column contains the masses, the second the charge.

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
        

        N = np.size(Positions[:,0])
        Force_LJ = np.zeros((N,3))

        for i in np.arange(N):
            #Find the Indices that should be used for sigma and eps (by adding the corresponding labels)
            index_LJ = (Labels[i,2]*np.ones(N) +Labels[:,2]).astype('int')

            #Create Arrays that contain the approriate Values for sigma and eps, depending on the interaction pair.
            sig6 =(Sigma*np.ones( (N,3) ) )[0,index_LJ]**6
            eps = (Epsilon*np.ones( (N,3) ) )[0,index_LJ]

            Positions_Difference = Positions - Positions[i,:] 
            
            Positions_Difference[:,0] = np.fmod(Positions_Difference[:,0],(L[0]/2))
            Positions_Difference[:,1] = np.fmod(Positions_Difference[:,1],(L[1]/2))
            Positions_Difference[:,2] = np.fmod(Positions_Difference[:,2],(L[2]/2))
            
            #Calculate the Distances
            dist = np.linalg.norm(Positions_Difference, axis = 1)
            dist_8 = dist**8
            dist_2 = dist**2
            for j in neighbours[i]:

                #Calculate the Value of the switch Function
                Switch = (switch_parameter[0]+switch_parameter[1]*dist[j]+switch_parameter[2]*dist_2[j]+switch_parameter[3]*dist[j]*dist_2[j])
                #Calculate the derivate of the switch Function
                d_Switch = (switch_parameter[1]/(2*np.sqrt(dist[j])) +switch_parameter[2] + 1.5*switch_parameter[3]*np.sqrt(dist[j]))

                #sigma to distance ratio
                sig6_dist_ratio = sig6[j] /dist[j]**6

                if dist[j] < r_switch:
                    #This is the Analytical Expression for the LJ force
                    Force_LJ[i,:] += 48.0 *eps[j] *sig6[j] *Positions_Difference[j,:] /dist_8[j] *(-sig6_dist_ratio +0.5)
                else:
                    #This is the derivative of the switch-function, used to cut the LJ-Potential
                    Force_LJ[i,:] += 8.0 *eps[j] *Positions_Difference[i,:] *sig6_dist_ratio*(
                    (6 /dist_2[j]*(-sig6_dist_ratio +0.5)*Switch ) +(d_Switch)*( sig6_dist_ratio-1)) 

        return Force_LJ
     
