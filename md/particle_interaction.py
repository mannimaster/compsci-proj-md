from abc import ABCMeta, abstractmethod, abstractproperty
from boxvectors import directions as directions
from neighbourlist import neighbourlist
import numpy as np
from scipy.special import erf
from scipy.special import erfc
from scipy.constants import epsilon_0
'''
Abstract class for implementing particle interaction
'''

class __particle_interaction(object):
    __metaclass__ = ABCMeta

    def __init__(self):
        raise NotImplementedError('You cannot create a particle_interaction object')
        return

    # The parameter input is missing
    @abstractmethod
    def compute_potential(self):
        pass


    @abstractmethod
    def compute_forces(self):
        pass
    # will this work in child classes ?
    def LinComb(self,n):
        return directions(n).get_directions()


    # if we have some class properties, we can implement it like below
    # @abstractproperty
    # def get_anything(self):

    #if we have a method that is the same for all child classes
    #@classmethod



class  coulomb(__particle_interaction):

    def __init__(self,std, n_boxes_short_range,k_max_long_range ):
        self.std = std
        self.n_boxes_short_range =n_boxes_short_range 
        self.k_max_long_range = k_max_long_range
        return


    def compute_potential(self,positions,box):
        # do anything
        return None


    def __short_range_potential(self):
        # do anything
        # use from super the result of neighbourlist
        return None


    def __long_range_potential(self):
        # do anything
        return None
    
    def compute_forces(self,Positions,R, Labels,L):
        Coulumb_forces = self.__short_range_forces(Positions,R, Labels,L) + self.__long_range_forces()
        return Coulumb_forces

    
    def __short_range_forces(self,Positions,R, Labels,L):
        ''' Calculate the Force resulting from the short range coulomb interaction between the Particles

        Parameters
        ---------------

        Positions: Nx3 Array
            Array with N rows and 3 Columns, each row i contains the x,y and z coordinates of particle i.

        R: Nx1 Array
            Array with N entries. Contains each Particles Distance to the coordinate origin.

        Labels: Nx? Array
            Array with N rows and ? Columns. The third Column should contain labels, that specify the chemical species of the Particles.
            Particle A should have the label 1 and Particle B should have the label 0.

        std: float
            Sigma that is mentioned in the coulomb forces.

        L:3x1 Array
            Array containg the Dimensions of the Simulation box

        Directions: Array
            Array that contains Pointers to neighbouring boxes

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
        charges = np.zeros((N,3))
        charges[:,0]=Labels[:,1]  
        charges[:,1]=Labels[:,1] 
        charges[:,2]=Labels[:,1] 
        for i in np.arange(N):

            dists_single_cell = Positions[i,:]-Positions
            #This Skript paralellizes the sum over j and executes the sum over vectors n within the loop
            for j in np.arange(k):
                dists = dists_single_cell+K[j]

                norm_dists = np.linalg.norm(dists)

                Force_short_range[i,:] += np.sum(charges*dists/norm_dists**2 
                *( erfc( dists/np.sqrt(2)/self.std )/norm_dists 
                + np.sqrt(2.0/np.pi)/self.std*np.exp(-norm_dists**2/(2*self.std**2) )) ,0)

        #Getting the Pre-factor right        
        Force_short_range = Force_short_range* charges /(8*np.pi*epsilon_0)
        return Force_short_range


    def __long_range_forces(self,Positions,Labels,L):
        ''' Calculate the Force resulting from the long range Coulomb interaction between the Particles

        Parameters
        ---------------

        Positions: Nx3 Array
            Array with N rows and 3 Columns, each row i contains the x,y and z coordinates of particle i.

        Labels: Nx? Array
            Array with N rows and ? Columns. The third Column should contain labels, that specify the chemical species of the Particles.
            Particle A should have the label 1 and Particle B should have the label 0.

        std: float
            Sigma that is mentioned in the coulomb forces.

        L:3x1 Array
            Array containg the Dimensions of the Simulation box
            
        k_max_long_range : float
            Maximal value of k1,k2,k3 which comes from cutoff-radius in k-space


        Returns
        --------------

        Force_long_range: Nx3 Array
            Array with N rows and 3 Columns. Each Row i contains the long-range-Force acting upon Particle i componentwise. 

        '''
        i = np.complex(0,1)
        
        # setup k-vector matrix
        k = self.LinComb(self.k_max_long_range)
        k *= (2*np.pi)/L[0]
        
        # setup data needed for calculation
        k_betqua = sum(np.transpose(k**2))          # |k|^2
        num_k = np.shape(k_betqua)[0]               # number of k-vectors
        num_Par = np.shape(Positions)[0]            # number of particles
        f_1 = np.zeros((num_k,3))                   # placeholder for Forces during summation over k outside of for loop
        f_2 = np.zeros((num_k))                     # placeholder for Forces during summation over k in for-loop
        Force_long_range = np.zeros((num_Par,3))    # placeholder for long range forces
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
            
            # result with k = (0,0,0) is nan --> set to 0
            f_1[np.isnan(f_1)] = 0   
            
            # actually sum over all k
            Force_long_range[h,:] = sum(f_1)   
        
        # get prefactor right                                      
        Force_long_range *= charges/(L[0]**3*epsilon_0)
        return Force_long_range


class lennard_jones(__particle_interaction):

    def __init__(self):
        return


    def compute_potential(self,positions,sigma,epsilon):
        # do anything
        # use from super the result of neighbourlist
        return None


    def compute_forces(self,Positions,R,Sigma, Epsilon, Labels,L, switch_parameter, r_switch):
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
            Particle A should have the label 0 and Particle B should have the label 1.

        switch_parameter: 4x1 Array
            Array that contains the coefficient values for the switch polynomial

        r_switch: float
            Distance where the switch function kicks in. 

        Returns
        --------------

        Force_LJ: Nx3 Array
            Array with N rows and 3 Columns. Each Row i contains the LJ-Force acting upon Particle i componentwise. 

        '''
        #Get Neighbour-List, naive_Dists should be replaced in the future by Nils' Work
        NB = neighbourlist()
        Neighbors = NB.naive_Dists(R)

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
            for j in Neighbors[i]:
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
                    (6 /dist**2 (-sig6_dist_ratio +0.5)*Switch ) +(d_Switch)*( sig6_dist_ratio-1)) 

            return Force_LJ
     
