from abc import ABCMeta, abstractmethod, abstractproperty
from boxvectors import directions as directions

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


    # if we have some class properties, we can implement it like below
    # @abstractproperty
    # def get_anything(self):

    #if we have a method that is the same for all child classes
    #@classmethod



class coloumb(__particle_interaction):

    def __init__(self):
        raise NotImplementedError('You cannot create a coloumb object, just use the methods')
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


    def compute_forces(self,positions,box):
        # compute forces
        return None


    def __short_range_forces(self,Positions,R,Sigma, Epsilon, Labels, std, L, Directions):
        ''' Calculate the Force resulting from the short range coulomb interaction between the Particles

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
        K = Directions
        K[:,0] *=L[0]
        K[:,1] *=L[1]
        K[:,2] *=L[2]

        N = np.size(Positions[:,0])
        Force_short_range = np.zeros((N,3))
        k = np.size(K[:,0])


        for i in np.arange(N):

            dists_single_cell = Positions[i,:]-Positions
            charges = np.resize( Labels[current_Neighbors,1] , np.shape(dists_single_cell))

            #This Skript paralellizes the sum over j and executes the sum over vectors n within the loop
            for j in np.arange(k):
                dists = dists_single_cell+K[j]

                norm_dists = np.linalg.norm(dists)

                Force_short_range[i,:] += np.sum(charges*dists/norm_dists**2 
                *( erfc( dists/np.sqrt(2)/std )/norm_dists 
                + np.sqrt(2.0/np.pi)/std*np.exp(-norm_dists**2/(2*std**2) )) ,0)

        #Getting the Pre-factor right        
        Force_short_range = Force_short_range* np.resize( Labels[:,1], np.shape(Positions) ) /(8*np.pi*Epsilon_0)
        return Force_short_range


    def __long_range_forces(self):
        # do anything
        return None


class lennard_jones(__particle_interaction):

    def __init__(self):
        raise NotImplementedError('You cannot create a lennard_jones object, just use the classmethods')
        return


    def compute_potential(self,positions,sigma,epsilon):
        # do anything
        # use from super the result of neighbourlist
        return None


    def compute_forces(self,Positions,R,Sigma, Epsilon, Labels, switch_parameter, r_switch):
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
        Neighbors = get_Neighbors(R)

        N = np.size(R)
        Force_LJ = np.zeros((N,3))

        for i in np.arange(N):
            #Find the Indices that should be used for sigma and eps (by adding the corresponding labels)
            index = (Labels[i,2]*np.ones(N) +Labels[:,2]).astype('int')

            #Create Arrays that contain the approriate Values for sigma and eps, depending on the interaction pair.
            sig6 =(Sigma*np.ones( (N,3) ) )[0,index]**6
            eps = (Epsilon*np.ones( (N,3) ) )[0,index]

            Positions_Difference = Positions[i,:] - Positions
            
            Positions_Difference[:,0] = Positions_Difference[:,0]%(L[0]/2)
            Positions_Difference[:,1] = Positions_Difference[:,0]%(L[1]/2)
            Positions_Difference[:,2] = Positions_Difference[:,0]%(L[2]/2)
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
                    Force_LJ[i,:] += 8.0 *eps[j] *Positions_Difference[i,:] *sig6_dist_ratio
                                    *( (6 /dist**2 (-sig6_dist_ratio +0.5)*Switch ) +(d_Switch)*( sig6_dist_ratio-1) 

        return Force_LJ
     
