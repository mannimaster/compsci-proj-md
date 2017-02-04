import numpy as np
from particle_interaction import coulomb
from particle_interaction import lennard_jones
from distribution import maxwellboltzmann
from neighbourlist import neighbourlist

class dynamics(object):

    def __init__(self):
        return

    def __velocity_verlet_integrator(self,
                                     Positions,
                                     Velocities, 
                                     Forces, 
                                     Labels,
                                     Sigma, 
                                     Epsilon ,
                                     dt,
                                     L, 
                                     std, 
                                     n_boxes_short_range,
                                     k_max_long_range,
                                     switch_parameter, 
                                     r_switch,
                                     k_cut,
                                     r_cut_coulomb):
        ''' The Verlocity Verlet Integrator
        '''
        
        N = np.size(Positions[:,0])
        R = np.sqrt(np.sum(Positions**2,1))
        

        Forces_old = Forces

        Positions_new = Positions +Velocities*dt + Forces_old/(np.outer(Labels[:,0],np.ones(3))) *dt**2 

        #Implement PBC
        Positions_new[:,0] = Positions_new[:,0]%L[0]
        Positions_new[:,1] = Positions_new[:,1]%L[1]
        Positions_new[:,2] = Positions_new[:,2]%L[2]
        
        #Compute neighbourlist
        neighbours, distances = neighbourlist().compute_neighbourlist(Positions_new, L[0], r_cut_coulomb)
        
        Forces_new = coulomb(
            std, 
            n_boxes_short_range,
            L,
            k_max_long_range, 
            k_cut).compute_forces(
            Positions_new, 
            Labels,
            L)+lennard_jones(
            ).compute_forces(
            Positions_new,
            Sigma, 
            Epsilon, 
            Labels,
            L, 
            switch_parameter, 
            r_switch,
            neighbours)
        
        Velocities_new = Velocities + (Forces_old+Forces_new)/(2*(np.outer(Labels[:,0],np.ones(3))))*dt

        return Positions_new, Velocities_new, Forces_new
    
    def Thermometer(self, Labels, Velocities):
        
        """ Compute the Instantaneos Temperature of a given Frame
        
        Parameters
        ------------
        Labels: Nx3 Array
            Array with one Row for each particle, containing masses, charges and chemical label repsectively
            
        Velocities: Nx3 Array
            Array with one Row for each particle, containing its velocities componentwise. 
           
        Returns
        -----------
        Temperatur: float
            The instantaneous Temperature of the System
        """
        m = Labels[:,0] #Masses
        N = np.size(m)
        M = np.sum(m) #Total Mass
        
        #Velocity of Center of Mass: (Sum_over_i (m_i*v_i) )/M
        CM = np.sum( np.reshape( np.repeat(m,3) ,(N,3) )*Velocities, 0)/M
        
        #remove components along external degrees of freedom
        internal_Velocities = Velocities-CM
        
        # T = 2/kB/N_df * <K>
        #<K> = Sum_over_i (m_i*v_i**2)/2
        # N_df = 3*N-3
        Temperature = np.sum(m*np.linalg.norm(internal_Velocities,axis = 1)**2 )/1.38064852e-23/(3*N-3) #Calculate the Actual Temperature

        return Temperature

    def compute_dynamics(self,
                         Positions,
                         Velocities,
                         Forces,
                         Labels,
                         Sigma, 
                         Epsilon ,
                         dt,
                         L,
                         std, 
                         n_boxes_short_range,
                         k_max_long_range,
                         switch_parameter,
                         p_rea,
                         T, 
                         r_switch,
                         k_cut,
                         r_cut_coulomb):
        
        """Propagates the System using Velocity Verlet Integrator and Andersen Thermostat"""
        
        # Calculate new Positions and Forces
        Positions_new, Velocities_new, Forces_new = self.__velocity_verlet_integrator(
            Positions,
            Velocities,
            Forces,
            Labels,
            Sigma, 
            Epsilon,
            dt,
            L,
            std, 
            n_boxes_short_range,
            k_max_long_range, 
            switch_parameter, 
            r_switch,
            k_cut,
            r_cut_coulomb)
        

        #Andersen Thermostat
        N = np.size(Positions[:,0])
        m = Labels[:,0]
        
        #Draw Random Number for every Particle
        Rand = np.random.uniform(size =N) 
        
        #Check wich random numbers are smaller than the reassingment probability p_rea
        indexes = np.where(Rand<p_rea) 
        if np.size(indexes) is not 0: 
            
            #Reassign a new Velocity to the Correspoding Particles
            Velocities[indexes] = maxwellboltzmann().sample_distribution(N = np.size(indexes), m = m[indexes], T=T)
            


        return Positions_new, Velocities_new, Forces_new
