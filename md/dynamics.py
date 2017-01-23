class dynamics(object):

    def __init__(self):
        return

    def __velocity_verlet_integrator(self,Positions,Velocities, Forces, Labels,Sigma, Epsilon ,dt):
        ''' The Verlocity Verlet Integrator
        '''
        N = np.size(Positions[:,0])
        R = np.sqrt(np.sum(Positions**2,1))

        Forces_old = Forces

        Positions_new = Positions +Velocities*dt + Forces_old/(np.outer(Labels[:,0],np.ones(3))) *dt**2 

        #Implement PBC
        Positions_new[:,0] = Positions_new[:,0]%L[0]
        Positions_new[:,1] = Positions_new[:,0]%L[1]
        Positions_new[:,2] = Positions_new[:,0]%L[2]

        # The function to calculate forces will change in the future and will include Coulumb interaciton
        Forces_new = Lennard_Jones_Force(Positions_new,R,Sigma, Epsilon, Labels)

        Velocities_new = Velocities + (Forces_old+Forces_new)/(2*(np.outer(Labels[:,0],np.ones(3))))*dt

        return Positions_new, Velocities_new, Forces_New
    
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
        CM = np.sum( np.reshape( np.repeat(m,3) ,(N,3) )*Velocities, 0)/M #Velocity of Center of Mass: (Sum_over_i (m_i*v_i) )/M
        internal_Velocities = Velocities-CM #remove components along external degrees of freedom
        Temperature = np.sum(m*np.linalg.norm(internal_Velocities,axis = 1)**2 )/1.38064852e-23/(3*N-3) #Calculate the Actual Temperature
        # T = 2/kB/N_df * <K>
        #<K> = Sum_over_i (m_i*v_i**2)/2
        # N_df = 3*N-3
        return Temperature

    def compute_dynamics(self,Positions,Velocities,Forces,Labels,Sigma, Epsilon ,dt, Steps):
        N = np.size(Positions[:,0])
        Trajectory = np.zeros((N,3,Steps+1))
        Trajectory[:,:,0] = Positions #initial Positions
        for i in np.arange(Steps):
        """Propagates the System using Velocity Verlet Integrator and Andersen Thermostat"""
            # Calculate new Positions and Forces
            Positions_new, Velocities_new, Forces_new = self.Velocity_Verlet(Positions,Velocities,Forces,Labels,Sigma, Epsilon ,dt)
            # Update Trajectory frame, Positions and Velocites
            Trajectory[:,:,i+1] = Positions_new
            Positions = Positions_new   
            Velocities = Velocities_new

            #Andersen Thermostat
            Rand = np.random.uniform(size =N) #Draw Random Number for every Particle
            indexes = np.where(Rand<p_rea) #Check wich random numbers are smaller than the reassingment probability p_rea
            if np.size(indexes) is not 0: 
                Velocities[indexes] = sample_maxwell_boltzmann(N = np.size(indexes), m = m[indexes], T=T)
        #Reassign a new Velocity to the Correspoding Particles
        
    
        return Trajectory
