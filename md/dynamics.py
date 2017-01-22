class dynamics(object):

    def __init__(self):
        raise NotImplementedError('You cannot create a dynamics object, just use the classmethods')
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

    def compute_dynamics(self, Steps):
        Trajectory = np.zeros((N,3,Steps))
        Trajectory[:,:,0] = Positions
        for i in np.arange(Steps):
        """Propagates the System using Velocity Verlet Integrator and Andersen Thermostat"""
            # Calculate new Positions and Forces
            Positions_new, Velocities_new, Forces_new = Velocity_Verlet(Positions,Velocities,Forces,Labels,Sigma, Epsilon ,dt)
            # Update Trajectory frame, Positions and Velocites
            Trajectory[:,:,i] = Positions_new
            Positions = Positions_new   
            Velocities = Velocities_new

            #Andersen Thermostat
            Rand = np.random.uniform(size =N) #Draw Random Number for every Particle
            indexes = np.where(Rand<p_rea) #Check wich random numbers are smaller than the reassingment probability p_rea
            if np.size(indexes) is not 0: 
                Velocities[indexes] = sample_maxwell_boltzmann(N = np.size(indexes), m = m[indexes], T=T)
        #Reassign a new Velocity to the Correspoding Particles
        
    
        return Trajectory
