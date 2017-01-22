import numpy as np
from abc import ABCMeta, abstractmethod, abstractproperty


class __distribution(object):
__metaclass__ = ABCMeta


    def __init__(self):
        raise NotImplementedError('You cannot create a distribution object')
        return

    @abstractmethod
    def sample_distribution(self):
        pass


class maxwellboltzmann(__distribution):

    def __init__(self):
        #probably implementable, otherwise a NotImplementedError as well ?
        return

    def sample_distribution(self,N,m,T,kB = 1.38064852e-23):
         ''' Sample velocities from the Maxwell-Boltzmann Distribution for a given chemical species


        Parameters
        ------------
        N: int
            Number of Particles

        m: Nx1 Array, or float
            mass of Particles, if m is a float, then all particles are assumed to have the same mass. 

        T: float
            Temperature in Kelvin

        kB: float
            Boltzmann Constant  = 1.38064852e-23

        Returns
        -----------
        v: N dimensional Array
            contains the velocities of each particle
        '''
        sigma = kB*T/m
        Velocities = np.zeros((N,3))

        v_x = np.sqrt(sigma)*np.random.normal(size=N)
        v_y = np.sqrt(sigma)*np.random.normal(size=N)
        v_z = np.sqrt(sigma)*np.random.normal(size=N)

    #     v_total = np.sqrt(v_x**2+v_y**2+v_z**2)

        Velocities[:,0]= v_x
        Velocities[:,1]= v_y
        Velocities[:,2]= v_z

        return Velocities
