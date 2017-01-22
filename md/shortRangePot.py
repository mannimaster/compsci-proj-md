import numpy as np
import scipy.special as scsp
import scipy.constants as scco
from neighbourlist import neighbourlist as nList

class PotentialShortCalc(object):
    constant = 1/(8*np.pi*scco.epsilon_0)
    
    def __init__(self):
        """
        self.constant = 1/(8*np.pi*scco.epsilon_0)
        """
        
        self.constant = 1/(8*np.pi*scco.epsilon_0)
        return
    
    def potentialShort(self, positions,  variance, labels, n_particles, neighbors, distances, sigma, epsilon):#distances should have the same format/order as the neighborlist
        """
        PROTOTYP
        potentialShort(self, positions,  variance, labels, n_particles, neighbors, distances)
        
        Calculates the short range potential using a neighbor list.
        
        Parameters
        ----------
        positions : np.array or single numbers
            array of the positionvectors with each line a position

        variance : scalar number
            variance of the erfc

        labels : np.array
            first column the masses, second the charge
            
        n_particles : scalar number
            number of particles
            
        neighbors : dictionary of lists
            all neighbors within given cutoff radius + skin radius
            
        distances : like neighbors
            All distances of all neighbors within given cutoff radius + skin radius. The position of the distance in the list is the same as its index in neighbors.
            
        sigma : np.array
            Contains the sigmas for the Lennard Jones potential.
            
        epsilon : np.array
            Contains the epsilons of the Lennard Jones potential
            
        Returns
        -------
        shortPotential : 1D np.array
            ..the short range potential at the position of each particle. The potentials within the array are in the same order as the positions in the positions array.
        """
        shortPotential = np.zeros(n_particles,dtype=float)
        
        for i in neighbors:
            for j,absDistance in zip(neighbors[i],distances[i]):#TIMEPROBLEM http://stackoverflow.com/questions/1663807/how-can-i-iterate-through-two-lists-in-parallel-in-python
                #print i,j,absDistance,labels[i,2]+labels[j,2]#, (sigma[labels[i,2]+labels[j,2]]/absDistance)**6
                sigmaPoSix = (sigma[labels[i,2]+labels[j,2]]/absDistance)**6
                shortPotential[i] += labels[j,1]/absDistance*scsp.erfc(absDistance/(np.sqrt(2)*variance))+4*epsilon[labels[i,2]+labels[j,2]]*(sigmaPoSix**2-sigmaPoSix)
                
        shortPotential *= self.constant
        return shortPotential