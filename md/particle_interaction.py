from abc import ABCMeta, abstractmethod, abstractproperty

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


    def __short_range_forces(self):
        # do anything
        # use from super the result of neighbourlist
        return None


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


    def compute_forces(self,positions,sigma,epsilon):
        # compute forces
        # use from super the result of neighbourlist
        return None