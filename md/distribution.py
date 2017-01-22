
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


    def generate_distribution(self):
        # do anything
        return None


    def sample_distribution(self):
        # do anything
        return None