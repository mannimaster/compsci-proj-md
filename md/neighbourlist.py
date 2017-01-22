import numpy as np

class neighbourlist(object):

    def __init__(self):
        #raise NotImplementedError('You cannot create a neighbourlist object, just use the classmethods')
        return


#    def compute_neighbourlist(self):
        # computeneighbourlist
#        return None

    def minDistance(self, position1, position2, box_length, half_box_length):
        """
        minDistance(position1,position2, box_length, half_box_length)

        minimal distancebetween two particles

        Parameters
        ----------
        position1, position2 : a skalar or np.array
            ..representing a position

        box_length : scalar number, positiv
            side length of the box otherwise denoted as "L"

        half_box_length : scalar number
            half the box_length

        Returns
        -------
        d : scalar number
            minimal distance between the inputed position considering periodic boundary conditions
        dr : vector distance
            minimal distance between the inputed position considering periodic boundary conditions
        """
        d = np.linalg.norm((position1 - position2 + half_box_length) % box_length - half_box_length)
        dr = ((position1 - position2 + half_box_length) % box_length - half_box_length)
        return d, dr

    def neighborList(self, n_particles, positions, r_cutoff, r_s, box_length):
        """
        neighborList(n_particles, positions, r_cutoff, r_s, box_length)

        creates a neighbor list

        Parameters
        ----------
        n_particles : scalar number
            number of particles

        positions : np.array or single numbers
            array of the positionvectors with each line a position

        r_cutoff : scalar number
            cutoff radius

        r_s : scalar number
            skin radius of the maximal possible movement in one time step

        box_length : scalar number, positiv
            side length of the box otherwise denoted as "L"

        Returns
        -------
        neighbor_list : dictionary of lists
            all neighbors within given cutoff radius + skin radius
        """
        neighbor_list = {}
        half_box_length = box_length / 2

        for i in range(0, n_particles):
            if i not in neighbor_list:
                # empty list of neighbors
                neighbor_list[i] = []
            for j in range(i + 1, n_particles):
                d, dr = minDistance(positions[i], positions[j], box_length, half_box_length)
                if d < r_cutoff + r_s:
                    # neighborship i => j
                    neighbor_list[i].append(j)
                    # neighborship j => i
                    if j not in neighbor_list:
                        neighbor_list[j] = [i]
                    else:
                        neighbor_list[j].append(i)
        return neighbor_list