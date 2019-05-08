# -*- coding: utf-8 -*-

import numpy as np


class Population:
    """
    A Population is a measure over the the cells at given timepoint.

    Parameters
    ----------
    time : int or float
        The time at which the cells where measured.
    p : 1-D array-like
        Measure over the cells at the given timepoint.
    name : str
        Optional population name.
    """

    def __init__(self, time, p, name=None):
        self.time = time
        self.p = np.asarray(p, dtype=np.float64)
        self.name = name

    def normalize(self):
        """
        Make the measure sum to 1, i.e. be a probability distribution over cells.
        """
        self.p = self.p / self.p.sum()

    def make_binary(self):
        """
        Set non-zero values to 1.
        """
        p = np.zeros(len(self.p))
        p[self.p > 0.0] = 1.0
        self.p = p

    @staticmethod
    def get_missing_population(*populations):
        initial_p_sum = np.array([pop.p for pop in populations]).T.sum(axis=1)
        missing_cells = np.where(initial_p_sum == 0)[0]
        if len(missing_cells) > 0:
            missing_cells_p = np.zeros_like(populations[0].p)
            missing_cells_p[missing_cells] = 1.0
            return Population(populations[0].time, missing_cells_p, 'Other')
        return None

    @staticmethod
    def copy(*populations, normalize=None, add_missing=False):
        populations_copy = []
        for pop in populations:
            pop_copy = Population(pop.time, pop.p, pop.name)
            populations_copy.append(pop_copy)
        populations = populations_copy
        if add_missing:
            # add "other" population if any cells are missing across all populations
            missing_pop = Population.get_missing_population(*populations)
            if missing_pop is not None:
                populations.append(missing_pop)
        if normalize or normalize == False:
            for pop in populations:
                if normalize:
                    pop.normalize()
                else:
                    pop.make_binary()
        return populations
