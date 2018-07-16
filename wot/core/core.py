# -*- coding: utf-8 -*-

from wot.population import Population
import wot.core
import numpy as np
import pandas as pd

class Core:
    """
    The Core takes care of computing and properly caching the transport maps
    needed when necessary. All computations using transport maps should be
    performed through the Core.

    Parameters
    ----------
    matrix : wot.Dataset
        The gene expression matrix for this core
    transport_maps_directory : str
        Path to the transport map directory, where already computed transport
        maps are stored and future transport maps may be cached.
    transport_maps_prefix : str, optional
        Prefix to use for the transport maps. This can highly speed up
        initialization if the directory is filled with other non-tmap files,
        and allows to have several transport maps not overriding each other.
        If None, all files named `{prefix}_{t1}_{t2}.{extension}` will be
        considered as transport maps.
        The default prefix for cached transport maps is 'tmaps'
    """

    default_tmap_prefix = "tmaps"

    def __init__(self, matrix, transport_maps_directory, transport_maps_prefix = None):
        self.matrix = matrix
        self.tmap_dir = transport_maps_directory
        self.tmap_prefix = transport_maps_prefix
        self.tmaps = wot.core.scan_transport_map_directory(self)
        self.timepoints = sorted(set(matrix.row_meta['day']))
        self.ot_config = {}

    def set_ot_config(self, **kwargs):
        """
        Set parameters for the Optimal Transport computation

        Parameters
        ----------
        **kwargs : dict
            Dictionnary of parameters. Will be inserted as is into OT configuration.

        Example
        -------
        core.set_ot_config(epsilon = .01)
        core.set_ot_config(lambda1 = 50, lambda2 = 80)
        """
        for k in kwargs.keys():
            self.ot_config[k] = kwargs[k]

    def compute_all_transport_maps(self, force = False):
        """
        Computes all required transport maps and caches everything for future use.

        Parameters
        ----------
        force : bool, optional, default : False
            Force recomputation of each transport map, after config update for instance.

        Returns
        -------
        None
            Only computes and caches all transport maps, does not return them.
        """
        t = self.timepoints
        for i in range(len(t) - 1):
            if force or self.tmaps.get((t[i], t[i+1]), None) is None:
                self.compute_transport_map(t[i], t[i+1])

    def compute_transport_map(self, t1, t2):
        """
        Computes the transport map from time t1 to time t2

        Parameters
        ----------
        t1 : float
            Source timepoint for the transport map
        t2 : float
            Destination timepoint for the transport map

        Returns
        -------
        None
            Only computes and caches the transport maps, does not return it.
        """
        if self.tmap_prefix is None:
            path = self.default_tmap_prefix
        else:
            path = self.tmap_prefix
        path += "_{}_{}.loom".format(t1, t2)
        config = { 't0': t1, 't1': t2 }
        tmap = wot.ot.OptimalTransportHelper.compute_single_transport_map(self.matrix, config)
        wot.io.write_dataset(tmap, path, output_format="loom", txt_full=False)
        self.tmaps[(t1, t2)] = path

    def transport_map(self, t1, t2):
        """
        Loads a transport map for a given pair of timepoints.

        Parameters
        ----------
        t1 : int or float
            Source timepoint of the transport map.
        t2 : int of float
            Destination timepoint of the transport map.

        Returns
        -------
        tmap : wot.Dataset
            The transport map from t1 to t2
        """
        return wot.core.load_transport_map(self, t1, t2)

    def push_forward(self, *populations, to_time = None):
        """
        Pushes the population forward through the computed transport maps

        Parameters
        ----------
        *populations : wot.Population
            Measure over the cells at a given timepoint to be pushed forward.
        to_time : int or float, optional
            Destination timepoint to push forward to.

        Returns
        -------
        result : wot.Population
            The push forward of the input population through the proper transport map.
            Array of populations if several populations were given as input.

        Raises
        ------
        ValueError
            If there is no further timepoint to push the population forward.
            If several populations are given as input but dot live in the same timepoint.

        Examples
        --------
        # Basic example
        core.push_forward(pop, to_time = 0)
        # Pushing several populations at once
        core.push_forward(pop1, pop2, pop3)
        # Pulling back after pushing forward
        core.pull_back(core.push_forward(pop))
        # Same, but several populations at once
        core.pull_back(* core.push_forward(pop1, pop2, pop3))
        """
        times = set([ pop.time for pop in populations ])
        if len(times) > 1:
            raise ValueError("Several populations were given, but they do not all live in the same timepoint")
        i = self.timepoints.index(list(times)[0])
        j = i + 1 if to_time is None else self.timepoints.index(to_time)

        if i == -1:
            raise ValueError("Timepoint not found")
        if j == -1:
            raise ValueError("Destination timepoint not found")
        if j >= len(self.timepoints):
            raise ValueError("No further timepoints. Unable to push forward")
        if i > j :
            raise ValueError("Destination timepoint is before source. Unable to push forward")

        p = np.vstack([ pop.p for pop in populations ])
        while i < j:
            t0 = self.timepoints[i]
            t1 = self.timepoints[i+1]
            tmap = self.transport_map(t0, t1)
            p = np.dot(p, tmap.x)
            i += 1

        result = [ Population(self.timepoints[i], p[k,:]) for k in range(p.shape[0]) ]
        if len(result) == 1:
            return result[0]
        else:
            return result

    def pull_back(self, *populations, to_time = None):
        """
        Pulls the population back through the computed transport maps

        Parameters
        ----------
        *populations : wot.Population
            Measure over the cells at a given timepoint to be pushed forward.
        to_time : int or float, optional
            Destination timepoint to pull back to.

        Returns
        -------
        result : wot.Population
            The pull back of the input population through the proper transport map.
            Array of populations if several populations were given as input.

        Raises
        ------
        ValueError
            If there is no previous timepoint to pull the population back.
            If several populations are given as input but dot live in the same timepoint.

        Examples
        --------
        # Basic example
        core.pull_back(pop, to_time = 0)
        # Pushing several populations at once
        core.pull_back(pop1, pop2, pop3)
        # Pulling back after pushing forward
        core.pull_back(core.push_forward(pop))
        # Same, but several populations at once
        core.pull_back(* core.push_forward(pop1, pop2, pop3))
        """
        times = set([ pop.time for pop in populations ])
        if len(times) > 1:
            raise ValueError("Several populations were given, but they do not all live in the same timepoint")
        i = self.timepoints.index(list(times)[0])
        j = i - 1 if to_time is None else self.timepoints.index(to_time)

        if i == -1:
            raise ValueError("Timepoint not found")
        if i == 0:
            raise ValueError("No previous timepoints. Unable to pull back")
        if j == -1:
            raise ValueError("Destination timepoint not found")
        if i < j :
            raise ValueError("Destination timepoint is after source. Unable to pull back")

        p = np.vstack([ pop.p for pop in populations ])
        while i > j:
            t1 = self.timepoints[i]
            t0 = self.timepoints[i-1]
            tmap = self.transport_map(t0, t1)
            p = np.dot(tmap.x, p.T).T
            i -= 1

        result = [ Population(self.timepoints[i], p[k,:]) for k in range(p.shape[0]) ]
        if len(result) == 1:
            return result[0]
        else:
            return result

    def ancestors(self, population, at_time=None):
        """
        Computes the ancestors of a given population by pulling back through transport maps

        Parameters
        ----------
        population : wot.Population
            Measure over the cells at a given timepoint to compute ancestors for.
        at_time : int or float, optional
            Timepoint for which to compute the ancestors.
            If None, compute ancestors for the previous available time point.

        Returns
        -------
        ancestors : wot.Population
            A population of cells, at the destination timepoint, most likely to
             be the ancestors of the input population.

        Raises
        ------
        ValueError
            If the selected destination timepoint does not exist.
            If the selected destination is after the original timepoint.

        Notes
        -----
        If population.time is 7 and at_time is 5, the Core would pull back through two transport maps.
        """
        return self.pull_back(population, to_time = at_time)

    def descendants(self, population, at_time=None):
        """
        Computes the descendants of a given population by pushing forward through transport maps

        Parameters
        ----------
        population : wot.Population
            Measure over the cells at a given timepoint to compute ancestors for.
        at_time : int or float, optional
            Timepoint for which to compute the ancestors.
            If None, compute ancestors for the previous available time point.

        Returns
        -------
        descendants : wot.Population
            A population of cells at the destination timepoint, most likely to
             be the descendants of the input population.

        Raises
        ------
        ValueError
            If the selected destination timepoint does not exist.
            If the selected destination is before the original timepoint.

        Notes
        -----
        If population.time is 5 and at_time is 7, the Core would push forward through two transport maps.
        """
        return self.push_forward(population, to_time = at_time)

    def population_from_ids(self, ids, at_time=None):
        """
        Constructs a population uniformly distributed among the ids given as input.

        Parameters
        ----------
        ids : list of str
            The list of cell ids that belong to that population.
        at_time : int or float, optional
            The time at which to construct the population.
            Cells that come from a different time point will be ignored.

        Returns
        -------
        population : wot.Population
            A population, uniformly distributed over the cells given as input.

        Raises
        ------
        ValueError
            If at_time is not specified and all cells do not live in the same timepoint.
            If the generated population would be empty.

        Notes
        -----
        The Population class is a measure over the cells at a given timepoint.
          It does not necessarily sum to 1. However, this method always returns
          a probability distribution over the cells of that time point.
        """
        day = at_time
        cell_inds = self.matrix.row_meta.index.get_indexer_for(ids)

        if at_time is None:
            day = self.matrix.row_meta.loc[ids[0], 'day']
            if not all(self.matrix.row_meta.iloc[cell_inds]['day'] == day):
                raise ValueError("All cells do not live in the same timepoint. Please choose one")

        day_query = self.matrix.row_meta['day'] == day
        all_inds = np.where(day_query)[0]

        p = [ 1 if id in cell_inds else 0 for id in all_inds ]
        p = np.asarray(p, dtype=np.float64)
        return Population(day, p / np.sum(p))
