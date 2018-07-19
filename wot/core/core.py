# -*- coding: utf-8 -*-

from multiprocessing import Process
from wot.population import Population
import wot.core
import os
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
        If None, all files named `{prefix}_{t0}_{t1}.{extension}` will be
        considered as transport maps.
        The default prefix for cached transport maps is 'tmaps'
    """

    default_tmap_prefix = "tmaps"

    def __init__(self, matrix, transport_maps_directory, transport_maps_prefix = None, max_threads = None):
        self.matrix = matrix
        self.tmap_dir = transport_maps_directory
        self.tmap_prefix = transport_maps_prefix
        self.tmaps = wot.core.scan_transport_map_directory(self)
        self.timepoints = sorted(set(matrix.row_meta['day']))
        self.ot_config = {}
        if max_threads is None:
            self.max_threads = 8
        else:
            self.max_threads = max_threads

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

    def compute_all_transport_maps(self, day_pairs = None, force = False):
        """
        Computes all required transport maps and caches everything for future use.

        Parameters
        ----------
        day_pairs : list of (float, float), optional
            Day pairs to compute the transport maps for.
            If None, maps for all consecutive days will be computed.
        force : bool, optional, default : False
            Force recomputation of each transport map, after config update for instance.

        Returns
        -------
        None
            Only computes and caches all transport maps, does not return them.
        """
        t = self.timepoints
        if day_pairs is None:
            day_pairs = [ (t[i], t[i+1]) for i in range(len(t) - 1) ]

        if not force:
            day_pairs = [ x for x in day_pairs if self.tmaps.get(x, None) is None ]

        m = self.max_threads

        if m > 1 :
            procs = []
            for s, d in day_pairs:
                p = Process(target=self.compute_transport_map, args=(s,d,))
                procs.append(p)

            for i in range(len(procs) + m):
                if i >= m :
                    procs[i - m].join()
                if i < len(procs):
                    procs[i].start()
            self.tmaps = wot.core.scan_transport_map_directory(self)
        else:
            for s, d in day_pairs :
                self.compute_transport_map(s, d)

    def compute_transport_map(self, t0, t1):
        """
        Computes the transport map from time t0 to time t1

        Parameters
        ----------
        t0 : float
            Source timepoint for the transport map
        t1 : float
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
        path += "_{}_{}.loom".format(t0, t1)
        config = { **self.ot_config, 't0': t0, 't1': t1 }
        tmap = wot.ot.OptimalTransportHelper.compute_single_transport_map(self.matrix, config)
        wot.io.write_dataset(tmap, os.path.join(self.tmap_dir, path),
                output_format="loom", txt_full=False)
        self.tmaps[(t0, t1)] = path

    def transport_map(self, t0, t1):
        """
        Loads a transport map for a given pair of timepoints.

        Parameters
        ----------
        t0 : int or float
            Source timepoint of the transport map.
        t1 : int of float
            Destination timepoint of the transport map.

        Returns
        -------
        tmap : wot.Dataset
            The transport map from t0 to t1
        """
        return wot.core.load_transport_map(self, t0, t1)

    def can_push_forward(self, *populations):
        """
        Checks if the populations can be pushed forward.

        Parameters
        ----------
        *populations : wot.Population
            Measure over the cells at a given timepoint to be pushed forward.

        Returns
        -------
        result : bool
            True if the populations can be pushed forward

        Raises
        ------
        ValueError
            If all populations are not in the same timepoint
        """
        return self.timepoints.index(wot.core.unique_timepoint(*populations)) \
                < len(self.timepoints) - 1

    def can_pull_back(self, *populations):
        """
        Checks if the populations can be pulled back.

        Parameters
        ----------
        *populations : wot.Population
            Measure over the cells at a given timepoint to be pulled back.

        Returns
        -------
        result : bool
            True if the populations can be pulled back.

        Raises
        ------
        ValueError
            If all populations are not in the same timepoint
        """
        return self.timepoints.index(wot.core.unique_timepoint(*populations)) > 0

    def push_forward(self, *populations, to_time = None, normalize=True):
        """
        Pushes the population forward through the computed transport maps

        Parameters
        ----------
        *populations : wot.Population
            Measure over the cells at a given timepoint to be pushed forward.
        to_time : int or float, optional
            Destination timepoint to push forward to.
        normalize : bool, optional, default: True
            Wether to normalize to a probability distribution or keep growth.

        Returns
        -------
        result : wot.Population
            The push forward of the input population through the proper transport map.
            Array of populations if several populations were given as input.

        Raises
        ------
        ValueError
            If there is no further timepoint to push the population forward.
        ValueError
            If several populations are given as input but dot live in the same timepoint.

        Examples
        --------
        >>> core.push_forward(pop, to_time = 2) # -> wot.Population
        Pushing several populations at once
        >>> core.push_forward(pop1, pop2, pop3) # -> list of wot.Population
        Pulling back after pushing forward
        >>> core.pull_back(core.push_forward(pop))
        Same, but several populations at once
        >>> core.pull_back(* core.push_forward(pop1, pop2, pop3))
        """
        i = self.timepoints.index(wot.core.unique_timepoint(*populations))
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
            if normalize:
                p = (p.T / np.sum(p, axis=1)).T
            i += 1

        result = [ Population(self.timepoints[i], p[k,:]) for k in range(p.shape[0]) ]
        if len(result) == 1:
            return result[0]
        else:
            return result

    def pull_back(self, *populations, to_time = None, normalize=True):
        """
        Pulls the population back through the computed transport maps

        Parameters
        ----------
        *populations : wot.Population
            Measure over the cells at a given timepoint to be pushed forward.
        to_time : int or float, optional
            Destination timepoint to pull back to.
        normalize : bool, optional, default: True
            Wether to normalize to a probability distribution or keep growth.

        Returns
        -------
        result : wot.Population
            The pull back of the input population through the proper transport map.
            Array of populations if several populations were given as input.

        Raises
        ------
        ValueError
            If there is no previous timepoint to pull the population back.
        ValueError
            If several populations are given as input but dot live in the same timepoint.

        Examples
        --------
        >>> core.pull_back(pop, to_time = 0) # -> wot.Population
        Pushing several populations at once
        >>> core.pull_back(pop1, pop2, pop3) # -> list of wot.Population
        Pulling back after pushing forward
        >>> core.pull_back(core.push_forward(pop))
        Same, but several populations at once
        >>> core.pull_back(* core.push_forward(pop1, pop2, pop3))
        """
        i = self.timepoints.index(wot.core.unique_timepoint(*populations))
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
            if normalize:
                p = (p.T / np.sum(p, axis=1)).T
            i -= 1

        result = [ Population(self.timepoints[i], p[k,:]) for k in range(p.shape[0]) ]
        if len(result) == 1:
            return result[0]
        else:
            return result

    def ancestors(self, *populations, at_time=None):
        """
        Computes the ancestors of a given population by pulling back through transport maps

        Parameters
        ----------
        *populations : wot.Population
            Measure over the cells at a given timepoint to compute ancestors for.
        at_time : int or float, optional
            Timepoint for which to compute the ancestors.
            If None, compute ancestors for the previous available time point.

        Returns
        -------
        ancestors : wot.Population or list of wot.Population
            A population of cells, at the destination timepoint, most likely to be the ancestors of the input population.
            List if several populations were given, single population otherwise.

        Raises
        ------
        ValueError
            If the selected destination timepoint does not exist.
        ValueError
            If the selected destination is after the original timepoint.

        Examples
        --------
        >>> core.ancestors(pop, at_time = 0) # -> wot.Population
        # Using several populations at once
        >>> core.ancestors(pop1, pop2, pop3) # -> list of wot.Population
        # Chaining ancestors and descendants
        >>> core.ancestors(core.descendants(pop))
        # Same, but several populations at once
        >>> core.ancestors(* core.descendants(pop1, pop2, pop3))

        Notes
        -----
        If population.time is 7 and at_time is 5, the Core would pull back through two transport maps.
        This method is only and alias to Core.pull_back
        """
        return self.pull_back(*populations, to_time = at_time)

    def descendants(self, *populations, at_time=None):
        """
        Computes the descendants of a given population by pushing forward through transport maps

        Parameters
        ----------
        *populations : wot.Population
            Measure over the cells at a given timepoint to compute ancestors for.
        at_time : int or float, optional
            Timepoint for which to compute the ancestors.
            If None, compute ancestors for the previous available time point.

        Returns
        -------
        descendants : wot.Population or list of wot.Population
            A population of cells at the destination timepoint, most likely to be the descendants of the input population.
            List if several populations were given, single population otherwise.

        Raises
        ------
        ValueError
            If the selected destination timepoint does not exist.
        ValueError
            If the selected destination is before the original timepoint.

        Examples
        --------
        >>> core.descendants(pop, at_time = 2) # -> wot.Population
        # Using several populations at once
        >>> core.descendants(pop1, pop2, pop3) # -> list of wot.Population
        # Chaining ancestors and descendants
        >>> core.ancestors(core.descendants(pop))
        # Same, but several populations at once
        >>> core.ancestors(* core.descendants(pop1, pop2, pop3))

        Notes
        -----
        If population.time is 5 and at_time is 7, the Core would push forward through two transport maps.
        This method is only and alias to Core.push_forward
        """
        return self.push_forward(*populations, to_time = at_time)

    def population_from_ids(self, *ids, at_time=None):
        """
        Constructs a population uniformly distributed among the ids given as input.

        Parameters
        ----------
        *ids : list of str
            The list of cell ids that belong to that population.
        at_time : int or float, optional
            The time at which to construct the population.
            Cells that come from a different time point will be ignored.

        Returns
        -------
        *populations : wot.Population
            A population, uniformly distributed over the cells given as input.
            List if several lists of ids were given, single population otherwise.

        Raises
        ------
        ValueError
            If at_time is not specified and all cells do not live in the same timepoint.
        ValueError
            If the generated population would be empty.

        Examples
        --------
        >>> cell_set = [ 'cell_1', 'cell_2', 'cell_3' ]
        >>> core.population_from_ids(cell_set) # -> wot.Population
        Multiple populations at once
        >>> multi_cell_sets = {
        >>>   'set_a': [ 'cell_a1', 'cell_a2'],
        >>>   'set_b': [ 'cell_b1', 'cell_b2'],
        >>> }
        >>> core.population_from_ids(* multi_cell_sets.values()) # -> list of wot.Population

        Notes
        -----
        The Population class is a measure over the cells at a given timepoint.
        It does not necessarily sum to 1. However, this method always returns a probability distribution over the cells of that time point.
        """
        day = at_time
        all_ids = [ i for ids_el in ids for i in ids_el ]
        cell_inds = self.matrix.row_meta.index.get_indexer_for(all_ids)

        if at_time is None:
            day = self.matrix.row_meta.loc[ids[0][0], 'day']
            if not all(self.matrix.row_meta.iloc[cell_inds]['day'] == day):
                raise ValueError("All cells do not live in the same timepoint. Please choose one")

        day_query = self.matrix.row_meta['day'] == day
        all_inds = np.where(day_query)[0]

        def get_population(ids_el):
            cell_inds = self.matrix.row_meta.index.get_indexer_for(ids_el)
            p = [ 1 if id in cell_inds else 0 for id in all_inds ]
            p = np.asarray(p, dtype=np.float64)
            return Population(day, p / np.sum(p))

        result = [ get_population(ids_el) for ids_el in ids ]
        if len(result) == 1:
            return result[0]
        else:
            return result

    def cell_ids(self, population):
        day = population.time
        return list(self.matrix.row_meta.index[self.matrix.row_meta['day'] == day])

    def population_census(self, cell_set_matrix, *populations):
        """
        Get a census for a population with respect to a given cell set matrix

        Parameters
        ----------
        cell_set_matrix : wot.Dataset
            Dataset of 0s and 1s denoting membership in each cell set.
            Cells as rows, cell sets as columns.
        *populations : wot.Population or list of wot.Population
            The population to be considered

        Returns
        -------
        census : 1D-array or list of 1D-array
            The census for the population.
            census[i] is the probabiliy that a cell from that population belongs to cell set number i from the cell_set_matrix.
            List of censuses if a several populations were given as input, single census otherwise.

        Notes
        -----
        If several populations are given, they must all live in the same timepoint.
        """
        day = wot.core.unique_timepoint(*populations)
        all_ids_at_t = self.matrix.row_meta.index[self.matrix.row_meta['day'] == day]
        inter_ids = cell_set_matrix.row_meta.index.intersection(all_ids_at_t)
        if len(inter_ids) == 0:
            census = [ [0] * cell_set_matrix.x.shape[1] ] * len(populations)
        else:
            pop_indexer = all_ids_at_t.get_indexer_for(inter_ids)
            csm_indexer = cell_set_matrix.row_meta.index.get_indexer_for(inter_ids)
            def get_census(p):
                return np.dot(p[pop_indexer], cell_set_matrix.x[csm_indexer,:])
            norm = lambda p : p if np.isclose(np.sum(p), 0) else p / np.sum(p)
            census = np.asarray([get_census(norm(pop.p)) for pop in populations],
                    dtype=np.float64)

        if len(census) == 1:
            return census[0]
        else:
            return census
