import wot.model

from wot.population import Population
import wot.model
import wot.io
import os
import glob
import numpy as np
import pandas as pd
import h5py


class TransportMapModel:
    def __init__(self, tmaps):
        self.day_pairs = set()
        unique_timepoints = set()
        self.timepoint_to_ids = {}
        for key in tmaps:
            t0 = key[0]
            t1 = key[1]
            self.day_pairs.add((t0, t1))
            unique_timepoints.add(t0)
            unique_timepoints.add(t1)
            load_t0 = self.timepoint_to_ids.get(t0) is None
            load_t1 = self.timepoint_to_ids.get(t1) is None
            if load_t0 or load_t1:
                path = tmaps[key]
                f = h5py.File(path, 'r')
                if load_t0:
                    ids = f['/row_attrs/id'][()].astype(str)
                    df = pd.DataFrame(index=ids, data={'day': t0})
                    self.timepoint_to_ids[t0] = df
                if load_t1:
                    ids = f['/col_attrs/id'][()].astype(str)
                    df = pd.DataFrame(index=ids, data={'day': t1})
                    self.timepoint_to_ids[t1] = df
                f.close()
        self.meta = None
        self.timepoints = sorted(unique_timepoints)
        for t in self.timepoints:
            df = self.timepoint_to_ids[t]
            if self.meta is None:
                self.meta = df
            else:
                self.meta = pd.concat((self.meta, df), copy=False)
        self.tmaps = tmaps

    @staticmethod
    def from_directory(tmap_out, with_covariates=False):
        """
        Creates a wot.TransportMapModel from an output directory
        :param tmap_out:
        :param with_covariates:
        :return:
        """
        tmap_dir, tmap_prefix = os.path.split(tmap_out)
        tmap_dir = tmap_dir or '.'
        tmap_prefix = tmap_prefix or "tmaps"
        tmaps = {}
        pattern = tmap_prefix
        if with_covariates:
            pattern += '_[0-9]*.[0-9]*_[0-9]*.[0-9]*_cv[0-9]*_cv[0-9]*.*loom'
        else:
            pattern += '_[0-9]*.[0-9]*_[0-9]*.[0-9]*.*loom'

        files = glob.glob(os.path.join(tmap_dir, pattern))

        for path in files:
            if not os.path.isfile(path):
                continue
            basename, ext = wot.io.get_filename_and_extension(path)
            tokens = basename.split('_')
            try:
                if with_covariates:
                    t1 = float(tokens[-4])
                    t2 = float(tokens[-3])
                    cv1 = int(tokens[-2][2:])
                    cv2 = int(tokens[-1][2:])
                    tmaps[(t1, t2, cv1, cv2)] = path
                else:
                    t1 = float(tokens[-2])
                    t2 = float(tokens[-1])
                    tmaps[(t1, t2)] = path
            except ValueError:
                continue
        if len(tmaps) is 0:
            raise ValueError('No transport maps found')
        return TransportMapModel(tmaps=tmaps)

    def get_transport_map(self, t0, t1, covariate=None):
        """
        Loads a transport map for a given pair of timepoints.

        Parameters
        ----------
        t0 : int or float
            Source timepoint of the transport map.
        t1 : int of float
            Destination timepoint of the transport map.
        covariate : None or (int, int), optional
            Restrict to certain covariate values. Do not restrict if None

        Returns
        -------
        tmap : wot.Dataset
            The transport map from t0 to t1
        """
        if t0 not in self.timepoints or t1 not in self.timepoints:
            raise ValueError("Timepoints {}, {} not found".format(t0, t1))

        atomic = (self.day_pairs is not None and (t0, t1) in self.day_pairs) \
                 or self.timepoints.index(t1) == self.timepoints.index(t0) + 1

        if not atomic and covariate is not None:
            raise ValueError("Covariate-restricted transport maps can only be atomic")

        if atomic:
            if covariate is None:
                return wot.io.read_dataset(self.tmaps.get((t0, t1)))
            else:
                cv0, cv1 = covariate
                return wot.io.read_dataset(self.tmaps.get((t0, t1, cv0, cv1)))
        else:
            path = wot.model.find_path(t0, t1, self.day_pairs, self.timepoints)
            return wot.model.chain_transport_maps(self, path)

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
        return self.timepoints.index(wot.model.unique_timepoint(*populations)) \
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
        return self.timepoints.index(wot.model.unique_timepoint(*populations)) > 0

    def push_forward(self, *populations, to_time=None, normalize=True, as_list=False):
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
        as_list : bool, optional, default: False
            Wether to return a list of length 1 when a single element is passed, or a Population

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
        >>> tmap_model.push_forward(pop, to_time = 2) # -> wot.Population
        Pushing several populations at once
        >>> tmap_model.push_forward(pop1, pop2, pop3) # -> list of wot.Population
        Pulling back after pushing forward
        >>> tmap_model.pull_back(tmap_model.push_forward(pop))
        Same, but several populations at once
        >>> tmap_model.pull_back(* tmap_model.push_forward(pop1, pop2, pop3))
        """
        i = self.timepoints.index(wot.model.unique_timepoint(*populations))
        j = i + 1 if to_time is None else self.timepoints.index(to_time)

        if i == -1:
            raise ValueError("Timepoint not found")
        if j == -1:
            raise ValueError("Destination timepoint not found")
        if j >= len(self.timepoints):
            raise ValueError("No further timepoints. Unable to push forward")
        if i > j:
            raise ValueError("Destination timepoint is before source. Unable to push forward")

        p = np.vstack([pop.p for pop in populations])
        while i < j:
            t0 = self.timepoints[i]
            t1 = self.timepoints[i + 1]
            tmap = self.get_transport_map(t0, t1)
            p = np.dot(p, tmap.x)
            if normalize:
                p = (p.T / np.sum(p, axis=1)).T
            i += 1

        result = [Population(self.timepoints[i], p[k, :]) for k in range(p.shape[0])]
        if len(result) == 1 and not as_list:
            return result[0]
        else:
            return result

    def pull_back(self, *populations, to_time=None, normalize=True, as_list=False):
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
        as_list : bool, optional, default: False
            Wether to return a listof length 1 when a single element is passed, or a Population

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
        >>> tmap_model.pull_back(pop, to_time = 0) # -> wot.Population
        Pushing several populations at once
        >>> tmap_model.pull_back(pop1, pop2, pop3) # -> list of wot.Population
        Pulling back after pushing forward
        >>> tmap_model.pull_back(tmap_model.push_forward(pop))
        Same, but several populations at once
        >>> tmap_model.pull_back(* tmap_model.push_forward(pop1, pop2, pop3))
        """
        i = self.timepoints.index(wot.model.unique_timepoint(*populations))
        j = i - 1 if to_time is None else self.timepoints.index(to_time)

        if i == -1:
            raise ValueError("Timepoint not found")
        if i == 0:
            raise ValueError("No previous timepoints. Unable to pull back")
        if j == -1:
            raise ValueError("Destination timepoint not found")
        if i < j:
            raise ValueError("Destination timepoint is after source. Unable to pull back")

        p = np.vstack([pop.p for pop in populations])
        while i > j:
            t1 = self.timepoints[i]
            t0 = self.timepoints[i - 1]
            tmap = self.get_transport_map(t0, t1)
            p = np.dot(tmap.x, p.T).T
            if normalize:
                p = (p.T / np.sum(p, axis=1)).T
            i -= 1

        result = [Population(self.timepoints[i], p[k, :]) for k in range(p.shape[0])]
        if len(result) == 1 and not as_list:
            return result[0]
        else:
            return result

    def ancestors(self, *populations, at_time=None, as_list=False):
        """
        Computes the ancestors of a given population by pulling back through transport maps

        Parameters
        ----------
        *populations : wot.Population
            Measure over the cells at a given timepoint to compute ancestors for.
        at_time : int or float, optional
            Timepoint for which to compute the ancestors.
            If None, compute ancestors for the previous available time point.
        as_list : bool, optional, default: False
            Wether to return a listof length 1 when a single element is passed, or a Population

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
        >>> tmap_model.ancestors(pop, at_time = 0) # -> wot.Population
        # Using several populations at once
        >>> tmap_model.ancestors(pop1, pop2, pop3) # -> list of wot.Population
        # Chaining ancestors and descendants
        >>> tmap_model.ancestors(tmap_model.descendants(pop))
        # Same, but several populations at once
        >>> tmap_model.ancestors(* tmap_model.descendants(pop1, pop2, pop3))

        Notes
        -----
        If population.time is 7 and at_time is 5, the OTModel would pull back through two transport maps.
        This method is only and alias to OTModel.pull_back
        """
        return self.pull_back(*populations, to_time=at_time, as_list=as_list)

    def descendants(self, *populations, at_time=None, as_list=False):
        """
        Computes the descendants of a given population by pushing forward through transport maps

        Parameters
        ----------
        *populations : wot.Population
            Measure over the cells at a given timepoint to compute ancestors for.
        at_time : int or float, optional
            Timepoint for which to compute the ancestors.
            If None, compute ancestors for the previous available time point.
        as_list : bool, optional, default: False
            Wether to return a listof length 1 when a single element is passed, or a Population

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
        >>> tmap_model.descendants(pop, at_time = 2) # -> wot.Population
        # Using several populations at once
        >>> tmap_model.descendants(pop1, pop2, pop3) # -> list of wot.Population
        # Chaining ancestors and descendants
        >>> tmap_model.ancestors(tmap_model.descendants(pop))
        # Same, but several populations at once
        >>> tmap_model.ancestors(* tmap_model.descendants(pop1, pop2, pop3))

        Notes
        -----
        If population.time is 5 and at_time is 7, the OTModel would push forward through two transport maps.
        This method is only and alias to OTModel.push_forward
        """
        return self.push_forward(*populations, to_time=at_time, as_list=as_list)

    def population_from_ids(self, *ids, at_time):
        """
        Constructs a population uniformly distributed among the ids given as input.

        Parameters
        ----------
        *ids : list of str
            The list of cell ids that belong to that population.
        at_time : int or float
            The time at which to construct the population.
            Cells that come from a different time point will be ignored.

        Returns
        -------
        *populations : list of wot.Population
            A population, uniformly distributed over the cells given as input.
            Returns None if the generated population would be empty



        Examples
        --------
        >>> cell_set = [ 'cell_1', 'cell_2', 'cell_3' ]
        >>> tmap_model.population_from_ids(cell_set) # -> wot.Population
        Multiple populations at once
        >>> multi_cell_sets = {
        >>>   'set_a': [ 'cell_a1', 'cell_a2'],
        >>>   'set_b': [ 'cell_b1', 'cell_b2'],
        >>> }
        >>> tmap_model.population_from_ids(* multi_cell_sets.values()) # -> list of wot.Population

        Notes
        -----
        The Population class is a measure over the cells at a given timepoint.
        It does not necessarily sum to 1. However, this method always returns a probability distribution over the cells of that time point.
        """

        day = float(at_time)
        df = self.timepoint_to_ids[day]

        def get_population(ids_el):
            cell_inds = df.index.get_indexer_for(ids_el)
            cell_inds = cell_inds[cell_inds > -1]
            if len(cell_inds) is 0:
                return None
            p = np.zeros(len(df), dtype=np.float64)
            p[cell_inds] = 1.0
            return Population(day, p / np.sum(p))

        result = [get_population(ids_el) for ids_el in ids]
        return result

    def population_from_cell_sets(self, cell_sets, at_time):
        """
        Similar to population_from_ids for cell sets

        Parameters
        ----------
        cell_sets : dict of str: list of str
            The dictionnary of ids
        at_time : float, optional
            The timepoint to consider

        Returns
        -------
        populations : dict of str: wot.Population
            The resulting populations
        """
        keys = list(cell_sets.keys())
        populations = self.population_from_ids(*[cell_sets[name] for name in keys], at_time=at_time)
        return {keys[i]: populations[i] for i in range(len(keys)) if populations[i] is not None}

    def cell_ids(self, population):
        day = population.time
        return self.timepoint_to_ids[day].index.values

    def compute_ancestor_census(self, cset_matrix, *populations):
        """
        Computes the census for the populations (for both ancestors and descendants).

        Parameters
        ----------
        tmap_model : wot.TransportMapModel
            The OTModel used to find ancestors and descendants of the population
        cset_matrix : wot.Dataset
            The cell set matrix, cells as rows, cell sets as columns. 1s denote membership.
        *populations : wot.Population
            The target populations
        """
        initial_populations = populations
        timepoints = []
        census = []

        def update(head, populations):
            x = 0 if head else len(census)
            timepoints.insert(x, wot.model.unique_timepoint(*populations))
            census.insert(x, self.population_census(cset_matrix, *populations))

        update(True, populations)
        while self.can_pull_back(*populations):
            populations = self.pull_back(*populations, as_list=True)
            update(True, populations)
        populations = initial_populations
        while self.can_push_forward(*populations):
            populations = self.push_forward(*populations, as_list=True)
            update(False, populations)
        census = np.asarray(census)
        if census.ndim == 3:
            # rearrange dimensions when more than one population is passed
            census = np.asarray([census[:, i, :] for i in range(census.shape[1])])
        return timepoints, census

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
            List of censuses

        Notes
        -----
        If several populations are given, they must all live in the same timepoint.
        """
        day = wot.model.unique_timepoint(*populations)
        df = self.timepoint_to_ids[day]
        inter_ids = cell_set_matrix.row_meta.index.intersection(df.index)
        if len(inter_ids) == 0:
            census = [[0] * cell_set_matrix.x.shape[1]] * len(populations)
        else:
            pop_indexer = df.index.get_indexer_for(inter_ids)
            csm_indexer = cell_set_matrix.row_meta.index.get_indexer_for(inter_ids)

            def get_census(p):
                return np.dot(p[pop_indexer], cell_set_matrix.x[csm_indexer, :])

            norm = lambda p: p if np.isclose(np.sum(p), 0) else p / np.sum(p)
            census = np.asarray([get_census(norm(pop.p)) for pop in populations],
                                dtype=np.float64)

        return census


