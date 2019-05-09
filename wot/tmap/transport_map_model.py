import os

import anndata
import h5py
import numpy as np
import pandas as pd

import wot.io
import wot.tmap
from wot.population import Population


class TransportMapModel:
    """
       Creates a transport map model for operating on pre-computed transport maps

       Parameters
       ----------
       tmaps : dict
           Maps day pairs to transport map path.
       meta : pandas.DataFrame
           Cell metadata with index cell ids and 'day'. Must match the order in the transport maps.
       timepoints : list
           Sorted list of cell timepoints
        day_pairs : list
            List of (t1,t2)
       """

    def __init__(self, tmaps, meta, timepoints=None, day_pairs=None, cache=False):
        self.tmaps = tmaps
        self.meta = meta
        self.cache = cache
        if timepoints is None:
            timepoints = sorted(meta['day'].unique())
        self.timepoints = timepoints

        if day_pairs is None:
            day_pairs = [(timepoints[i], timepoints[i + 1]) for i in range(len(timepoints) - 1)]
        self.day_pairs = day_pairs

    def fates(self, populations):
        """
        Computes fates for each population

        Parameters
        ----------
        self : wot.TransportMapModel
            The TransportMapModel used to find fates
        populations : list of wot.Population
            The target populations such as ones from self.population_from_cell_sets. The populations must be from the same time.
        Returns
        -------
        fates : anndata.AnnData
            Rows : all cells, Columns : populations index. At point (i, j) : the probability that cell i belongs to population j
        """
        start_day = wot.tmap.unique_timepoint(*populations)  # check for unique timepoint
        populations = Population.copy(*populations, normalize=False, add_missing=True)
        pop_names = [pop.name for pop in populations]
        results = []

        results.insert(0, np.array([pop.p for pop in populations]).T)
        while self.can_pull_back(*populations):
            populations = self.pull_back(*populations, as_list=True, normalize=False)
            results.insert(0, np.array([pop.p for pop in populations]).T)

        X = np.concatenate(results)
        X /= X.sum(axis=1, keepdims=1)
        obs = self.meta.copy()
        obs = obs[obs['day'] <= start_day]
        return anndata.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pop_names))

    def transition_table(self, start_populations, end_populations):
        """
       Computes a transition table from the starting populations to the ending populations

       Parameters
       ----------
       self : wot.TransportMapModel
           The TransportMapModel
       start_populations : list of wot.Population
           The target populations such as ones from self.population_from_cell_sets. THe populations must be from the same time.

       Returns
       -------
       transition table : anndata.AnnData
           Rows : starting populations, Columns : ending populations.
       """
        # add "other" population if any cells are missing across all populations
        start_time = wot.tmap.unique_timepoint(*start_populations)
        start_populations = Population.copy(*start_populations, normalize=False, add_missing=True)
        end_populations = Population.copy(*end_populations, normalize=False, add_missing=True)
        wot.tmap.unique_timepoint(*end_populations)  # check for unique timepoint
        populations = end_populations
        results = []
        results.insert(0, np.array([pop.p for pop in populations]).T)
        while self.can_pull_back(*populations) and wot.tmap.unique_timepoint(*populations) > start_time:
            populations = self.pull_back(*populations, as_list=True, normalize=False)

        end_p = np.vstack([pop.p for pop in populations])
        start_p = np.vstack([pop.p for pop in start_populations])
        p = (start_p @ end_p.T)
        p = p / p.sum()
        return anndata.AnnData(X=p, obs=pd.DataFrame(index=[p.name for p in start_populations]),
                               var=pd.DataFrame(index=[p.name for p in end_populations]))

    def trajectories(self, populations):
        """
        Computes a trajectory for each population

        Parameters
        ----------
        self : wot.TransportMapModel
            The TransportMapModel used to find ancestors and descendants of the population
        populations : list of wot.Population
            The target populations such as ones from self.population_from_cell_sets. THe populations must be from the same time.

        Returns
        -------
        trajectories : anndata.AnnData
            Rows : all cells, Columns : populations index. At point (i, j) : the probability that cell i is an
            ancestor/descendant of population j
        """
        wot.tmap.unique_timepoint(*populations)  # check for unique timepoint
        trajectories = []
        populations = Population.copy(*populations, normalize=True, add_missing=False)
        initial_populations = populations

        def update(head, populations_to_update):
            idx = 0 if head else len(trajectories)
            trajectories.insert(idx, np.array([pop.p for pop in populations_to_update]).T)

        update(True, populations)
        while self.can_pull_back(*populations):
            populations = self.pull_back(*populations, as_list=True)
            update(True, populations)
        populations = initial_populations
        while self.can_push_forward(*populations):
            populations = self.push_forward(*populations, as_list=True)
            update(False, populations)

        return anndata.AnnData(X=np.concatenate(trajectories), obs=self.meta.copy(),
                               var=pd.DataFrame(index=[p.name for p in populations]))

    def get_coupling(self, t0, t1, covariate=None):
        """
        Loads a coupling for a given pair of timepoints.

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
        tmap : anndata.AnnData
            The transport map from t0 to t1
        """
        if t0 not in self.timepoints or t1 not in self.timepoints:
            raise ValueError("Day pair {}, {} not found".format(t0, t1))

        atomic = (self.day_pairs is not None and (t0, t1) in self.day_pairs) \
                 or self.timepoints.index(t1) == self.timepoints.index(t0) + 1

        if not atomic and covariate is not None:
            raise ValueError("Covariate-restricted transport maps can only be atomic")

        if atomic:
            if covariate is None:
                key = (t0, t1)
            else:
                cv0, cv1 = covariate
                key = (t0, t1, str(cv0), str(cv1))
            ds_or_path = self.tmaps.get(key)
            if ds_or_path is None:
                raise ValueError('No transport map found for {}', key)
            if type(ds_or_path) is anndata.AnnData:
                return ds_or_path
            ds = wot.io.read_dataset(ds_or_path)
            if self.cache:
                self.tmaps[key] = ds
            return ds

        else:
            path = wot.tmap.find_path(t0, t1, self.day_pairs, self.timepoints)
            return wot.tmap.chain_transport_maps(self, path)

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
        return self.timepoints.index(wot.tmap.unique_timepoint(*populations)) \
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
        return self.timepoints.index(wot.tmap.unique_timepoint(*populations)) > 0

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
        >>> self.push_forward(pop, to_time = 2) # -> wot.Population
        Pushing several populations at once
        >>> self.push_forward(pop1, pop2, pop3) # -> list of wot.Population
        Pulling back after pushing forward
        >>> self.pull_back(self.push_forward(pop))
        Same, but several populations at once
        >>> self.pull_back(* self.push_forward(pop1, pop2, pop3))
        """
        i = self.timepoints.index(wot.tmap.unique_timepoint(*populations))
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
            tmap = self.get_coupling(t0, t1)
            p = p @ tmap.X
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
        >>> self.pull_back(pop, to_time = 0) # -> wot.Population
        Pushing several populations at once
        >>> self.pull_back(pop1, pop2, pop3) # -> list of wot.Population
        Pulling back after pushing forward
        >>> self.pull_back(self.push_forward(pop))
        Same, but several populations at once
        >>> self.pull_back(* self.push_forward(pop1, pop2, pop3))
        """
        i = self.timepoints.index(wot.tmap.unique_timepoint(*populations))
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
            tmap = self.get_coupling(t0, t1)
            p = (tmap.X @ p.T).T
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
        >>> self.ancestors(pop, at_time = 0) # -> wot.Population
        # Using several populations at once
        >>> self.ancestors(pop1, pop2, pop3) # -> list of wot.Population
        # Chaining ancestors and descendants
        >>> self.ancestors(self.descendants(pop))
        # Same, but several populations at once
        >>> self.ancestors(* self.descendants(pop1, pop2, pop3))

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
        >>> self.descendants(pop, at_time = 2) # -> wot.Population
        # Using several populations at once
        >>> self.descendants(pop1, pop2, pop3) # -> list of wot.Population
        # Chaining ancestors and descendants
        >>> self.ancestors(self.descendants(pop))
        # Same, but several populations at once
        >>> self.ancestors(* self.descendants(pop1, pop2, pop3))

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
        >>> self.population_from_ids(cell_set) # -> wot.Population
        Multiple populations at once
        >>> multi_cell_sets = {
        >>>   'set_a': [ 'cell_a1', 'cell_a2'],
        >>>   'set_b': [ 'cell_b1', 'cell_b2'],
        >>> }
        >>> self.population_from_ids(* multi_cell_sets.values()) # -> list of wot.Population

        Notes
        -----
        The Population class is a measure over the cells at a given timepoint.
        It does not necessarily sum to 1. However, this method always returns a probability distribution over the cells of that time point.
        """

        day = float(at_time)
        df = self.meta[self.meta['day'] == day]

        def get_population(ids_el):
            cell_indices = df.index.get_indexer_for(ids_el)
            cell_indices = cell_indices[cell_indices > -1]

            if len(cell_indices) is 0:
                return None
            p = np.zeros(len(df), dtype=np.float64)
            p[cell_indices] = 1.0

            return Population(day, p)

        result = [get_population(ids_el) for ids_el in ids]
        return result

    def population_from_cell_sets(self, cell_sets, at_time):
        """
        Similar to population_from_ids for cell sets

        Parameters
        ----------
        cell_sets : dict of str: list of str
            The dictionary of ids
        at_time : float, optional
            The timepoint to consider

        Returns
        -------
        populations : list of wot.Population
            The resulting populations
        """
        keys = list(cell_sets.keys())
        populations = self.population_from_ids(*[cell_sets[name] for name in keys], at_time=at_time)
        filtered_populations = []
        for i in range(len(keys)):
            if populations[i] is not None:
                populations[i].name = keys[i]
                filtered_populations.append(populations[i])

        return filtered_populations

    def cell_ids(self, population):
        day = population.time
        df = self.meta[self.meta['day'] == day]
        return df.index.values

    def ancestor_census(self, cset_matrix, *populations):
        """
        Computes the census for the populations (for both ancestors and descendants).

        Parameters
        ----------
        self : wot.TransportMapModel
            The OTModel used to find ancestors and descendants of the population
        cset_matrix : anndata.AnnData
            The cell set matrix, cells as rows, cell sets as columns. 1s denote membership.
        *populations : wot.Population
            The target populations
        """
        initial_populations = populations
        timepoints = []
        census = []

        def update(head, populations):
            x = 0 if head else len(census)
            timepoints.insert(x, wot.tmap.unique_timepoint(*populations))
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
        cell_set_matrix : anndata.AnnData
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
        day = wot.tmap.unique_timepoint(*populations)
        df = self.meta[self.meta['day'] == day]
        inter_ids = cell_set_matrix.obs.index.intersection(df.index)
        if len(inter_ids) == 0:
            census = [[0] * cell_set_matrix.X.shape[1]] * len(populations)
        else:
            pop_indexer = df.index.get_indexer_for(inter_ids)
            csm_indexer = cell_set_matrix.obs.index.get_indexer_for(inter_ids)

            def get_census(p):
                return np.dot(p[pop_indexer], cell_set_matrix.X[csm_indexer, :])

            norm = lambda p: p if np.isclose(np.sum(p), 0) else p / np.sum(p)
            census = np.asarray([get_census(norm(pop.p)) for pop in populations],
                                dtype=np.float64)

        return census

    def to_json(self, path):
        import json
        meta = self.meta.to_dict(orient='list')
        meta['id'] = self.meta.index.values.tolist()
        paths = []
        day_pairs = list(self.day_pairs)
        for i in range(len(day_pairs)):
            paths.append(self.tmaps[day_pairs[i]])
        d = {'day_pairs': day_pairs, 'paths': paths, 'timepoints': list(self.timepoints),
             'meta': meta}
        with open(path, 'w') as f:
            json.dump(d, f, ensure_ascii=False)

    @staticmethod
    def from_json(index_path):
        import json
        delete_index = False
        if index_path.startswith('gs://'):
            index_path = wot.io.download_gs_url(index_path)
            delete_index = True

        with open(index_path, 'r') as r:
            tmap_info = json.load(r)
        if delete_index:
            os.remove(index_path)
        meta = pd.DataFrame(index=tmap_info['meta']['id'], data={'day': tmap_info['meta']['day']})
        timepoints = tmap_info['timepoints']
        paths = tmap_info['paths']
        day_pairs = tmap_info['day_pairs']
        tmaps = {}
        for i in range(len(day_pairs)):
            tmaps[tuple(day_pairs[i])] = paths[i]
        return TransportMapModel(tmaps=tmaps, meta=meta, timepoints=timepoints, day_pairs=day_pairs)

    @staticmethod
    def from_directory(tmap_out, with_covariates=False, cache=False):
        """
        Creates a wot.TransportMapModel from an output directory.

        Parameters
        ----------
        :param tmap_out:
        :param with_covariates:
        :return: TransportMapModel instance
        """
        tmap_dir, tmap_prefix = os.path.split(tmap_out)
        tmap_dir = tmap_dir or '.'
        tmap_prefix = tmap_prefix or "tmaps"
        tmaps = {}
        import re
        files = os.listdir(tmap_dir)
        day_regex = '([0-9]*\.?[0-9]+)'
        if with_covariates:
            pattern = re.compile(
                tmap_prefix + '_{}_{}_cv([a-zA-Z0-9]+)_cv([a-zA-Z0-9]+)[\.h5ad|\.loom]'.format(day_regex, day_regex))
        else:
            pattern = re.compile(tmap_prefix + '_{}_{}[\.h5ad|\.loom]'.format(day_regex, day_regex))
        for f in files:
            path = os.path.join(tmap_dir, f)
            if os.path.isfile(path):
                m = pattern.match(f)
                if m is not None:

                    try:
                        t1 = float(m.group(1))
                        t2 = float(m.group(2))
                        if with_covariates:
                            cv1 = m.group(3)
                            cv2 = m.group(4)
                            tmaps[(t1, t2, cv1, cv2)] = path
                        else:
                            tmaps[(t1, t2)] = path
                    except ValueError:
                        print('Unable to find day pair for ' + f)
                        pass

        if len(tmaps) is 0:
            raise ValueError('No transport maps found in ' + tmap_dir + ' with prefix ' + tmap_prefix)
        day_pairs = set()
        timepoints = set()
        tmap_keys = list(tmaps.keys())
        tmap_keys.sort(key=lambda x: x[0])
        meta = None
        for i in range(len(tmap_keys)):
            key = tmap_keys[i]
            t0 = key[0]
            t1 = key[1]
            day_pairs.add((t0, t1))
            timepoints.add(t0)
            timepoints.add(t1)
            if not with_covariates:
                path = tmaps[key]
                f = h5py.File(path, 'r')

                if path.endswith('.loom'):
                    rids = f['/row_attrs/id'][()].astype(str)
                    cids = f['/col_attrs/id'][()].astype(str) if i == len(tmap_keys) - 1 else None
                else:
                    rids = f['/obs']['index'][()].astype(str)
                    cids = f['/var']['index'][()].astype(str) if i == len(tmap_keys) - 1 else None
                rdf = pd.DataFrame(index=rids, data={'day': t0})
                cdf = pd.DataFrame(index=cids, data={'day': t1}) if cids is not None else None
                if meta is None:
                    meta = rdf
                else:
                    meta = pd.concat((meta, rdf), copy=False) if cdf is None else pd.concat((meta, rdf, cdf),
                                                                                            copy=False)

                f.close()
        timepoints = sorted(timepoints)
        return TransportMapModel(tmaps=tmaps, meta=meta, timepoints=timepoints, day_pairs=day_pairs, cache=cache)
