# -*- coding: utf-8 -*-


def split_anndata(dataset, metadata):
    """
    Split AnnData into sub-datasets according to a obs metadata

    Parameters
    ----------
    metadata : str
        The metadata to use for the split

    Returns
    -------
    splits : dict of t: anndata.AnnData
        Dictionary of datasets. t is the type of the 'metadata' column.
        Each cell in splits[k] has its 'metadata' column constant to k.

    Raises
    ------
    ValueError
        If the metadata is not present
    """

    if metadata not in dataset.obs.columns:
        raise ValueError("Cannot split on '{}' : column not present".format(metadata))

    def extract(group):
        return dataset[group.index]

    return {name: extract(group) for name, group in dataset.obs.groupby(metadata)}
