## Version 1.0.8.post2, August 24, 2020
- Use anndata for reading and writing datasets

## Version 1.0.8.post1, April 21, 2020
- Support anndata 0.7

## Version 1.0.8, October 16, 2019
- Compute t-test for diff_exp. Store results in long format.

## Version 1.0.7, October 7, 2019
- Fixed invalid import statement

## Version 1.0.6, October 7, 2019
- Use [pegasus](https://pegasus.readthedocs.io/) for reading and writing datasets.
- Removed neighborhood_graph command (recommended to use pegasus instead)

## Version 1.0.5, June 27, 2019
- Fixed trajectory error when starting from a middle (not start or end) time point
- Fixed bug when reading transport maps with only two time points

## Version 1.0.4, June 11, 2019
- Ensure distance matrix is C-contiguous for earth mover's distance

## Version 1.0.3, June 4, 2019
- Output gene set scores in one file in gene_set_scores command
- Compute distance matrix in parallel
- Updated arguments in convert_matrix command
- Fixed writing gct or complete text file with empty columns

## Version 1.0.2, May 16, 2019
- Compute fraction expressed ratio in diff_exp tool and updated tool arguments

## Version 1.0.1, May 13, 2019 
- Fixed argument names in gene_set_scores and compute_all_transport_maps functions

## Version 1.0.0, May 13, 2019 
- First official wot release


