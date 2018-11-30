import anndata

""" Reads in a gct file .

The main method is parse. parse_into_3 creates the row
metadata, column metadata, and data dataframes, while the
assemble_multi_index_df method in GCToo.py assembles them.

1) Example GCT v1.3:
        ----- start of file ------
        #1.3
        96 36 9 15
        ---------------------------------------------------
        |id|        rhd          |          cid           |
        ---------------------------------------------------
        |  |                     |                        |
        |c |                     |                        |
        |h |      (blank)        |      col_metadata      |
        |d |                     |                        |
        |  |                     |                        |
        ---------------------------------------------------
        |  |                     |                        |
        |r |                     |                        |
        |i |    row_metadata     |          data          |
        |d |                     |                        |
        |  |                     |                        |
        ---------------------------------------------------
        ----- end of file ------

        Notes:
        - line 1 of file ("#1.3") refers to the version number
        - line 2 of file ("96 36 9 15") refers to the following:
                -96 = number of data rows
                -36 = number of data columns
                -9 = number of row metadata fields (+1 for the 'id' column -- first column)
                -15 = number of col metadata fields (+1 for the 'id' row -- first row)
        - Once read into a DataFrame, col_metadata_df is stored as the transpose of how it looks in the gct file.
                That is, col_metadata_df.shape = (num_cid, num_chd).

2) Example GCT v1.2

        ----- start of file ------
        #1.2
        96 36
        -----------------------------------------------
        |"NAME" |"Description"|          cid           |
        -----------------------------------------------
        |   r   |             |                        |
        |   i   |             |                        |
        |   d   |row_metadata |         data           |
        |       |             |                        |
        |       |             |                        |
        -----------------------------------------------

        ----- end of file ------
        Notes:
        - line 1 of file ("#1.3") refers to the version number
        - line 2 of file ("96 36 9 15") refers to the following:
                -96 = number of data rows
                -36 = number of data columns

"""

import numpy as np
import pandas as pd

__author__ = "Lev Litichevskiy, Oana Enache"
__email__ = "lev@broadinstitute.org"


def read_gct(file_path):
    """ The main method.

    Args:
        - file_path (string): full path to gct file you want to parse

    Returns:
        anndata.AnnData

    """

    # Read version and dimensions
    (version, num_data_rows, num_data_cols,
     num_row_metadata, num_col_metadata) = read_version_and_dims(file_path)

    # Read in metadata and data
    (row_metadata, col_metadata, data) = parse_into_3(
        file_path, num_data_rows, num_data_cols,
        num_row_metadata, num_col_metadata)

    row_metadata.index.name = None
    col_metadata.index.name = None
    row_metadata.columns.name = None
    col_metadata.columns.name = None
    return anndata.AnnData(data, row_metadata, col_metadata)


def read_version_and_dims(file_path):
    # Open file
    f = open(file_path, "r")

    # Get version from the first line
    version = f.readline().strip().lstrip("#")

    if version not in ["1.3", "1.2"]:
        err_msg = ("Only GCT1.2 and 1.3 are supported. The first row of the GCT " +
                   "file must simply be (without quotes) '#1.3' or '#1.2'")

        raise Exception(err_msg.format(version))

    # Convert version to a string
    version_as_string = "GCT" + str(version)

    # Read dimensions from the second line
    dims = f.readline().strip().split("\t")

    # Close file
    f.close()

    # Check that the second row is what we expect
    if version == "1.2" and len(dims) != 2:
        error_msg = "GCT1.2 should have 2 dimension-related entries in row 2. dims: {}"

        raise Exception(error_msg.format(dims))
    elif version == "1.3" and len(dims) != 4:
        error_msg = "GCT1.3 should have 4 dimension-related entries in row 2. dims: {}"

        raise Exception(error_msg.format(dims))

    # Explicitly define each dimension
    num_data_rows = int(dims[0])
    num_data_cols = int(dims[1])
    if len(dims) == 4:
        num_row_metadata = int(dims[2])
        num_col_metadata = int(dims[3])
    else:
        num_row_metadata = 1
        num_col_metadata = 0

        # Return version and dimensions
    return version_as_string, num_data_rows, num_data_cols, num_row_metadata, num_col_metadata


def parse_into_3(file_path, num_data_rows, num_data_cols, num_row_metadata, num_col_metadata):
    # Read the gct file beginning with line 3
    full_df = pd.read_csv(file_path, sep="\t", header=None, skiprows=2, dtype=str)

    # Check that full_df is the size we expect
    assert full_df.shape == (num_col_metadata + num_data_rows + 1,
                             num_row_metadata + num_data_cols + 1), (
        ("The shape of full_df is not as expected: data is {} x {} " +
         "but there are {} row meta fields and {} col fields").format(
            num_data_rows, num_data_cols, num_row_metadata, num_col_metadata))

    # Assemble metadata dataframes
    row_metadata = assemble_row_metadata(full_df, num_col_metadata, num_data_rows, num_row_metadata)
    col_metadata = assemble_col_metadata(full_df, num_col_metadata, num_row_metadata, num_data_cols)

    # Assemble data dataframe
    data = assemble_data(full_df, num_col_metadata, num_data_rows, num_row_metadata, num_data_cols)

    # Return 3 dataframes
    return row_metadata, col_metadata, data


def assemble_row_metadata(full_df, num_col_metadata, num_data_rows, num_row_metadata):
    # Extract values
    row_metadata_row_inds = range(num_col_metadata + 1, num_col_metadata + num_data_rows + 1)
    row_metadata_col_inds = range(1, num_row_metadata + 1)
    row_metadata = full_df.iloc[row_metadata_row_inds, row_metadata_col_inds]

    # Create index from the first column of full_df (after the filler block)
    row_metadata.set_index(full_df.iloc[row_metadata_row_inds, 0], inplace=True)
    # Create columns from the top row of full_df (before cids start)
    row_metadata.columns = full_df.iloc[0, row_metadata_col_inds].values

    # Rename the index name and columns name
    # row_metadata.index.name = row_index_name
    # row_metadata.columns.name = row_header_name

    # Convert metadata to numeric if possible
    row_metadata = row_metadata.apply(lambda x: pd.to_numeric(x, errors="ignore"))

    return row_metadata


def assemble_col_metadata(full_df, num_col_metadata, num_row_metadata, num_data_cols):
    # Extract values
    col_metadata_row_inds = range(1, num_col_metadata + 1)
    col_metadata_col_inds = range(num_row_metadata + 1, num_row_metadata + num_data_cols + 1)
    col_metadata = full_df.iloc[col_metadata_row_inds, col_metadata_col_inds]

    # Transpose so that samples are the rows and headers are the columns
    col_metadata = col_metadata.T

    # Create index from the top row of full_df (after the filler block)
    col_metadata.set_index(full_df.iloc[0, col_metadata_col_inds], inplace=True)

    # Create columns from the first column of full_df (before rids start)
    col_metadata.columns = full_df.iloc[col_metadata_row_inds, 0].values

    # # Rename the index name and columns name
    # col_metadata.index.name = 'id'
    # col_metadata.columns.name = 'description'

    # Convert metadata to numeric if possible
    col_metadata = col_metadata.apply(lambda x: pd.to_numeric(x, errors="ignore"))

    return col_metadata


def assemble_data(full_df, num_col_metadata, num_data_rows, num_row_metadata, num_data_cols):
    # Extract values
    data_row_inds = range(num_col_metadata + 1, num_col_metadata + num_data_rows + 1)
    data_col_inds = range(num_row_metadata + 1, num_row_metadata + num_data_cols + 1)
    data = full_df.iloc[data_row_inds, data_col_inds].values
    return data.astype(np.float64)
