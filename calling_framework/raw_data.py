"""module provides functions to get raw data from nanopore sequencing.
Usually, data is provided in form of a .fast5-file; however sometimes
one might encounter Robjects or other data-types from foreign sources.
"""

import numpy as np
import h5py
import rpy2.robjects as robjects

def raw_from_fast5(path, channel, start, end):
    """ read the raw data from a fast5 file.

    Args:
       path: path to file
       channels: list with channel-ids to extracts

    Returns:
        dict[channel_id] -> np.array(raw_values)
    """

    f = h5py.File(path, r)
    raw_data = dict()
    for i in channels:
        raw_data[i] = np.array(f["/Raw/Channel_{0}/Signal".format(i)])

    return raw_data


