"""datafiles

This is a testable script to find datasets in files.

Usage:
  datafiles.py

Options:
  -h --help   Print this screen.
"""

import os.path
import logging
import glob
import numpy as np
import h5py

logger=logging.getLogger("datafiles")


def files_matching(condition, directory):
    """
    Given a directory and a conditional on the HDF5 group
    for run0000 in that each file, return a list of open HDF5
    file handles.
    """
    chosen_files=list()
    hdf_files=glob.glob(os.path.join(directory, "*.hdf5"))
    for f in sorted(hdf_files):
        try:
            infile=h5py.File(f, "r")
            grp=infile["run0000"]
            if condition(grp):
                logger.debug("files_matching appending {0}".format(f))
                chosen_files.append(infile)
            else:
                infile.close()
        except OSError as e:
            logger.debug("Could not open {0}".format(f))
    return chosen_files


def accumulate_array(h5stream_iter, array_name):
    array=None
    for h5f in h5stream_iter:
        for run in [x for x in h5f.keys() if x.startswith("run")]:
            grp=h5f[run]
            if array_name in grp.keys():
                if array is None:
                    shape=grp[array_name].shape
                    array=np.zeros(shape, dtype=grp[array_name].dtype)
                if array.shape!=grp[array_name].shape:
                    logger.error("File {0} has shape {1}".format(
                            h5f.filename, grp[array_name].shape))
                array+=grp[array_name]
            else:
                pass # skip missing arrays
    if array is None:
        logger.error("Couldn't find {0} in any file.".format(array_name))
    return array




if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    naadsmish=lambda x: (x.attrs["algorithm"]=="naadsmish" and
            "config" in x.attrs and
            x.attrs["config"]=="rising")
    for f in files_matching(naadsmish, "."):
        print(f.filename)
        print(f["run0000"].attrs["runtime"])
    print(accumulate_array(files_matching(naadsmish, "."), "duration"))