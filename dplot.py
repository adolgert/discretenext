"""dplot

Plots discrete time data.

Usage:
  dplot.py [-v] [-q] evince <aglob> <bglob>

Options:
  -h, --help  Show this screen.
  -v          Verbose logging.
  -q          Quiet logging.
"""
import logging
import os
import glob
import numpy as np
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
import docopt
import h5py

import datafiles
import pdfedit

logger=logging.getLogger("dplot")


def evince(file_globs):
    """
    Plot h_{ij} for given distributions.
    """

    ifire=list()
    rfire=list()
    for fg in file_globs:
        files=glob.glob(fg)
        fstream=list()
        for hf in files:
            try:
                fstream.append(h5py.File(hf, "r"))
            except:
                logger.error("Could not open {0}".format(hf))
        ifire.append(datafiles.accumulate_array(fstream, "Ifire"))
        rfire.append(datafiles.accumulate_array(fstream, "Rfire"))

    for arrays, name in zip([ifire, rfire], ["Infection", "Recovery"]):
        fig=plt.figure(1, figsize=(3, 2))
        ax=fig.add_subplot(111)
        ax.set_title("Firing {0}".format(name))
        dt=0.01
        end=0
        for arrlen in arrays:
            maxarr=np.max(arrlen)/100
            end=max(end, np.where(arrlen>maxarr)[0][-1])
        ax.set_xlim(0, (end+10)*dt)
        x=dt*np.array(range(arrays[0].shape[0]), dtype=np.double)
        for arr in arrays:
            plt.plot(x, arr/np.sum(arr))
        plt.show()
        plt.savefig("Firing{0}".format(name))
        plt.clf()


if __name__ == "__main__":
    arguments = docopt.docopt(__doc__, version="dplot 1.0")
    if arguments["-v"]:
        logging.basicConfig(level=logging.DEBUG)
    elif arguments["-q"]:
        logging.basicConfig(level=logging.ERROR)
    else:
        logging.basicConfig(level=logging.INFO)

    if arguments["evince"]:
        evince([arguments["<aglob>"], arguments["<bglob>"]])