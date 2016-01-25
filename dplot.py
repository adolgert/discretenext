"""dplot

Plots discrete time data.

Usage:
  dplot.py [-v] [-q] evince <aglob> <bglob>
  dplot.py [-v] [-q] interrupt <aglob> <bglob>

Options:
  -h, --help  Show this screen.
  -v          Verbose logging.
  -q          Quiet logging.
"""
import logging
import math
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


def interrupt(file_globs):
    interrupt_n(file_globs, 4)
    interrupt_n(file_globs, 5)

def interrupt_n(file_globs, n):
    summary=list()
    for fg in file_globs:
        files=glob.glob(fg)
        fstream=list()
        for hf in files:
            try:
                fstream.append(h5py.File(hf, "r"))
            except:
                logger.error("Could not open {0}".format(hf))
        summary.append(datafiles.accumulate_array(fstream, "{0}fire".format(n)))
        fsum=np.sum(summary[-1])
        dsum=np.sum(datafiles.accumulate_array(fstream, "{0}disable".format(n)))
        total=fsum+dsum
        p=fsum/total
        variance=total*p*(1-p)
        logger.info("fire ratio {0} N {1} sqrt(var)/total {2} name {3}".format(
                p, total, math.sqrt(variance)/total, fg))

    dt=0.01
    fig=plt.figure(1, figsize=(8, 5))
    ax=fig.add_subplot(111)
    ax.set_title("Holding Time")
    end=0
    for arr in summary:
        cutoff=np.max(arr)/100
        end=max(end, dt*np.where(arr>cutoff)[0][-1])
    x=dt*np.array(range(summary[0].shape[0]), np.int)
    ax.set_xlim(0, 1.1*end)
    for sums in summary:
        plt.plot(x, sums/np.sum(sums), ".", ms=1)
    plt.tight_layout()
    fname="interrupt{0}.png".format(n)
    logger.info("writing "+fname)
    plt.savefig(fname, format="png")
    plt.clf()


def evince(file_globs):
    """
    Plot h_{ij} for given distributions.
    """

    ifire=list()
    rfire=list()
    summary=list()
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
        summary.append(datafiles.accumulate_array(fstream, "summary"))

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
        plt.tight_layout()
        plt.savefig("Firing{0}.pdf".format(name), format="pdf")
        plt.clf()

    fig=plt.figure(1, figsize=(3, 2))
    ax=fig.add_subplot(111)
    ax.set_title("Occupancy Fraction")
    x=np.array(range(summary[0].shape[0]), np.int)
    for sums in summary:
        plt.plot(x, sums/np.sum(sums))
    plt.tight_layout()
    plt.show()
    plt.savefig("occupancy.pdf", format="pdf")
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
    if arguments["interrupt"]:
        interrupt([arguments["<aglob>"], arguments["<bglob>"]])