import numpy as np
import argparse

from model.logger import Logger
from model.plots_cleaned import Plot
from plotting.param_classifier import PlotParameters

def valid_figsize(s):
    if ":" not in s:
        raise argparse.ArgumentTypeError(f"not a valid figsize format: {s!r}")
    return tuple([float(i) for i in s.split(":")])

def main():
    p = argparse.ArgumentParser()
    p.add_argument("method", type = str)
    p.add_argument("-skip", default = 1, type = int)
    p.add_argument("-f1", "--file1", type = str, required=True)
    p.add_argument("-f2", "--file2", type = str)
    p.add_argument("-fs", "--fig_size", default = (20, 20), type=valid_figsize, help="Figsize - format H:W")
    p.add_argument("-fts", "--fontsize", type = int, default = 17)
    args = p.parse_args()
    plot(**vars(args))



def plot(**kwargs):
    plotmethod = kwargs.pop("method")
    skip = kwargs.pop("skip")
    params = PlotParameters(kwargs)
    p = Plot(**params.object_init)
    lons = np.loadtxt(f"lons" + str(skip) + ".txt")
    lats = np.loadtxt(f"lats" + str(skip) + ".txt")
    getattr(p, plotmethod)(lons = lons, lats = lats, **params.method)


if __name__ == '__main__':
    main()