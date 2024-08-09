import numpy as np
import argparse

from model.logger import Logger
from model.plots import Plot

def valid_figsize(s:str):
    if ":" not in s:
        raise argparse.ArgumentTypeError(f"not a valid figsize format: {s!r}")
    return tuple([float(i) for i in s.split(":")])

def configure_title(s:str):
    return s.replace("_", " ")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("method", type = str)
    p.add_argument("-f1", "--file1", type = str, required=True)
    p.add_argument("-f2", "--file2", type = str)
    p.add_argument("-lons", type = str, help="lons.txt file from simulation.")
    p.add_argument("-lats", type = str, help="lats.txt file from simulation.")
    p.add_argument("-cmap", type = str, required=False)
    p.add_argument("-t", "--title", type = str, required=False,
                   help = "Figure title - format aaa_bbb_ccc")
    p.add_argument("-abs", "--absolute", action="store_true",
                   help = "PLot absolute diffrence.")
    p.add_argument("-diff", action="store_false",
                   help="If set, don't plot difference of two files.")
    p.add_argument("-o", "--outfile", type = str)
    p.add_argument("-fs", "--figsize", default = (20, 20), type=valid_figsize,
                   help="Figsize - format H:W")
    p.add_argument("-fts", "--fontsize", type = int, default = 17)
    p.add_argument("-fwgt", "--fontweight", type = str, default="normal")
    p.add_argument("-d", "--depth", type = float, default=-200,
                   help = "Depth at which calculation takes placce. If >= 5000, it is set to sea floor.")
    p.add_argument("-clip", type = str, help = "Clip values - format Vmin:Vmax. Use m for minus sign.")
    p.add_argument("-sh", "--shrink", type = float, default=1)

    args = p.parse_args()
    plot(**vars(args))

def plot(**kwargs):
    plotmethod = kwargs.pop("method")

    lons = kwargs.pop("lons")
    lats = kwargs.pop("lats")

    kwargs["lons"] = np.loadtxt(lons) if lons is not None else None
    kwargs["lats"] = np.loadtxt(lats) if lats is not None else None

    if kwargs["title"] is not None:
        kwargs["title"] = configure_title(kwargs["title"])

    p = Plot(**kwargs)
    getattr(p, plotmethod)()


if __name__ == '__main__':
    main()