import numpy as np
import argparse

from carbondrift.models.logger import Logger
from carbondrift.models.plots import Plot

def valid_figsize(s:str):
    if ":" not in s:
        raise argparse.ArgumentTypeError(f"not a valid figsize format: {s!r}")
    return tuple([float(i) for i in s.split(":")])

def configure_title(s:str):
    return s.replace("_", " ")

def configure_suptitles(s:str):
    return s.split(",")

def configure_locations(s:str):
    split = s.split(",")
    proper_format_locs = []
    for loc in split:
        i, j = loc.split(":")
        if "m" in i:
            lon = -float(i.strip("m"))
        else:
            lon = float(i)
        if "m" in j:
            lat = -float(j.strip("m"))
        else:
            lat = float(j)
        proper_format_locs.append((lon, lat))
    return proper_format_locs

def configure_clipping(s:str):
    if ":" not in s:
        raise argparse.ArgumentTypeError(f"not a valid clip format: {s!r}")
    clipmin, clipmax = s.split(":")
    if "m" in clipmin:
        clipmin = -float(clipmin.strip("m"))
    else:
        clipmin = float(clipmin)
    if "m" in clipmax:
        clipmax = -float(clipmax.strip("m"))
    else:
        clipmax = float(clipmax)
    return [clipmin,clipmax]

def configure_linestyles(s:str):
    return s.split(",")

def configure_ax_limits(s:str):
    if ":" not in s:
        raise argparse.ArgumentTypeError(f"not a valid xlim/ylim format: {s!r}")
    x1, x2 = s.split(":")
    if "m" in x1:
        x1 = -float(x1.strip("m"))
    else:
        x1 = float(x1)
    if "m" in x2:
        x2 = -float(x2.strip("m"))
    else:
        x2 = float(x2)
    return [x1, x2]

def configure_color(s:str):
    return s.split(",")

def configure_legend(s:str):
    return s.split(",")

def main():
    p = argparse.ArgumentParser(fromfile_prefix_chars='@')
    p.add_argument("method", type = str, help = "Ploting method to call.")
    p.add_argument("files", nargs = "+", type = str, help = "NetCDF filepaths to simulations.")
    p.add_argument("-lons", type = str, help=".txt file of longitude seed. If None -180:180:1")
    p.add_argument("-lats", type = str, help=".txt file of latitude seed. If None -90:90:1")
    p.add_argument("-cmap", type = str, required=False)
    p.add_argument("-t", "--title", type = str, required=False,
                   help = "Figure title - format aaa_bbb_ccc")
    p.add_argument("-abs", "--absolute", action="store_true",
                   help = "Plot absolute difference.")
    p.add_argument("-diff", action="store_true",
                   help="Plot difference of first two files.")
    p.add_argument("-add", action="store_true",
                   help="Plot sum of all files.")
    p.add_argument("-o", "--outfile", type = str)
    p.add_argument("-fs", "--figsize", default = (20, 20), type=valid_figsize,
                   help="Figsize - format W:H")
    p.add_argument("-fts", "--fontsize", type = int, default = 17)
    p.add_argument("-fwgt", "--fontweight", type = str, default="normal")
    p.add_argument("-d", "--depth", type = float, default=-200,
                   help = "Depth at which calculation takes place. If >= 5000 (<=-5000), it is set to sea floor.")
    p.add_argument("-clip", type = configure_clipping, help = "Clip values - format Vmin:Vmax. Use m for minus sign.")
    p.add_argument("-sh", "--shrink", type = float, default=1, help = "Shrink colorbar value.")
    p.add_argument("-loc", "--locations", type = str,
                   help = "Specific plotting locations - format lon1:lat1,lon2:lat2...")
    p.add_argument("-p1", "--prop1", type = str, help="Property to plot on x axis.")
    p.add_argument("-p2", "--prop2", type = str, help="Property to plot on y axis.")
    p.add_argument("-cblabel", "--colorbarlabel", type = str)
    p.add_argument("-xlabel", type = str)
    p.add_argument("-ylabel", type = str)
    p.add_argument("-xlim", type = configure_ax_limits, help = "Xlim - format xmin:xmax")
    p.add_argument("-ylim", type = configure_ax_limits, help = "Ylim - format ymin:ymax")
    p.add_argument("-lw", "--linewidth", type = float)
    p.add_argument("-legend", type = configure_legend, help="Legend - format label1,label2")
    p.add_argument("-color", type = configure_color, help = "Plot color - format color1,color2,color3...")
    p.add_argument("-ls", "--linestyle", type = configure_linestyles, help = "Plot linestlyes - format ls1,ls2,ls3...")
    p.add_argument("-bins", type = int, help = "Histogram bins.")
    p.add_argument("-group", type = str, default="none", 
                   choices=["none", "biome", "lonmean", "latmean"], help = "Group cells by.")
    p.add_argument("-diffidx", type = int)
    p.add_argument("-mc", "--martincurve", type = float,
                   help="Martin curve coefficient. If None, it will not be plotted.")
    p.add_argument("-supt", "--suptitle", type = str,
                   help = "Titles for specific axes. Format ax1title,ax2title...")
    p.add_argument("-mhwpath", "--mhwintesitypath", type = str)
    p.add_argument("-dpi", default=300, type=int)

    args = p.parse_args()
    return plot(**vars(args))

def plot(**kwargs):
    plotmethod = kwargs.pop("method")
    if not hasattr(Plot, plotmethod):
        methods_list = [method for method in dir(Plot) if callable(getattr(Plot, method)) and not method.startswith("__")]
        closest = methods_list[np.argmin([len(i.strip(plotmethod)) for i in methods_list])]
        raise ValueError(f"Module Plot doesn't have specified method {plotmethod}. Did you mean {closest}?")
    
    files = kwargs.pop("files", [])
    if len(files) == 0:
        raise ValueError(f"No simulation files have been provided.")

    lons = kwargs.pop("lons")
    lats = kwargs.pop("lats")

    kwargs["lons"] = np.loadtxt(lons) if lons is not None else None
    kwargs["lats"] = np.loadtxt(lats) if lats is not None else None

    if kwargs["title"] is not None:
        kwargs["title"] = configure_title(kwargs["title"])
    
    if kwargs["suptitle"] is not None:
        kwargs["suptitle"] = configure_suptitles(kwargs["suptitle"])
    
    if kwargs["locations"] is not None:
        kwargs["locations"] = configure_locations(kwargs["locations"])

    p = Plot(*files, **kwargs)
    return getattr(p, plotmethod)()


if __name__ == '__main__':
    print(main())