from netCDF4 import Dataset
import numpy as np
from global_land_mask import globe


def rectangle_seed(pathtofile, save = False, cut = 1):
    X = Dataset(pathtofile)
    lat = X["lat"][:]
    lon = X["lon"][:]

    number = len(lat) * len(lon)
    latmin, latmax = np.min(lat), np.max(lat)
    lonmin, lonmax = np.min(lon), np.max(lon)

    lons = np.linspace(lonmin, lonmax, len(lon) // cut)
    lats = np.linspace(latmin, latmax, len(lat) // cut)

    number = len(lons) * len(lats)
    
    if save:
        import os
        curdir = os.getcwd()
        os.chdir("/home/perharic/Documents")
        np.savetxt("lons{}.txt".format(cut), lons)
        np.savetxt("lats{}.txt".format(cut), lats)
        os.chdir(curdir)

    lons, lats = np.meshgrid(lons, lats)

    return lons, lats, number

def mask_seed(lons, lats):
    for i in range(lons.shape[0]):
        for j in range(lons.shape[1]):
            if globe.is_land(lon = lons[i, j], lat = lats[i, j]):
                lons[i, j] = np.nan
                lats[i, j] = np.nan
    
    return lons, lats

