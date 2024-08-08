from netCDF4 import Dataset
import numpy as np
from global_land_mask import globe
import os

class RectangleSeed:

    def __init__(self, ncfile, skip = 1):
        X = Dataset(ncfile)
        self.lat = X["lat"][:]
        self.lon = X["lon"][:]
        self.s = skip
        
    
    def grid(self):
        if self.s > 1:
            latmin, latmax = np.min(self.lat), np.max(self.lat)
            lonmin, lonmax = np.min(self.lon), np.max(self.lon)

            lons = np.linspace(lonmin, lonmax, len(self.lon) // self.s)
            lats = np.linspace(latmin, latmax, len(self.lat) // self.s)
        elif self.s == 1:
            lons, lats = self.lon, self.lat
        else:
            raise ValueError("Paramater skip must be >=1.")
        
        # self.N = len(lons) * len(lats)
        
        if not os.path.isfile("lons{}.txt".format(self.s)):
            np.savetxt("lons{}.txt".format(self.s), lons)
            np.savetxt("lats{}.txt".format(self.s), lats)

        lons, lats = np.meshgrid(lons, lats)
        
        return lons, lats
    @staticmethod
    def mask_grid(lons, lats):
        lon,lat = [], []
        for i in range(lons.shape[0]):
            for j in range(lons.shape[1]):
                if globe.is_ocean(lon = lons[i, j], lat = lats[i, j]):
                    lon.append(lons[i, j])
                    lat.append(lats[i, j])
            
        return lon, lat

class LineSeed:

    def __init__(self, lon, lat):
        lonmin, lonmax, lonnum = [float(i) for i in lon.split("-")]
        latmin, latmax, latnum = [float(i) for i in lat.split("-")]
        lonnum = int(lonnum)
        latnum = int(latnum)
        if lonnum != latnum:
            num = min(lonnum, latnum)
            raise Warning("Seeding lonnum and latnum are different, taking the smaller value.")
        else:
            num = lonnum
        
        self.lons = np.linspace(lonmin, lonmax, num)
        self.lats = np.linspace(latmin, latmax, num)

class Seed:

    def __init__(self, seedtype, ncfile = None, skip = 1, lon = None, lat  = None):
        if seedtype == "rectangle":
            if ncfile is None:
                raise NameError("Missing netCDF file for seeding.")
            
            o = RectangleSeed(ncfile, skip)
            self.lon, self.lat = o.mask_grid(*o.grid())
        elif seedtype == "line":
            if lon is None or lat is None:
                raise ValueError("Missing lon, lat values for seeding.")
            
            o = LineSeed(lon, lat)
            self.lon, self.lat = o.lons, o.lats

