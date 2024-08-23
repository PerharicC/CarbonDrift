from netCDF4 import Dataset
import numpy as np
from global_land_mask import globe
import os
import pickle


def load_data(path):
    with open(path, "rb") as f:
        data = pickle.load(f)
    return data

def generate_masses(lon, lat, species:str, path:str):
    data = load_data(path)
    keys = list(data.keys())
    found = False
    for i in range(len(keys)):
        if species.lower() == keys[i].lower():
            species = keys[i]
            found = True
            break
    
    if not found:
        raise NameError("Specified species is not in given data.")
        
    lons, lats, mass = data[species]["lons"], data[species]["lats"], data[species]["mass"]
    M = []

    for lam, phi in zip(lon, lat):
        i, j = np.argwhere((lons == lam) & (lats == phi)).ravel()
        M.append(mass[i, j])
    
    return M


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
        lons = lon.split(":")
        lats= lat.split(":")
        for i in range(3):
            if "m" in lons[i]:
                lons[i] = -float(lons[i].strip("m"))
            else:
                lons[i] = float(lons[i])
            if "m" in lats[i]:
                lats[i] = -float(lats[i].strip("m"))
            else:
                lats[i] = float(lats[i])
        lonmin, lonmax, lonnum = lons
        latmin, latmax, latnum = lats
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

    def __init__(self, seedtype, species, massdata, ncfile = None, skip = 1, lon = None, lat  = None, z = 0):
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
        self.mass = generate_masses(self.lon, self.lat, species, massdata)
        self.z = -np.abs(z)