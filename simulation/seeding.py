from netCDF4 import Dataset
import numpy as np
import pickle
import pandas as pd
import json
from opendrift.readers import reader_global_landmask


def load_data(path):
    with open(path, "rb") as f:
        data = pickle.load(f)
    return data

def load_json_data(file):
    with open(file, "r") as f:
        data = json.load(f)
    return data

def luo_seed(bathymetrypath, areapath, lon, lat, phylum, biomegridpath, initialmassdata, poctype):
    bmt = Dataset(bathymetrypath)
    area = np.load(areapath)
    biomes = np.load(biomegridpath)
    biome_name_map = ["HCSS", "LC", "HCPS", "COAST"]
    biome_grid_cell_area = [[], [], [], []]
    
    for phi in np.arange(-90, 90, 1):
        for lam in np.arange(-180, 180, 1):
            i = np.where(bmt["lat"][:] == phi)[0][0]
            j = np.where(bmt["lon"][:] == lam)[0][0]
            if np.isnan(biomes[i, j]):
                continue
            biome_grid_cell_area[int(biomes[i, j])].append(area[i, j])
    biome_areas = [sum(i) for i in biome_grid_cell_area]
    # print(np.asarray(biome_areas) / 10 ** 14)

    z = []
    m = []
    biome_seed = []
    biome_counter = [0, 0, 0, 0]
    for phi, lam in zip(lat, lon):
        i = np.where(bmt["lat"][:] == phi)[0][0]
        j = np.where(bmt["lon"][:] == lam)[0][0]
        depth = np.abs(bmt["topo"][i, j])
        b = int(biomes[i, j])
        biome_seed.append(b)
        biome = biome_name_map[b]
        if type(initialmassdata) == float:
            m.append(1)
        else:
            m.append(initialmassdata[phylum.lower()][poctype][biome] * biome_grid_cell_area[b][biome_counter[b]] / biome_areas[b])
        if depth < 30:
            z.append(0)
        elif depth < 60:
            if phylum.lower() == "chordata":
                z.append(-25)
            else:
                z.append(-10)
        else:
            if phylum.lower() == "chordata":
                z.append(-50)
            else:
                z.append(-20)
        biome_counter[b] += 1
    
    return np.asarray(m), np.asarray(z), np.asarray(biome_seed, dtype = np.int32)

class RectangleSeed:

    def __init__(self, latmnin = -90, latmax = 90, lonmin = -180, lonmax = 180, dx = 1, dy = 1):
        self.lat = np.arange(latmnin, latmax, dx)
        self.lon = np.arange(lonmin, lonmax, dy)
    
    @staticmethod
    def mask_grid(lons, lats, bathymetrypath, biomegridpath):
        depth = Dataset(bathymetrypath)["topo"][:]
        # tmp = Dataset(temperaturepath)["thetao"][0, :, :, :]
        biomegrid = np.load(biomegridpath)
        lon,lat = [], []
        lm = reader_global_landmask.Reader()
        # land = lm.get_variables("land_binary_mask", y =lats, x = lons)["land_binary_mask"]
        for i in range(len(lats)):
            for j in range(len(lons)):
                land = lm.get_variables("land_binary_mask", y =np.asarray([lats[i]]), x =np.asarray([lons[j]]))["land_binary_mask"]
                k = np.where(np.arange(-90, 90, 1) == lats[i])[0][0]
                l = np.where(np.arange(-180, 180, 1) == lons[j])[0][0]
                if not land:
                    if depth[k, l] > 0:
                        continue
                    if np.isnan(biomegrid[k, l]):
                        continue
                    # if tmp.mask[0, i, j]:
                    #     continue
                    # if np.any(tmp.data[:, i, j][np.invert(tmp.mask[:, i, j])] <= 0):
                    #     continue
                    lon.append(lons[j])
                    lat.append(lats[i])
            
        return lon, lat

class Seed:

    def __init__(self, latmin = -90, latmax = 90, lonmin = -180, lonmax = 180, dx = 1, dy = 1,
                 bathymetrypath = None, poctype = None, phylum = None, areapath = None,
                 outfile = None, biomegridpath = None, initialmassdata = None):

        o = RectangleSeed(latmin, latmax, lonmin, lonmax, dx, dy)
        self.lon, self.lat = o.mask_grid(o.lon, o.lat, bathymetrypath, biomegridpath)
        
        if initialmassdata is None:
            mass = 1.
        else:
            mass = load_json_data(initialmassdata)
        self.mass, self.z, self.biome = luo_seed(bathymetrypath, areapath, self.lon, self.lat, phylum, biomegridpath, mass, poctype)
        
        self.z = -np.abs(self.z)
        if outfile is not None:
            self.save_seed(outfile)
    
    def save_seed(self, outfile):
        data = {}
        data["lon"] = self.lon
        data["lat"] = self.lat
        data["mass"] = self.mass
        data["z"] = self.z
        if self.biome is not None:
            data["origin_marker"] = self.biome
        df = pd.DataFrame(data)
        df.to_pickle(outfile)

class SeedFromFile:
    def __init__(self, file):
        data = load_data(file)
        for keys, values in data.items():
            setattr(self, keys, values.values)