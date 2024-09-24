from netCDF4 import Dataset
import numpy as np
from global_land_mask import globe
import os
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

def seed_like_Luo(bathymetrypath, lon, lat, phylum, biomegridpath, initialmassdata, poctype):
    bmt = Dataset(bathymetrypath)
    biomes = np.load(biomegridpath)
    biome_name_map = ["HCSS", "LC", "HCPS", "COAST"]
    biome_N_count = [0, 0, 0, 0]
    for phi, lam in zip(lat, lon):
        i = np.where(bmt["lat"][:] == phi)[0][0]
        j = np.where(bmt["lon"][:] == lam)[0][0]
        biome_N_count[int(biomes[i, j])]+=1
    z = []
    m = []
    biome_seed = []
    for phi, lam in zip(lat, lon):
        i = np.where(bmt["lat"][:] == phi)[0][0]
        j = np.where(bmt["lon"][:] == lam)[0][0]
        depth = np.abs(bmt["topo"][i, j])
        b = int(biomes[i, j])
        biome_seed.append(b)
        biome = biome_name_map[b]
        m.append(initialmassdata[phylum.lower()][poctype][biome] / biome_N_count[b])
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
    return np.asarray(m), np.asarray(z), np.asarray(biome_seed, dtype = np.int32)

class RectangleSeed:

    def __init__(self, skip = 1):
        self.lat = np.arange(-90, 90, 1)
        self.lon = np.arange(-180, 180, 1)
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
        
        # if not os.path.isfile("lons{}.txt".format(self.s)):
        #     np.savetxt("lons{}.txt".format(self.s), lons)
        #     np.savetxt("lats{}.txt".format(self.s), lats)
        # lats, lons = np.meshgrid(lats, lons)
        return lons, lats
    @staticmethod
    def mask_grid(lons, lats, bathymetrypath, temperaturepath, biomegridpath):
        depth = Dataset(bathymetrypath)["topo"][:]
        tmp = Dataset(temperaturepath)["thetao"][0, :, :, :]
        biomegrid = np.load(biomegridpath)
        lon,lat = [], []
        lm = reader_global_landmask.Reader()
        # land = lm.get_variables("land_binary_mask", y =lats, x = lons)["land_binary_mask"]
        for i in range(len(lats)):
            for j in range(len(lons)):
                land = lm.get_variables("land_binary_mask", y =np.asarray([lats[i]]), x =np.asarray([lons[j]]))["land_binary_mask"]
                if not land:
                    if depth[i, j] > 0:
                        continue
                    if np.isnan(biomegrid[i, j]):
                        continue
                    # if tmp.mask[0, i, j]:
                    #     continue
                    # if np.any(tmp.data[:, i, j][np.invert(tmp.mask[:, i, j])] <= 0):
                    #     continue
                    lon.append(lons[j])
                    lat.append(lats[i])
            
        return lon, lat

class Seed:

    def __init__(self, seedtype, seeddata = None, ncfile = None, skip = 1, bathymetrypath = None, poctype = None,
                 phylum = None, lon = None, lat = None, z = None, outfile = None, biomegridpath = None, initialmassdata = None):
        
        if seedtype == "rectangle":
            if ncfile is None:
                raise NameError("Missing netCDF file for seeding.")
            
            o = RectangleSeed(skip)
            self.lon, self.lat = o.mask_grid(*o.grid(), bathymetrypath, ncfile, biomegridpath)
        elif seedtype == "line":
            if lon is None or lat is None:
                raise ValueError("Missing lon, lat values for seeding.")
            
            self.lon, self.lat = lon, lat
        
        if z is None:
            self.mass, self.z, self.biome = seed_like_Luo(bathymetrypath, self.lon, self.lat, phylum, biomegridpath, load_json_data(initialmassdata), poctype)
        else:
            self.biome = None
            self.mass = np.ones(len(self.lon))
            self.z = z
        
        
        if seeddata is None:
            # self.mass = np.ones(len(self.lon))
            self.r0 = np.ones(len(self.lon))
        # elif microbialdecaytype == "mass":
        #     self.mass = generate_mass_from_file(self.lon, self.lat, phylum, seeddata)
        #     self.r0 = None
        # else:
        #     mass, r = generate_mass_from_random_genus(self.lon, self.lat, phylum, seeddata, weightspath, biomegridpath)
        #     self.mass = mass
        #     self.r0 = r
        
        self.z = -np.abs(self.z)
        if outfile is not None:
            self.save_seed(outfile)
    
    def save_seed(self, outfile):
        data = {}
        data["lon"] = self.lon
        data["lat"] = self.lat
        data["mass"] = self.mass
        data["r0"] = self.r0
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















"""

def generate_mass_from_file(lon, lat, phylum:str, path:str): #Mass in units mg/m*2 Y
    data = load_data(path)
    keys = list(data.keys())
    found = False
    for i in range(len(keys)):
        if phylum.lower() == keys[i].lower():
            phylum = keys[i]
            found = True
            break
    
    if not found:
        raise NameError("Specified phylum is not in given data.")
        
    lons, lats, mass = data[phylum]["lons"], data[phylum]["lats"], data[phylum]["mass"]
    M = []

    for lam, phi in zip(lon, lat):
        i, j = np.argwhere((lons == lam) & (lats == phi)).ravel()
        M.append(mass[i, j])
    
    return np.asarray(M)

def get_radius(mass, genus): #Units of mass in mg
    
    def log_biom_eq(m, a, b):
        return (m / (10 ** a))**(1/b)
    
    def power_biom_eq(m, a, b):
        return (m / a) ** (1 / b)
    
    if genus == "Dolioletta":
        L = power_biom_eq(mass * 1000, 0.45, 2.3)
    elif genus == "Pryosoma":
        L = log_biom_eq(mass/1000/0.392, -1.8, 1.7)
    elif genus == "Brooksia":
        L = power_biom_eq(mass, 0.00011, 2)
    elif genus == "Cyclosalpa":
        L = power_biom_eq(mass, 0.005, 1.8)
    elif genus == "Iasis":
        L = power_biom_eq(mass, 0.0006, 2.3)
    elif genus == "Ihlea":
        L = power_biom_eq(mass, 0.0005, 2.4)
    elif genus == "Pegea":
        L = log_biom_eq(mass*1000, -0.11, 1.9)
    elif genus == "Salpa":
        L = log_biom_eq(mass * 1000, -0.2, 2.4)
    elif genus in ['Thalia', 'Thetys']:
        L = log_biom_eq(mass * 1000, 0.5, 1.2)
    elif genus in ['Craspedacusta', 'Clytia', 'Tiaropsis', 'Eutima', 'Catablema',
                   'Calycopsis', 'Halicreas', 'Octophialucium', 'Cosmetira',
                   'Aegina', 'Mitrocomella', 'Pegantha', 'Cunina',
                   'Cunina', 'Polyorchis', 'Calycopsis', 'Larsonia', 'Liriope',
                   'Eirene', 'Mitrocomella', 'Cyclocana', 'Laodicia', 'Eutonina',
                   'Melicertum', 'Cunina', 'Sarsia', 'Pandea', 'Bougainvillia', 'Sarsia',
                   'Dipleurosoma', 'Dichotomia', 'Aglaura', 'Amphinemum', 'Sarsia', 'Euphysa',
                   'Amphinema', 'Lizzia', 'Bougainvillea', 'Solmaris', 'Bougainvillea',
                   'Hybodcon', 'Rathkea', 'Solmaris', 'Nemopsis', 'Gossea', 'Eutima',
                   'Bougainvillea', 'Geryonia', 'Zygocanna', 'Aequorea', 'Tima', 'Solmissus',
                   'Ptychogena', 'Tima', 'Aequorea', 'Aequorea', 'Aequorea', 'Staurophora',
                   'Rhacostoma', 'Chiropsalmus', 'Eperetmus']:
        L = log_biom_eq(mass, -0.3, 0.06)
    elif genus in ['Linuche', 'Atolla', 'Atolla', 'Atorella', 'Paraphyllina']:
        L = power_biom_eq(mass / 0.1557, 0.33, 2.5)
    elif genus in ['Periphylla', 'Nausithoe', 'Tetraplatia']:
        L = (mass / 1000 / 0.196 / np.exp(-4.7)) ** (1 / 3.1)
    elif genus in ['Cotylorhiza', 'Rhopilema', 'Cassiopea', 'Stomolophus', 'Rhizostoma', 'Catostylus']:
        L = power_biom_eq(mass / 1000 / 0.0054, 0.08, 3.1) * 10
    elif genus == "Phyllorhiza":
        L = power_biom_eq(mass / 1000 / 0.336, 0.003, 2.6)
    elif genus in ['Aurelia', 'Phacellophora']:
        L = (mass / 0.03 - 1.4) / 0.5
    elif genus in ['Cyanea', 'Chrysaora', 'Cyanea', 'Drymonema']:
        L = power_biom_eq(mass / 1000 / 0.005, 0.19, 2.8) * 10
    elif genus in ['Pelagia', 'Chrysaora']:
        L = (mass * 1000  * np.exp(1.06)) ** (1 / 2.96)
    elif genus == "Beroe":
        L = power_biom_eq(mass / 0.0015, 1.77, 2.23)
    elif genus in ['Pleurobrachia', 'Callianira', 'Haeckelia', 'Mertensia', 'Lampea']:
        L = power_biom_eq(mass / 0.043 / 0.043, 0.7, 2.5)
    elif genus == "Bolinopsis":
        L = power_biom_eq(mass / 0.015, 0.13, 2.2)
    elif genus in ['Eurhamphaea', 'Deiopea', 'Leucothea']:
        L = power_biom_eq(mass / 0.009, 0.045, 2.3)
    elif genus == "Mnemiopsis":
        L = power_biom_eq(mass / 0.0015, 1.1, 2.8)
    elif genus == "Ocyropsis":
        L = power_biom_eq(mass / 0.023, 0.12, 2.2)
    else:
        raise NameError(genus + "DOES NOT HAVE A SPECIFIED BIOMETRIC EQUATION.")

    return L / 2 / 10 ** 3 #Return half of diameter in meters.

def get_CW_WW_conversion(phylum, genus):
    if phylum == "Chordata":
        return 0.58 / 100
    elif phylum == "Ctenophora":
        return 0.17 / 100
    else:
        if genus in ['Craspedacusta', 'Clytia', 'Tiaropsis', 'Eutima', 'Catablema',
                   'Calycopsis', 'Halicreas', 'Octophialucium', 'Cosmetira',
                   'Aegina', 'Mitrocomella', 'Pegantha', 'Cunina',
                   'Cunina', 'Polyorchis', 'Calycopsis', 'Larsonia', 'Liriope',
                   'Eirene', 'Mitrocomella', 'Cyclocana', 'Laodicia', 'Eutonina',
                   'Melicertum', 'Cunina', 'Sarsia', 'Pandea', 'Bougainvillia', 'Sarsia',
                   'Dipleurosoma', 'Dichotomia', 'Aglaura', 'Amphinemum', 'Sarsia', 'Euphysa',
                   'Amphinema', 'Lizzia', 'Bougainvillea', 'Solmaris', 'Bougainvillea',
                   'Hybodcon', 'Rathkea', 'Solmaris', 'Nemopsis', 'Gossea', 'Eutima',
                   'Bougainvillea', 'Geryonia', 'Zygocanna', 'Aequorea', 'Tima', 'Solmissus',
                   'Ptychogena', 'Tima', 'Aequorea', 'Aequorea', 'Aequorea', 'Staurophora',
                   'Rhacostoma', 'Chiropsalmus', 'Eperetmus']:
            return 0.51 / 100
        else:
            return 0.47 / 100

def generate_mass_from_random_genus(lon, lat, phylum:str, path:str, weights_path:str, biomegridpath):
    lons = np.arange(-180, 180, 1)
    lats = np.arange(-90, 90, 1)
    m1 = load_data(path)
    weights = load_data(weights_path)
    biomes = np.loadtxt(biomegridpath, dtype=np.int32)
    M = []
    r = []
    for lam, phi in zip(lon, lat):
        i = np.where(lons == lam)[0][0]
        j = np.where(lats == phi)[0][0]
        biome = biomes[j, i] #Find in which biome is (phi, lam).
        W = weights[phylum][str(biome)] #List of two lists, first element is genus name, second element is its probability.
        genus = np.random.choice(W[0], p = W[1])
        mass = m1[phylum][str(biome)] #Assign geometric mean mass of given biome.
        
        M.append(mass * 10 ** (-6) / get_CW_WW_conversion(phylum, genus)) # Factor to account for C -> WW conversion and mg -> kg.
        r.append(get_radius(mass, genus))
    
    return np.asarray(M), np.asarray(r)

"""