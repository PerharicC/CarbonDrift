from netCDF4 import Dataset
import numpy as np
from global_land_mask import globe
import os
import pickle


def load_data(path):
    with open(path, "rb") as f:
        data = pickle.load(f)
    return data

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
    """Biometric equations compiled for different genuses from Lucas et al 2011."""
    
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
        
        M.append(mass * 8) # Factor 8 to account for C -> DW conversion
        r.append(get_radius(mass, genus))
    
    return np.asarray(M), np.asarray(r)

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

    def __init__(self, seedtype, phylum, massdata, ncfile = None, skip = 1,
                 lon = None, lat = None, z = 0, massgen_type = None, microbialdecaytype = "mass",
                 weightspath = None, biomegridpath = None):
        
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
        
        if massgen_type is None:
            if microbialdecaytype == "mass":
                massgen_type = "from_file"
            elif microbialdecaytype =="area":
                massgen_type = "random_genus"
        
        if massgen_type == "from_file":
            self.mass = generate_mass_from_file(self.lon, self.lat, phylum, massdata)
            self.r0 = None
        elif massgen_type == "random_genus":
            mass, r = generate_mass_from_random_genus(self.lon, self.lat, phylum, massdata, weightspath, biomegridpath)
            self.mass = mass
            self.r0 = r
        else:
            self.mass = np.ones(len(self.lon))
            self.r0 = None
        
        self.z = -np.abs(z)