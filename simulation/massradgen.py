import numpy as np
import numpy.random as rnd
import pickle

def load_data(path):
    with open(path, "rb") as f:
        data = pickle.load(f)
    return data

def generate_masses(lon:np.array, species:str, path:str):
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
        
    mass = data[species]["mass"]
    mass[np.isnan(lon)] = np.nan
    return mass

def generate_radii(lon, lat, mass, species:str):
    new_mass = []
    new_lons = []
    new_lats =  []
    radius = []
    sample_function = None # TODO Create random function for sampling for different species.
    for i, M in enumerate(mass):
        m = 0
        while m < M:
            r = sample_function()
            mi = k * r
            m += mi
            new_lats.append(lat[i])
            new_lons.append(lon[i])
            new_mass.append(mi)
            radius.append(r)
        if m > M:
            m -= mi
            mi = M - m
            new_mass[-1] = mi
            radius[-1] = mi / k
    return new_mass, new_lons, new_lats, radius



#TODO Generate masses and radii from a distribution of radii.