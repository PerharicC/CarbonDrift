import os

from carbondrift.models.areadecay.carbondrift import *
import numpy as np
from datetime import datetime, timedelta
from copy import copy
import xarray as xr
from netCDF4 import Dataset
import numpy.ma as ma

# import model.massdecay.carbondriftopen



class GridRun:
    """Run the CarbonDrift module in sections of initial lon and lat seperately and then merge the simulations together."""

    def __init__(self, **kwargs):

        """
        Parameters
        ----------
        starttime: datetime object
            Beginning of simulation (default datetime.now())
        reader: list of readers
            The readers to be used in the simulation (defualt None)
        lon: 1D/2D ndarray or float
            Initial longtitude of seeded particles (default 0)
        lat: 1D/2D ndarray or float
            Initial latitude of seeded particles (default 0)
        configure: dict
            Dictionary of all properties, which should be configured.
        deactivate_fragmentation: bool
            Boolean statement if fragmentation should be deactivated (default False)
        deactivate_horizontal_advection: bool
            Boolean statement if horizontal advection should be deactivated (default False)
        split_factor: int
            Number of splits (default 1)
        merge: bool
            Merge and delete all individual files (default True)
        + All key arguments for CarbonDrift.
        """

        self.time = kwargs.pop("starttime", datetime.now())
        self.reader = kwargs.pop("reader", None)
        self.lon = kwargs.pop("lon", 0)
        self.lat = kwargs.pop("lat", 0)
        self.config = kwargs.pop("configure", {})
        self.deactivate_frag = kwargs.pop("deactivate_fragmentation", False)
        self.deactivate_ha = kwargs.pop("deactivate_horizontal_advection", False)
        self.split_factor = kwargs.pop("split_factor", 1)
        self.merge = kwargs.pop("merge", True)
        self.o_kwargs = kwargs
    
    def run(self, **kwargs):

        if self.split_factor <= 1:
            o = CarbonDrift(**self.o_kwargs)
            if self.reader is not None:
                o.add_reader(self.reader)
            
            if self.deactivate_frag:
                o.deactivate_fragmentation()
            
            if self.deactivate_ha:
                o.deactivate_horizontal_advection()
            
            for key, value in self.config.items():
                o.set_config(key, value)
            try:
                mass = np.ones(self.lon.shape[0] * self.lat.shape[1])
            except AttributeError:
                mass = np.ones(len(self.lon))
            o.seed_elements(lon = self.lon, lat = self.lat, z = 0, time = self.time, mass = mass)

            o.run(**kwargs)
        else:    
            lons, lats = self.split()
            self.filename = copy(kwargs["outfile"]).split(".")[0]

            for i in range(len(lons)):
                o = CarbonDrift(**self.o_kwargs)

                if self.reader is not None:
                    o.add_reader(self.reader)
                
                if self.deactivate_frag:
                    o.deactivate_fragmentation()
                if self.deactivate_ha:
                    o.deactivate_horizontal_advection()
                
                for key, value in self.config.items():
                    o.set_config(key, value)

                o.seed_elements(lon = lons[i], lat = lats[i], z = 0, time = self.time, mass = np.ones(len(lons[i])))

                kwargs["outfile"] = self.filename + str(i) + ".nc"
                o.run(**kwargs)
                del o
            
            if self.merge:
                self.merge_files(len(lons))

    def split(self):
        num = self.lon.shape[0] * self.lon.shape[1]
        idx = num // self.split_factor

        lons = []
        lats = []
        lon = self.lon.flatten()
        lat = self.lat.flatten()
        for i in range(0, len(lon), idx):
            j = i + idx
            if j > len(lon):
                j = len(lon)
            
            if np.all(np.isnan(lon[i:j])) or np.all(np.isnan(lat[i:j])):
                continue
            lons.append(lon[i:j])
            lats.append(lat[i:j])
        return lons, lats
    
    def merge_files(self, num):
        from tqdm import trange
        time_dim = 0
        trajectory_dim = 0
        file_max = 0

        for i in range(num):
            file_name = self.filename + str(i) + ".nc"
            data = Dataset(file_name)
            time_dim = max(time_dim, len(data.dimensions["time"]))
            if time_dim == len(data.dimensions["time"]):
                file_max = i
            trajectory_dim += len(data.dimensions["trajectory"])
            if i != num -1:
                data.close()
        
        variables = {"trajectory": np.arange(1, trajectory_dim + 1, 1)}
        variable_properties = {}
        
        new_dataset = Dataset(self.filename + ".nc", "w")
        new_dataset.createDimension("time", time_dim)
        new_dataset.createDimension("trajectory", trajectory_dim)

        for variable in data.variables.keys():
            for i in trange(num, desc="Now writing: " + variable):
                file_name = self.filename + str(i) + ".nc"
                data = Dataset(file_name)
                value = data.variables[variable]
                if i == file_max:
                    variable_properties[variable] = {"dtype": value.dtype, "dim": value.dimensions, "other":value.__dict__}
                if variable == "time":
                    if len(value) == time_dim:
                        variables["time"] = value[:]
                
                elif variable == "trajectory":
                    continue

                else:
                    if len(value[0]) < time_dim:
                        missing_len = time_dim - len(value[0])
                        value = ma.append(value, ma.MaskedArray(np.zeros((len(value),missing_len)), mask = True), axis = 1)
                    
                    if variable not in variables.keys():
                        variables[variable] = value[:]
                        
                    else:
                        variables[variable] = ma.append(variables[variable], value, axis = 0)
                
                data.close()
            
            new_var = new_dataset.createVariable(variable, variable_properties[variable]["dtype"], variable_properties[variable]["dim"])
            new_var.setncatts(variable_properties[variable]["other"])
            new_var[:] = variables[variable]
            variable_properties.pop(variable)
            variables.pop(variable)

        for i in range(num):
            file_name = self.filename + str(i) + ".nc"
            os.remove(file_name)
        
        new_dataset.close()


        # for varname, var in variable_properties.items():
        #     new_dataset = Dataset(self.filename + "_" + varname + ".nc", "w")
        #     new_dataset.createDimension("time", time_dim)
        #     new_dataset.createDimension("trajectory", trajectory_dim)
        #     new_var = new_dataset.createVariable(varname, var["dtype"], var["dim"])
        #     new_var.setncatts(var["other"])

        #     new_var[:] = variables[varname]
        #     new_dataset.close()

