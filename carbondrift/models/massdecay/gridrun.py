import os

from carbondrift.models.massdecay.carbondrift import *
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
        self.time = kwargs.pop("starttime")
        self.reader = kwargs.pop("reader", [])
        self.lon = kwargs.pop("lon")
        self.lat = kwargs.pop("lat")
        self.z = kwargs.pop("z")
        self.m0 = kwargs["m0"]
        self.origin_marker = kwargs.pop("origin_marker", None)
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
            o.seed_elements(lon = self.lon, lat = self.lat, z = self.z,
                            time = self.time, mass = self.m0, origin_marker = self.origin_marker)

            o.run(**kwargs)
        else:    
            lons, lats, z, mass, origin_marker = self.split()
            self.filename = copy(kwargs["outfile"]).split(".")[0]

            for i in range(len(lons)):
                self.o_kwargs["m0"] = mass[i]
                o = CarbonDrift(**self.o_kwargs)

                if self.reader is not None:
                    o.add_reader(self.reader)
                
                if self.deactivate_frag:
                    o.deactivate_fragmentation()
                if self.deactivate_ha:
                    o.deactivate_horizontal_advection()
                
                for key, value in self.config.items():
                    o.set_config(key, value)

                o.seed_elements(lon = lons[i], lat = lats[i], z = z[i],
                                time = self.time, mass = mass[i], origin_marker = origin_marker[i])
                kwargs["outfile"] = self.filename + str(i) + ".nc"
                o.run(**kwargs)
                del o
            
            if self.merge:
                self.merge_files(len(lons))

    def split(self):
        array_sizes = [len(self.lon) // self.split_factor] * self.split_factor
        for i in range(len(self.lon) % self.split_factor):
            array_sizes[i] += 1
        
        indices = np.cumsum(array_sizes[:-1])
        lons = np.split(self.lon, indices)
        lats = np.split(self.lat, indices)
        z = np.split(self.z, indices)
        mass = np.split(self.m0, indices)
        if self.origin_marker is not None:
            origin_marker = np.split(self.origin_marker, indices)
        else:
            origin_marker = None
        return lons, lats, z, mass, origin_marker
    
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

