import numpy as np
from copy import copy
import model.massdecay.gridrun as mcdgrid
import model.areadecay.gridrun as acdgrid
import model.massdecay.carbondrift as mcd
import model.areadecay.carbondrift as acd
from simulation.param_classifier import Parameters
from opendrift.readers import reader_global_landmask
from opendrift.readers import reader_netCDF_CF_generic
from simulation.seeding import *

from model.logger import Logger

log = Logger("CarbonDrift.simulation.prepare_run")
logger = log.LOGGER


class PrepareSimulation:

    def __init__(self, simulation_type, **kwargs):
        move_to_ocean = kwargs.pop("oceanonly")
        self.microbialdecaytype = kwargs.pop("microbialdecaytype")
        seeddata = kwargs.pop("seeddata")
        sf = kwargs.pop("splitfactor")
        logger.debug("Classifying parameters.")
        self.p = Parameters(kwargs)
        logger.debug("Defining initial position of particles.")
        self.seeder = SeedFromFile(seeddata)
        logger.debug("Preparing readers.")
        self.readers = self.prepare_readers(**self.p.readers)
        logger.debug("Setting config.")
        self.config = {'drift:advection_scheme':"runge-kutta", 'general:use_auto_landmask': False, 'seed:ocean_only':move_to_ocean}
        logger.debug("Creating CD object instance.")
        if simulation_type == "normal":
            self.initialize_normal_run()
        elif simulation_type == "grid":
            self.initialize_grid_run(sf)
        else:
            raise NameError(simulation_type + " has not yet been implemented.")
    
    def prepare_readers(self, temperature, bathymetry, readerbuffer, current = None):
        readers = []
        bat = reader_netCDF_CF_generic.Reader(bathymetry, standard_name_mapping={"topo":"depth"})
        readers.append(bat)
        tmp = reader_netCDF_CF_generic.Reader(temperature)
        # if self.microbialdecaytype == "area":
        #     sea_surface_temperature = tmp.get_variables_interpolated(["sea_water_temperature"], "sea_water_temperature", -1, time = self.p.object_init["starttime"],
        #                                         lon = self.seeder.lon, lat = self.seeder.lat, z = -0.5)[0]["sea_water_temperature"].data
        #     self.p.object_init["sea_surface_temperature"] = sea_surface_temperature
        tmp.verticalbuffer = readerbuffer
        tmp.always_valid = True
        readers.append(tmp)
        if current is not None:
            curr = reader_netCDF_CF_generic.Reader(current)
            curr.verticalbuffer = readerbuffer
            # curr.always_valid = True
            readers.append(curr)
        readers.append(reader_global_landmask.Reader())
        return readers
    
    def initialize_grid_run(self, sf):
        self.p.object_init["reader"] = self.readers
        self.p.object_init["configure"] = self.config
        self.p.object_init["lon"] = self.seeder.lon
        self.p.object_init["lat"] = self.seeder.lat
        self.p.object_init["m0"] = copy(self.seeder.mass)
        self.p.object_init["z"] = self.seeder.z
        self.p.object_init["split_factor"] = sf
        if hasattr(self.seeder, "origin_marker"):
            self.p.object_init["origin_marker"] = self.seeder.origin_marker
        if self.microbialdecaytype == "mass":
            self.obj = mcdgrid.GridRun(**self.p.object_init)
        elif self.microbialdecaytype == "area":
            # self.p.object_init["r0"] = self.seeder.r0
            self.obj = acdgrid.GridRun(**self.p.object_init)
    
    def initialize_normal_run(self):
        params = self.p
        advection = params.object_init.pop("deactivate_horizontal_advection")
        fragmentation = params.object_init.pop("deactivate_fragmentation")
        time = params.object_init.pop("starttime")
        params.object_init["m0"] = copy(self.seeder.mass)
        if self.microbialdecaytype == "mass":
            self.obj = mcd.CarbonDrift(**params.object_init)
        elif self.microbialdecaytype == "area":
            # params.object_init["r0"] = self.seeder.r0
            self.obj = acd.CarbonDrift(**params.object_init)

        if advection:
            self.obj.deactivate_horizontal_advection()
        if fragmentation:
            self.obj.deactivate_fragmentation()
        self.obj.add_reader(self.readers)

        for key, value in self.config.items():
            self.obj.set_config(key, value)
        
        if hasattr(self.seeder, "origin_marker"):
                origin_marker = self.seeder.origin_marker
        else:
                origin_marker = 0
        if self.microbialdecaytype == "mass":
            self.obj.seed_elements(lon=self.seeder.lon, lat = self.seeder.lat, z=self.seeder.z,
                                   mass = self.seeder.mass, time=time, origin_marker = origin_marker)
        elif self.microbialdecaytype == "area":
            self.obj.seed_elements(lon=self.seeder.lon, lat = self.seeder.lat, z=self.seeder.z,
                                   mass = self.seeder.mass, time=time, origin_marker = origin_marker)