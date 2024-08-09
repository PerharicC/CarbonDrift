import numpy as np
from model.massdecay.gridrun import GridRun
from model.massdecay.carbondrift import *
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
        logger.debug("Classifying parameters.")
        self.p = Parameters(kwargs)
        logger.debug("Preparing readers.")
        self.readers = self.prepare_readers(**self.p.readers)
        logger.debug("Setting config.")
        self.config = {'drift:advection_scheme':"runge-kutta", 'general:use_auto_landmask': False, 'seed:ocean_only':move_to_ocean}
        logger.debug("Defining initial position of particles.")
        self.seeder = Seed(**self.p.seeding)
        logger.debug("Creating CD object instance.")
        if simulation_type == "normal":
            self.initialize_normal_run()
        elif simulation_type == "grid":
            self.initialize_grid_run()
        else:
            raise NameError(simulation_type + " has not yet been implemented.")
    
    def prepare_readers(self, temperature, bathymetry, temperaturebuffer, current = None):
        readers = []
        bat = reader_netCDF_CF_generic.Reader(bathymetry, standard_name_mapping={"topo":"depth"})
        readers.append(bat)
        tmp = reader_netCDF_CF_generic.Reader(temperature)
        tmp.verticalbuffer = temperaturebuffer
        tmp.always_valid = True
        readers.append(tmp)
        if current is not None:
            curr = reader_netCDF_CF_generic.Reader(current)
            readers.append(curr)
        readers.append(reader_global_landmask.Reader())
        return readers
    
    def initialize_grid_run(self):
        self.p.object_init["reader"] = self.readers
        self.p.object_init["configure"] = self.config
        self.p.object_init["lon"] = self.seeder.lon
        self.p.object_init["lat"] = self.seeder.lat
        self.obj = GridRun(**self.p.object_init)
    
    def initialize_normal_run(self):
        params = self.p
        advection = params.object_init.pop("deactivate_horizontal_advection")
        fragmentation = params.object_init.pop("deactivate_fragmentation")
        time = params.object_init.pop("starttime")
        self.obj = CarbonDrift(**params.object_init)

        if advection:
            self.obj.deactivate_horizontal_advection()
        if fragmentation:
            self.obj.deactivate_fragmentation()
        self.obj.add_reader(self.readers)

        for key, value in self.config.items():
            self.obj.set_config(key, value)
        
        self.obj.seed_elements(lon=self.seeder.lon, lat = self.seeder.lat, z=0, mass = np.ones(len(self.seeder.lon)), time=time)