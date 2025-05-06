import numpy as np
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_global_landmask
from datetime import datetime, timedelta
from netCDF4 import Dataset
import os

import carbondrift
from carbondrift.models.areadecay import carbondrift as cdrift

def run(outfile):
    #Import readers
    sup_dir = f"{os.path.dirname(os.path.dirname(carbondrift.__file__))}/supplementary_data"
    btm = reader_netCDF_CF_generic.Reader(os.path.join(sup_dir, "etopo2_remaped1deg.nc"), standard_name_mapping={"topo":"depth"})
    tmp = reader_netCDF_CF_generic.Reader(os.path.join(sup_dir, "tmp_luo_remaped1deg.nc"))
    lm = reader_global_landmask.Reader()
    #Increase buffer, since svertical speeds are large!
    tmp.verticalbuffer = 100
    tmp.always_valid = True
    #Define start time
    time = datetime(year=1993, month = 1, day = 1)
    #Set up configures
    config = {'drift:advection_scheme':"runge-kutta", 'general:use_auto_landmask': False, 'seed:ocean_only':False}
    #Initial mass array
    mass =  np.ones(2)
    #Initialize Carbon Drift object with area decay
    o = cdrift.CarbonDrift(loglevel = 0, initial_velocity = -0.009259, m0 = mass, decay_type = "linear")
    #Deactivate horizontal advection and fragmentation
    o.deactivate_horizontal_advection()
    o.deactivate_fragmentation()
    #Set readers and configures
    o.add_reader([btm, tmp, lm])
    for key, value in config.items():
        o.set_config(key, value)
    #Seed elements at specified locations
    o.seed_elements(lon = [146, -45], lat = [58, 30], time = time, mass = mass, z = -50)
    #Run the simulation
    o.run(outfile=outfile, steps = 5, time_step=timedelta(minutes = 30), time_step_output=timedelta(hours = 2))

if __name__ == "__main__":
    outfile = "test.nc"
    run(outfile)
    #Inspect file
    data = Dataset(outfile)
    data.close()
