import os
import numpy as np
from datetime import datetime, timedelta
import sys

sys.path.append("/home/perharic/Documents/model")

from runbysubgrids_split_vars import GridRun
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_global_landmask
from useful_functions import *


temperature_file, current_file, outfile = sys.argv[1], sys.argv[2], sys.argv[3]

# os.chdir("/home/perharic/Documents")

bat = reader_netCDF_CF_generic.Reader('etopo2_cmems_med.nc', standard_name_mapping={"topo":"depth"})
depthmax = np.abs(np.min(Dataset("etopo2_cmems_med.nc")["topo"][:]))

os.chdir("/home/perharic/Documents/2003/readers")
os.chdir("/home/perharic/Documents/climatology")

temp = reader_netCDF_CF_generic.Reader(temperature_file)
temp.verticalbuffer = 100
temp.always_valid = True
# time = temp.start_time
time = datetime(year = 2003, month=6, day = 12, hour = 12)

# curr = reader_netCDF_CF_generic.Reader(current_file)

reader_landmask = reader_global_landmask.Reader()

lons, lats, num = rectangle_seed(temperature_file, cut = 3)
lons, lats = mask_seed(lons, lats)

os.chdir("/home/perharic/Documents/paper")

configure = {'drift:advection_scheme':"runge-kutta",
             'general:use_auto_landmask': False}

o = GridRun(loglevel = 0, initial_velocity = -0.01,  distribution = "mass_sqrt", mass_threshold = 0.1,
            decay_type = "linear", deactivate_horizontal_advection = True, deactivate_fragmentation = True, starttime = time, reader = [bat, temp, reader_landmask],
            configure = configure, split_factor = 1, lon = lons, lat = lats)


o.run(steps = 200, time_step=timedelta(minutes = 30), outfile=outfile)
