import os
import numpy as np
from datetime import datetime, timedelta
import sys

sys.path.append("/home/perharic/Documents/model")

from model.massdecay.divideandconquer import DivideAndConquerCarbonDrift
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_global_landmask

from global_land_mask import globe


temperature_file, current_file, outfile = sys.argv[1], sys.argv[2], sys.argv[3]

# os.chdir("/home/perharic/Documents")

bat = reader_netCDF_CF_generic.Reader('etopo2_cmems_med.nc', standard_name_mapping={"topo":"depth"})
# depthmax = np.abs(np.min(Dataset("etopo2_cmems_med.nc")["topo"][:]))
os.chdir("/home/perharic/Documents/2003/readers")

temp = reader_netCDF_CF_generic.Reader(temperature_file)
# time = temp.start_time
time = datetime(year = 2003, month=6, day = 12, hour = 12)
# curr = reader_netCDF_CF_generic.Reader(current_file)
reader_landmask = reader_global_landmask.Reader()

lons, lats, num = rectangle_seed(temperature_file, cut = 2)
lons, lats = mask_seed(lons, lats)

mass = np.ones(num)

os.chdir("/home/perharic/Documents/paper")

configure = {'drift:advection_scheme':"runge-kutta",
             'general:use_auto_landmask': False}

o = DivideAndConquerCarbonDrift(loglevel = 0, initial_velocity = -0.01,  distribution = "mass_sqrt", mass_threshold = 0.1, decay_type = "linear", configure = configure,
                                mass = mass, lon = lons, lat = lats, starttime = time, reader = [bat, temp, reader_landmask], max_num = num//2, deactivate_horizontal_advection = True)

# o.economic_run(steps = 200, time_step = timedelta(minutes = 30), outfile = outfile)

o.tracker = "jun_2003_f2_tracker.txt"
o.recursive_run_conquer()
o.delete_in_between_files()
# o.plot(linecolor="z", fast = True)
