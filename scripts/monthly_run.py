
#%%
import os
import numpy as np
from datetime import datetime, timedelta
import sys

sys.path.append("/home/perharic/Documents/model")

from carbondrift import *
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_global_landmask
from useful_functions import *

from global_land_mask import globe

#Read files
temperature_file, current_file, outfile = sys.argv[1], sys.argv[2], sys.argv[3]

#Import readers
bat = reader_netCDF_CF_generic.Reader('etopo2_cmems_med.nc', standard_name_mapping={"topo":"depth"})

#Calculate max depth for depth normalization in z dependent fragmentation distribution
depthmax = np.abs(np.min(Dataset("etopo2_cmems_med.nc")["topo"][:]))

# os.chdir("/home/perharic/Documents/2003/readers")
os.chdir("/home/perharic/Documents/climatology")

temp = reader_netCDF_CF_generic.Reader(temperature_file)#, standard_name_mapping={"thetao_avg":"sea_water_temperature"})
temp.verticalbuffer = 100
temp.always_valid = True
# time = temp.start_time
time = datetime(year = 2003, month=6, day = 12, hour = 12)

# curr = reader_netCDF_CF_generic.Reader(current_file)

reader_landmask = reader_global_landmask.Reader()

#Create seeding grid.
lons, lats, num = rectangle_seed(temperature_file, cut = 3)
#Mask in land elements.
lons, lats = mask_seed(lons, lats)

mass = np.ones(num)

os.chdir("/home/perharic/Documents/paper")


o = CarbonDrift(loglevel = 0, initial_velocity = -0.01,  distribution = "mass_sqrt", z_threshold = 0.1, depth_norm = depthmax, mass_threshold = 0.1, decay_type = "linear")
o.deactivate_fragmentation()
o.deactivate_horizontal_advection()

o.add_reader([bat, temp, reader_landmask])
o.set_config('drift:advection_scheme', 'runge-kutta')
o.set_config('general:use_auto_landmask', False)

# o.environment_profiles = {'sea_water_temperature': temp}
# o.required_profiles_z_range = [-5000, 0] 

o.seed_elements(lon=lons, lat = lats, z=0, mass = mass, time=time)
o.run(steps = 250, time_step=timedelta(minutes = 30), outfile=outfile)

# o.plot(linecolor="z", fast = True)
