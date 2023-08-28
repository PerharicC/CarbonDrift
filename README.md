# CarbonDrift
Lagrangian 3d  tracking of carbon particles

## Run a simulation

Start by importing all the necessary libraries/modules.

```python
import numpy as np
from datetime import datetime, timeldelta

from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_global_landmask

from runbysubgrids_split_vars import GridRun
from useful_functions import *
```

Next we define the file locations for the readers and the outfile filename.

```python
temperature_file = f"path_to_temperature_netCDF_file"
bat_file = f"path_to_bathymetry_netCDF_file"
outfile = "netCDF_outfile_filename"
```
We then initialise the opendrift Reader objects.

```python
tmp = reader_netCDF_CF_generic.Reader(temperature_file)
bat = reader_netCDF_CF_generic.Reader(bat_file)
reader_landmask = reader_global_landmask.Reader()
```

This next step only applies to climatology temperature files, which have an unlimited time dimension. Thus we have to tell opendrift, that it should read from this file for all times.

```python
tmp.always_valid = True #Default False
```

We now ensure that the buffer is large enough, so that the environment variables are corroectly updated. FUrthermore we define the simulation starting time.


```python
tmp.verticalbuffer = 100 #Deault 1
time = datetime(year = year, month = motnh, day = day, hour = hour)
```

Next we must define the lon & lat values for seeding the elements. This can be easily achieved with two simple 1D numpy arrays. For a grid seeding, however, one should write the following code

```python
#Get a meshgrid array for lons and lats corresponding to the grid of the given netCDF fie.
#Num is the number of cells, i.e. num = lons.shape[0] * lons.shape[1]
lons, lats, num = rectangle_seed(temperature_file)
#Mask out cells which are on land, by filling them with nans.
lons, lats = mask_seed(lons, lats)
```

If this is the first time running a simulation with this grid, you should save the 1D grid arrays by calling
```python
lons, lats, num = rectangle_seed(temperature_file, save = True)
```
Furthermore if you would like to reduce the number of cells, you can add a key argument **cut** to the rectangle_seed function. For example a cut value of 3 will reduce the number of rows and columns by three, i.e. num&rarr;num // 9.

Next we make a dictionary of configures to be set and a list of all readers to be read.

```python
configure = {'drift:advection_scheme':"runge-kutta",
             'general:use_auto_landmask': False}
readers = [bat, tmp, reader_landmask]
```

We can now finally initialise our GridRun object

```python
loglevel = 0 #Minimal information display in terminal.

split_factor = 1#Number of seperate simulations.
#If 1, no splitting will occur and the run will be equivalent
#to simply using the carbondrift.CarbonDrift module directly.

decay_type = "linear" #Or exponentital.

o = GridRun(loglevel = loglevel, initial_velocity = -0.01, decay_type = decay_type,
            deactivate_fragmentation = True, starttime = time, reader = readers,
            configure = configure, split_factor = split_factor, lon = lons, lat = lats)
```
Finally we run the simulation by calling the run method.

```python
simulation_steps = 200
minutes = 30
time_step = timedelta(minutes = minutes)
o.run(steps = simulation_steps, time_step = time_step, outfile = outfile)
```

## Plotting a simulation
Now that we have successfully ran and saved a simulation to a netCDF file, we can look at how to plot the results. We do this through the Plot object in plots_cleaned. We begin again with importing all the necessary libraries.

```python
import numpy as np
from plots_cleaned import Plot
```

Most of the methods in Plot are written for comparing two different netCDF files, for example a simulation of carbon decay using the climatology data
and a simulation using data from a marine heatwave. We thus define the two file location paths.

```python
clim_data = "path_to_climatology_simulation_netCDF_file"
mhw_file = "path_to_marine_heatwave_simulation_netCDF_file"
```

For some plots, properties at two specific locations are plotted for comparison. Therefore we must define the coordinates of the the locations as a list of tuples.

```python
locations = [(x1, y1), (x2, y2)]
```
