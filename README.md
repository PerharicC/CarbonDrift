# CarbonDrift
Lagrangian 3d  tracking of carbon particles

## Download

Open the terminal and clone CarbonDrift.

```console
~$ git clone <path>
```

## Run a simulation

Run the following comand in the root directory

```console
~/CarbonDrift$ python -m simulation.run -tmp <temperature netcdf filepath> -b <bathimetry netcdf filepath> -s <starttime> -o <out netcdf filepath> 
```

For further information on other arguments use header in terminal.

**Note:** To avoid long terminal inputs, one can save paramers, which rarely change, to a .txt file of format

```
-tmp
<temperatuure netcdf filepath>
-b
<bathimetry netcdf filepath>
-s
2010-01-01-0
```

The simulation can then be run as:

```console
~/CarbonDrift$ python -m simulation.run @parameters.txt -o <out netcdf filepath> 
```
where one can now add the changing parameters next to the .txt file.

## Plotting a simulation

Now that we have successfully ran and saved a simulation to a netCDF file, we can look at how to plot the results. In the root directory type in the terminal

```console
~/CarbonDrift$ python -m plotting.plot_run <plot method, e.g mass_map> -f1 <file1 location, e.g. clim.nc> -f2 <file2 location, e.g. MHW.nc (not necessary if not plotting difference)> -lons <lon grid txt filepath (created in simulation)> -lats <lat grid txt filepath (created in simulation)> -o <outfile.pdf filepath> -d <depth at which mass is summed> -clip <min:max (use m for minus)>
```

Similarly one can save ceertain/all paramaers to a .txt file and in a terminal simply run the command

```console
~/CarbonDrift$ python -m plotting.plot_run @parameters.txt
```

### Plot locations of specified drifters

```console
~/CarbonDrift$ python -m plotting.plot_run drifter_locations -f1 <file1 location, e.g. clim.nc> -loc <locations in format lon1:lat1,lon2:lat2 etc> -o <outfile.pdf filepath> -lines
```

Add argument lines only if you would like to display gridlines to the drifter locations.

**Example:**

```console
~/CarbonDrift$ python -m plotting.plot_run drifter_locations -f1 <file1.nc> -loc m135:0,20:m45,90:m20 -lines
```

**Output:**

![](/images/locations.png)

### Plot some drifter properties

To make a 2D plot of some property vs another property fo two different simulations type in terminal:

```console
~/CarbonDrift$ python -m plotting.plot_run drifter_properties -f1 <file1 location, e.g. clim.nc> -f2<file2 location, e.g. MHW.nc> -loc <locations in format lon1:lat1,lon2:lat2 etc> -o <outfile.pdf filepath> -p1 <property on x axis> -p2 <property on y axis>
```

One can choose between properties


```python
["time", "mass", "z", "temperature", "x_sea_water_velocity", "y_sea_water_velocity"]
```

**Example1:**

```console
~/CarbonDrift$ python -m plotting.plot_run drifter_properties -f1 <file1.nc> -f2<file2.nc> -loc m135:0,20:m45,90:m20 -p1 time -p2 mass
```

**Output1:**

![](/images/mass_t.png)

**Example2:**

```console
~/CarbonDrift$ python -m plotting.plot_run drifter_properties -f1 <file1.nc> -f2<file2.nc> -loc m135:0,20:m45,90:m20 -p1 mass -p2 z
```

**Output2:**

![](/images/z_mass.png)
