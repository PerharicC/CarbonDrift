# CarbonDrift
Lagrangian 3d tracking of carbon particles

## Download

Open the terminal and clone CarbonDrift.

```console
~$ git clone <path>
```

## Run a simulation
### Standard opendrift simulation
To run a simulation one can follow the standard procedureof any opendrift submodule:
1. Create an object instance.
2. Add readers, configures etc. 
3. Seed elements at desired locations, depths, with appropriate masses etc.
4. Run and save the simulation by calling the run() method.

### Automized method
To run a CarbonDrift simulation from the terminal one can simply run the following comand in the root directory

```console
~/CarbonDrift$ python -m simulation.run -tmp <temperature netcdf filepath> -b <bathimetry netcdf filepath> -s <starttime> -o <out netcdf filepath> 
```

For further information on other arguments use header in terminal.

**Note:** To avoid long terminal inputs, one can save parameters to a parameters.txt file of format

```
-tmp
<temperatuure netcdf filepath>
-b
<bathimetry netcdf filepath>
-s
2010-01-01-0
--steps
420
```

The simulation can then be run as:

```console
~/CarbonDrift$ python -m simulation.run @parameters.txt -o <out netcdf filepath> 
```


## Plotting a simulation

CarbonDrift provides some built in plotting methods in the model.plots submodule. To run any one of these methods open the terminal in the root directory and type

```console
~/CarbonDrift$ python -m plotting.plot_run <plot method, e.g mass_map> -f1 <file1path> -f2 <file2path (not necessary if not plotting difference)> -o <outfile.pdf filepath> -d <depth at which mass is summed> -clip <min:max (use m for minus)>
```

Similarly one can save certain/all parameters to a .txt file and in a terminal simply run the command

```console
~/CarbonDrift$ python -m plotting.plot_run @parameters.txt
```
### Plotting examples
In this section we will show some plotting methods and their outputs.
1. **Plot locations of specific drifters**

Save the following parameters to a .txt file:

```
drifter_locations
-f1
<simualtion filepath>
-loc
m48:30,5:0
-t
Locations_of_drifters
-legend
$L_1$,$L_2$
-o
./images/locations.png
```

Next run the following command in the terminal

```console
~/CarbonDrift$ python -m plotting.plot_run @location_params.txt
```
And the output will be:

![](/images/locations.png)

2. **Plot some drifter properties**

To make a 2D plot of some property vs another property for up to four different simulations the parameter file should contain the following arguments:

```
drifter_properties
-f1
<file1 path>
-f2
<file2 path>
-f3
<file3 path>
-f4
<file4 path>
-fs
8:10
-loc
m48:30
-p1
mass
-p2
z
-lines
-t
LIN_decay,_z0_=_-100_m,_w0_=_500_m/d
-xlabel
$m/m_0$
-ylabel
$z\,\mathrm{[m]}$
-lw
3
-xlim
0:1
-ylim
m4500:0
-legend
$\frac{\mathrm{d}m}{\mathrm{d}t}=-km\text{  - constant speed}$,$\frac{\mathrm{d}m}{\mathrm{d}t}=-kS\text{  - constant speed}$,$\frac{\mathrm{d}m}{\mathrm{d}t}=-km\text{  - variable speed}$,$\frac{\mathrm{d}m}{\mathrm{d}t}=-kS\text{  - variable speed}$ 
-o
./images/z_mass.png
```

The output should look something like this:

![](/images/z_mass.png)

To plot other properties one can choose between

```python
["time", "mass", "z", "sea_water_temperature", "x_sea_water_velocity", "y_sea_water_velocity"]
```

3. **Cartopy map of mass at given depth**

To plot the mass of each grid cell at a given depth the following parameters may be provided:

```
mass_map
-f1
<file1 path>
-fs
12:8
-cblabel
$m/m_0\,\,\text{at sea floor}$
-t
EXP_decay:_z0_=_-100_m,_w0_=_500_m/d
-d
-6000
-cmap
tab20b
--shrink
0.8
-o
./images/mass_map.png
```

**NOTE:** If one adds the -abs argument, the current mass will be plotted, otherwise the mass fraction $\frac{m}{m_0}(\varphi, \lambda)$ will be plotted for each grid cell.

**NOTE:** If one would like to plot the difference between two simulations, file2 should also be provided in the parameters, together with the -diff argument. One can then impose the -abs argument to plot the absolute difference file1 - file2. If -abs is not added, the plotted result will be (file1 - file2) / file1.

**NOTE:** If the depth argument -d is set to less than -5000, depth will be set to sea floor.

A possible output will look something like:

![](/images/mass_map.png)

4. **Total mass at depth**

To calculate the mass sum of all grid cells at a given depth, one can run the simulation with parameters

```
get_mass_sum_at_depth
-f1
<file1 path>
-d
-100
-abs
```

The -abs argument insures, that the absolute sum is returned. Alternatively the fraction
$\frac{\sum_i m_i (z = d)}{\sum_i m_i (t = 0)}$ is plotted.

To get biome specific masses the method *get_mass_sum_at_depth* should be replaced by *get_biome_weighted_mass_at_depth* in the parameter.txt file.

**WORK IN PROGRESS**
