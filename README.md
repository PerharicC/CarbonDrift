# CarbonDrift
Lagrangian 3d  tracking of carbon particles

## Run a simulation

Run the following comand in the root directory

```console
~/CarbonDrift$ python -m simulation.run -tmp <temperature netcdf filepath> -b <bathimetry netcdf filepath> -s <starttime> -o <out netcdf filepath> 
```

For further information on other arguments use header in terminal.

## Plotting a simulation

Now that we have successfully ran and saved a simulation to a netCDF file, we can look at how to plot the results. In the root directory type in the terminal

```console
~/CarbonDrift$ python -m plotting.plot_run <plot method, e.g mass_map> -f1 <file1 location, e.g. clim.nc> -f2 <file2 location, e.g. MHW.nc (not necessary if not plotting difference)> -lons <lon grid txt filepath (created in simulation)> -lats <lat grid txt filepath (created in simulation)> -o <outfile.pdf filepath> -d <depth at which mass is summed> -clip <min:max (use m for minus)>
```
# Old not yet updated...
### Plot the area of where the elements were initially seeded

To plot the area of the grid where the elements were initially seeded, together with the two marked locations, we call the method

```python
p.area(lons, lats)
```

![](/images/fig1_loc2.png)


### Plot the z to fraction of mass dependance

To plot $z(\Delta m/m_0)$ at the two specified locations, we call the method

```python
p.mass_z(mhw_file)
```

![](/images/fig3_dm_m0_exp_loc2.png)


### Plot the z to time dependance

To plot $z(t)$ at the two specified locations, we call the method

```python
p.time_z(mhw_file)
```

![](/images/fig4_exp_loc2.png)


### Plot the fraction of mass, which crossed a certain sea depth

To plot the fraction of mass, which had crossed the euphotic and twilight zones or had reached the sea floor we call the method

```python
p.mass_map(lons, lats, mhw_file, write_file = "some_text.txt")
```

We use the *write_file* argument to store the calculations into a .txt file. The next time we run the code we simply change the argument *write_file*&rarr;*read_file*.

Example of an output can be found [here](/images/fig5_dm_m0_exp.pdf)




