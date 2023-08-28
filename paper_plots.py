import sys
import os
import numpy as np
sys.path.append("/home/perharic/Documents/model")

from plots_cleaned import Plot

os.chdir("/home/perharic/Documents/paper")
locations = [(5.45, 42.02), (10.79, 40.02)]
locations = [(10.33, 41.56), (7.5, 43.2)]
locations = [(9.96, 39.27), (9.75, 41.8125)]
p = Plot(sys.argv[1], distribution="mass_sqrt", fig_size=(15, 10), locations = locations)
os.chdir("/home/perharic/Documents")

lons = np.loadtxt("lons.txt")
lats = np.loadtxt("lats.txt")
os.chdir("/home/perharic/Documents/paper")

file2 = sys.argv[2]

# p.area(lons, lats)

# p.figsize = (12, 9)
# p.mass_z(file2)

# p.plot_3D(k = 10, markevery = 20)
# p.figsize = (14, 7)
# p.time_z(file2)

# p.figsize = (10, 15)
# p.mass_map(lons, lats, file2, pathtofile2 = "/home/perharic/Documents/paper", read_file = "mass_exp.txt")

# p.zone_mass_distribution(lons, lats, file2 = file2, pathtofile2 = "/home/perharic/Documents/paper")

# p.figsize = (15, 10) 
# p.m_T_correlations(markevery = 15)

p.mass_at_depth_time(file2)
