import sys
import os
import numpy as np
sys.path.append("/home/perharic/Documents/model")

from plots import Plot

os.chdir("/home/perharic/Documents/paper")
p = Plot(sys.argv[1], distribution="mass_sqrt", fig_size=(20, 10))
os.chdir("/home/perharic/Documents")

lons = np.loadtxt("lons.txt")
lats = np.loadtxt("lats.txt")

os.chdir("/home/perharic/Documents/paper/images")

if len(sys.argv) >= 3:
    title = " ".join(sys.argv[2].split("-"))
else:
    title = None
if len(sys.argv) > 3:
    file2 = sys.argv[3]
else:
    file2 = None
if len(sys.argv) > 4:
    temp_file = sys.argv[4]
else:
    temp_file = None

# file = "mass.txt"
# p.mass_map(lons, lats, title, file2 = file2, pathtofile2= "/home/perharic/Documents/paper")
# p.reset_figure()
# file = ["mass1.txt", "mass2.txt"]
# p.zone_mass_distribution(lons, lats, title, file2, pathtofile2="/home/perharic/Documents/frag_diff/no_advection", read_file = file, labels = ["2000", "2003"], m0 = 31101)
# p.animate_3D(k = 10, outfile="jul.gif")#, markevery=20, plot_floor=False)
# p.reset_figure()
# p.vertical_mass_distribution(bins = 100, file2 = file2, title = title)
# p.mass_z(file2)