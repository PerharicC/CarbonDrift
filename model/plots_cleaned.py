import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, CheckButtons
from matplotlib import animation
from matplotlib.gridspec import GridSpec

import cartopy
from cartopy import config
import cartopy.crs as ccrs

from carbondrift import *
import carbondriftopen as cdo
from netCDF4 import Dataset, num2date

from numba import jit
from datetime import datetime, timedelta
import cftime


# plt.rcParams.update({
#     "text.usetex": True,
#     "font.family": "Helvetica"
# })


class MidpointNormalize(mpl.colors.Normalize):
    def __init__(self, vmin, vmax, midpoint=0, clip=False):
        self.midpoint = midpoint
        mpl.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):

        if self.midpoint == self.vmax:
            normalized_min = 0
        else:
            normalized_min = max(0, 1 / 2 * (1 - abs((self.midpoint - self.vmin) / (self.midpoint - self.vmax))))
        
        if self.midpoint == self.vmin:
            normalized_max = 1
        else:
            normalized_max = min(1, 1 / 2 * (1 + abs((self.vmax - self.midpoint) / (self.midpoint - self.vmin))))
        
        normalized_mid = 0.5
        x, y = [self.vmin, self.midpoint, self.vmax], [normalized_min, normalized_mid, normalized_max]
        return np.ma.masked_array(np.interp(value, x, y))


class Open:
    """The alternative for opening .nc files, as I have found that opendrift sometimes has problems with understanding masked datasets."""

    def __init__(self, filename):

        """Parameters
        -------------
        filename: str
            Name of .nc file to be imported."""
        
        self.data = Dataset(filename)
        self.steps_output = len(self.data["time"][:])
    
    def get_property(self, property):
        return np.ma.transpose(self.data[property][:])
    
    def get_time_array(self):
        time = self.data.variables["time"]
        units = time.units
        return num2date(time[:], units)
    

class Plot:

    "Class for plotting saved .nc CarbonDrift simulations."

    def __init__(self, filename, distribution, locations = None, fig_size = (20, 20), fontsize = 17):

        """
        Parameters
        ----------
        filename: str
            .nc file to be ploted
        distribution: str
            Name of distribution used in .nc file
        fig_size: tuple of two ints
            Size of figure (default: (20, 20))
        fontsize: int
            Fontsize (default: 20)"""
        
        self.filename = filename
        # self.o = cdo.open(filename, distribution = distribution)
        self.o = Open(filename)
        self.time = self.o.get_time_array()
        self.distribution = distribution

        plt.rcParams.update({'font.size': fontsize})

        fig, ax = plt.subplots(1, 1, figsize = fig_size)
        self.fig = fig
        self.ax = ax
        self.figsize = fig_size
        self.loc = locations

        self.steps_output = self.o.steps_output
    
    @staticmethod
    def get_cmap(x, map):
        cmap = plt.get_cmap(map, len(x))
        c = np.arange(1., 1 + len(x))
        norm = mpl.colors.BoundaryNorm(np.arange(len(c)+1)+0.5,len(c))
        sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        return cmap, sm, c
    
    @staticmethod
    def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
        '''
        Function to offset the "center" of a colormap. Useful for
        data with a negative min and positive max and you want the
        middle of the colormap's dynamic range to be at zero.

        Input
        -----
        cmap : The matplotlib colormap to be altered
        start : Offset from lowest point in the colormap's range.
            Defaults to 0.0 (no lower offset). Should be between
            0.0 and `midpoint`.
        midpoint : The new center of the colormap. Defaults to 
            0.5 (no shift). Should be between 0.0 and 1.0. In
            general, this should be  1 - vmax / (vmax + abs(vmin))
            For example if your data range from -15.0 to +5.0 and
            you want the center of the colormap at 0.0, `midpoint`
            should be set to  1 - 5/(5 + 15)) or 0.75
        stop : Offset from highest point in the colormap's range.
            Defaults to 1.0 (no upper offset). Should be between
            `midpoint` and 1.0.
        '''
        cdict = {
            'red': [],
            'green': [],
            'blue': [],
            'alpha': []
        }

        # regular index to compute the colors
        reg_index = np.linspace(start, stop, 257)

        # shifted index to match the data
        shift_index = np.hstack([
            np.linspace(0.0, midpoint, 128, endpoint=False), 
            np.linspace(midpoint, 1.0, 129, endpoint=True)
        ])

        for ri, si in zip(reg_index, shift_index):
            r, g, b, a = cmap(ri)

            cdict['red'].append((si, r, r))
            cdict['green'].append((si, g, g))
            cdict['blue'].append((si, b, b))
            cdict['alpha'].append((si, a, a))

        newcmap = mpl.colors.LinearSegmentedColormap(name, cdict)
        plt.register_cmap(cmap=newcmap)

        return newcmap

    @staticmethod
    def get_colormap_midpoint(arr):
        nans = np.invert(np.logical_or(np.isnan(arr), np.isinf(arr)))
        m = np.min(arr[nans])
        if m > 0:
            m = 0
        M = np.max(arr[nans])
        return m, 1 - M / (M +np.abs(m)), M
    
    def reset_figure(self):
        fig, ax = plt.subplots(1, 1, figsize = self.figsize)
        self.fig = fig
        self.ax = ax
    
    def vertical_mass_distribution(self, maxdepth=None, bins=None, maxnum=None, file2 = None, log = True, title = None):
        """
        Plot the mass sum in each specified interval of z.
        
        Parameters
        ----------
        maxdepth: float
            Maxdepth of plotting. If None, the maxdepth reached by the particle is used.
        bins: int
            Number of bins for plotting. If None, the number is calculated from length of interval, 
            which, if not configured, will be set to 1 m.
        maxnum: float
            The high limit of bars. If None, it is calculated as the sum of the initial masses.
        file2: str
            File name of as second simulation, where the fragmentation was deactivated.
        """

        plt.close()
        fig = plt.figure(figsize = self.figsize)

        mainplot = fig.add_axes([.15, .3, .8, .5])
        sliderax = fig.add_axes([.15, .08, .75, .05])

        tslider = Slider(sliderax, 'Timestep', 0, self.steps_output-1,
                         valinit=0, valfmt='%0.0f', valstep = 1, color = "black")
        if title is None:
            title = "Vertical mass distribution"
        
        if file2 is not None:
            add2ax = fig.add_axes([0.7, 0.85, 0.2, 0.03])
            check = CheckButtons(
                                ax=add2ax,
                                labels=['Show file2'],
                                label_props={'color': "black"},
                                frame_props={'edgecolor': "black"},
                                check_props={'facecolor': "red"},)

        if file2 is not None:
            o2 = Open(file2)

            z2 = o2.get_property('z')
            z2 = np.ma.filled(z2, np.nan)

            mass2 = o2.get_property("mass")
            mass2 = np.ma.filled(mass2, 0)

            lon2 = self.o.get_property("lon")
            lon2 = np.ma.filled(lon2, np.nan)

            self.show2 = False
        else:
            z2 = [0]

        z = self.o.get_property('z')
        z = np.ma.filled(z, np.nan)

        mass = self.o.get_property("mass")
        mass = np.ma.filled(mass, 0)

        lon = self.o.get_property("lon")
        lon = np.ma.filled(lon, np.nan)

        if maxdepth is not None:
            maxrange = -np.abs(maxdepth)
        else:
            maxrange = min(np.min(z[np.invert(np.isnan(z))]), np.min(z2[np.invert(np.isnan(z2))]))
        
        if bins is not None:
            dz = -maxrange/bins
        else:
            try:
                dz = self.o.get_config('vertical_mixing:verticalresolution')
            except:
                dz = 1.
        
        @jit(nopython = True)
        def mass_dist(mass, z, m0, dz=dz, maxrange = maxrange):
            dist = []
            zmax = 0
            zmin = -dz
            while zmin > maxrange:
                dist.append(np.sum(mass[np.where((z > zmin) & (z<=zmax) & np.invert(np.isnan(z))& np.invert(np.isnan(mass)))[0]]) / m0)
                zmax = zmin
                zmin -= dz
            dist.append(np.sum(mass[np.where((z >= zmin) & (z<zmax)& np.invert(np.isnan(z))& np.invert(np.isnan(mass)))[0]]) / m0)
            return dist[::-1]
        
        @jit(nopython = True)
        def mass_in_zone(mass, z, m0):
            mass200 = np.sum(mass[np.where(z >= -200& np.invert(np.isnan(z)))[0]]) / m0
            mass1000 = np.sum(mass[np.where((z < -200) & (z >= -1000)& np.invert(np.isnan(z)))[0]]) / m0
            massdepth = np.sum(mass[np.where(z < -1000& np.invert(np.isnan(z)))[0]]) / m0
            return [mass200, mass1000, massdepth]
        
        @jit(nopython=True)
        def precalculate(mass, z, lon):
            dists = []
            mass_zone = []
            at_sea = np.invert(np.isnan(lon))
            for timestep in range(mass.shape[0]):
                mi = mass[timestep, at_sea]
                zi = z[timestep, at_sea]
                if timestep > 0:
                    idx = idx = np.where((mi < 1) &np.invert(np.isnan(mi)) &(zi <0))[0]
                    mi = mi[idx]
                    zi = zi[idx]
                else:
                    m0 = np.sum(mi)
                dists.append(mass_dist(mi, zi, m0))
                mass_zone.append(mass_in_zone(mi, zi, m0))
            return dists, mass_zone
        
        print("Precalculating values to plot for faster transitions.")
        dists, mass_zone = precalculate(mass, z, lon[0, :])
        print("Finsihed precalculating.")

        if file2 is not None:
            print("Precalculating second file values to plot for faster transitions.")
            dists2, mass_zone2 = precalculate(mass2, z2, lon2[0, :])
            print("Finsihed precalculating.")
        
        if maxnum is None:
            maxnum = np.sum(dists[0])

        yaxis = np.linspace(maxrange+ dz/2, -dz/2, len(dists[0]))

        def update(val):
            tindex = int(tslider.val)
            mainplot.cla()
            mainplot.grid()
            mainplot.barh(yaxis, dists[tindex], color = "blue", label = "file1", height = dz)
            mainplot.axhline(-200, linestyle = "dashed", color = "green", label = "z = -200 m")
            mainplot.axhline(-1000, linestyle = "dashed", color = "black", label = "z = -1000 m")
            mainplot.text(maxnum / 2, -100, r"$\sum_i \log m_i = {}$".format(round(mass_zone[tindex][0], 2)), fontsize = 20)
            mainplot.text(maxnum / 2, -600, r"$\sum_i \log m_i = {}$".format(round(mass_zone[tindex][1], 2)), fontsize = 20)
            mainplot.text(maxnum / 2, np.mean([-1000, maxrange]), r"$\sum_i m_i = {}$".format(round(mass_zone[tindex][2], 2)), fontsize = 20)
            mainplot.set_ylim([maxrange, 0])
            # mainplot.set_xlim([0, maxnum])
            mainplot.set_xlabel('log(mass)')
            mainplot.set_ylabel("z")
            mainplot.set_title(title, fontweight = "bold")
            fig.canvas.draw_idle()
            mainplot.legend()
        
        def update2(val):
            tindex = int(tslider.val)
            mainplot.cla()
            mainplot.grid()
            mainplot.barh(yaxis, dists[tindex], color = "blue", label = "file1", height = dz)
            if self.show2:
                mainplot.barh(yaxis, dists2[tindex], color = "red", alpha = 0.7, label = "file2", height = dz)
            mainplot.axhline(-200, linestyle = "dashed", color = "green", label = "z = -200 m")
            mainplot.axhline(-1000, linestyle = "dashed", color = "black", label = "z = -1000 m")
            if self.show2:
                mainplot.text(maxnum / 2, -100, r"file1: $\sum_i \log m_i = {}$; file2: $\sum_i \log m_i = {}$".format(round(mass_zone[tindex][0], 2), round(mass_zone2[tindex][0], 2)), fontsize = 20)
                mainplot.text(maxnum / 2, -600, r"file1: $\sum_i \log m_i = {}$; file2: $\sum_i \log m_i = {}$".format(round(mass_zone[tindex][1], 2), round(mass_zone2[tindex][1], 2)), fontsize = 20)
                mainplot.text(maxnum / 2, np.mean([-1000, maxrange]), r"file1: $\sum_i \log m_i = {}$; file2: $\sum_i \log m_i = {}$".format(round(mass_zone[tindex][2], 2), round(mass_zone2[tindex][2], 2)), fontsize = 20)
            else:
                mainplot.text(maxnum / 2, -100, r"file1: $\sum_i \log m_i = {}$".format(round(mass_zone[tindex][0], 2)), fontsize = 20)
                mainplot.text(maxnum / 2, -600, r"file1: $\sum_i \log m_i = {}$".format(round(mass_zone[tindex][1], 2)), fontsize = 20)
                mainplot.text(maxnum / 2, np.mean([-1000, maxrange]), r"file1: $\sum_i \log m_i = {}$".format(round(mass_zone[tindex][2], 2)), fontsize = 20)    
            mainplot.set_ylim([maxrange, 0])
            # mainplot.set_xlim([0, maxnum])
            mainplot.set_xlabel('log(mass)')
            mainplot.set_ylabel("z")
            mainplot.set_title(title, fontweight = "bold")
            fig.canvas.draw_idle()
            mainplot.legend()
        
        def show_file2(val):
            self.show2 = not self.show2
            update2(0)
        
        update(0)
        
        if file2 is not None:
            tslider.on_changed(update2)
            check.on_clicked(show_file2)
        else:
            tslider.on_changed(update)
        
        plt.show()
        
        return fig, mainplot, sliderax, tslider
    
    def plot_3D(self, k=1, markevery = 1, plot_floor = False):

        """
        A 3d visualization of the simulation.
        
        Parameters
        ----------
        k: float
            The scale factor for markersize (default: 1)
        markevery: int
            Slice number of particles for plotting (default: 1)
        plot_floor: bool
            Plot a wireframe of sea floor (default: False)"""

        plt.close()

        cmap, sm, c = self.get_cmap(self.time, "jet")

        fig = plt.figure(figsize=self.figsize)
        ax = fig.add_subplot(projection='3d')
        
        lon = self.o.get_property('lon')
        lon = np.ma.filled(lon, np.nan)

        lat = self.o.get_property("lat")
        lat = np.ma.filled(lat, np.nan)

        z = self.o.get_property('z')
        z = np.ma.filled(z, np.nan)

        mass = self.o.get_property("mass")
        mass = np.ma.filled(mass, 0)

        for i in range(len(self.time)):
            if i == 0:
                idx = np.where((mass[1, :] < 1) & (z[1, :] < 0))[0]
                ax.scatter(lon[i, idx][::markevery], lat[i, idx][::markevery], z[i, idx][::markevery], s=mass[i, idx][::markevery] * k, c = cmap(i), zorder = 3)
            else:
                idx = np.where((mass[i, :] < 1) & (z[i, :] < 0))[0]
                ax.scatter(lon[i, idx][::markevery], lat[i, idx][::markevery], z[i, idx][::markevery], s=mass[i, idx][::markevery] * k, c = cmap(i), zorder = 2)
        
        if plot_floor:
            cwd = os.getcwd()
            os.chdir("/home/perharic/Documents")
            depth = Dataset("etopo2_cmems_med.nc")
            os.chdir(cwd)
            lons, lats = depth["lon"][:], depth["lat"][:]
            lons, lats = np.meshgrid(lons, lats)
            depth = depth["topo"][:] 
            depth[depth >= 0] = np.nan
            ax.plot_wireframe(lons, lats, depth, alpha = 0.2, color = "black", zorder = 1)
            


        cb1 = plt.colorbar(sm, shrink = 0.5, ax = ax)
        cb1.set_ticks([1, len(self.time)])
        cb1.set_ticklabels([self.time[0], self.time[-1]])
        cb1.set_label("time")

        ax.yaxis.labelpad=20
        ax.xaxis.labelpad=20
        ax.zaxis.labelpad=10
        ax.set_xlabel("lon")
        ax.set_ylabel("lat")
        ax.set_zlabel("z")
        plt.show()
    
    def animate_3D(self, k=1, outfile = None):
        """
        A 3d animation of the simulation.
        
        Parameters
        ----------
        k: float
            The scale factor for markersize (default: 1)
        outfile: str
            File name for saved animation (default: None)"""
        
        plt.close()

        fig = plt.figure(figsize=self.figsize)
        ax = fig.add_subplot(projection='3d')
        
        lon = self.o.get_property('lon')
        lon = np.ma.filled(lon, np.nan)

        lat = self.o.get_property("lat")
        lat = np.ma.filled(lat, np.nan)

        z = self.o.get_property('z')
        z = np.ma.filled(z, np.nan)

        mass = self.o.get_property("mass")
        mass = np.ma.filled(mass, np.nan)

        ax.scatter(lon[0, :], lat[0, :], z[0, :], s=mass[0, :] * k, c = "blue")

        ax.yaxis.labelpad=30
        ax.xaxis.labelpad=30
        ax.zaxis.labelpad=30
        ax.set_xlabel("lon")
        ax.set_ylabel("lat")
        ax.set_zlabel("z")
        xlim = [np.min(lon[np.invert(np.isnan(lon))]), np.max(lon[np.invert(np.isnan(lon))])]
        ylim = [np.min(lat[np.invert(np.isnan(lat))]), np.max(lat[np.invert(np.isnan(lat))])]
        zlim = [np.min(z[np.invert(np.isnan(z))]), np.max(z[np.invert(np.isnan(z))])]
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_zlim(*zlim)
        ax.set_title(self.time[0])
        ax.view_init(elev=10, azim=45, roll=0)

        def update(i):
            ax.cla()
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
            ax.set_zlim(*zlim)
            ax.yaxis.labelpad=30
            ax.xaxis.labelpad=30
            ax.zaxis.labelpad=30
            ax.set_xlabel("lon")
            ax.set_ylabel("lat")
            ax.set_zlabel("z")
            idx = np.invert((z[1, :] >= 0) & (mass[1, :] == 1))
            ax.scatter(lon[i, idx], lat[i, idx], z[i, idx], s=mass[i, idx] * k, c = "blue")
            ax.set_title(self.time[i])
            ax.view_init(elev=10, azim=45, roll=0)

        frames = len(self.time)

        anim=animation.FuncAnimation(fig, update, blit=False, frames = frames, interval=100)

        if outfile is not None:
            writer = animation.PillowWriter(fps=10)
            anim.save(outfile, writer = writer)

        plt.show()
    
    def number_of_decays(self):

        """Plot the total number of decays and the decay number per particle at each time step."""

        mass = self.o.get_property("mass")
        mass = np.ma.filled(mass, np.nan)

        @jit(nopython = True)
        def get_decays(mass, time):
            decay_count = []
            specific_dc = []
            last = len(np.where(mass[0, :] > 0)[0])
            for i in range(1, len(time)):
                new = len(np.where(mass[i, :] > 0)[0])

                decay_count.append(new - last)
                specific_dc.append((new - last) / last)

                last = new
            return decay_count, specific_dc
        
        decay_count, specific_dc = get_decays(mass, self.time)
        self.ax.plot(self.time[1:], decay_count, color = "blue", label = "total decay number")
        ax2 = plt.twinx(self.ax)
        ax2.plot(self.time[1:], specific_dc, color = "black", linestyle = "dashed", label = "specific decay number")

        self.ax.set_ylabel("Total number of decays")
        ax2.set_ylabel("Number of decays per particle")
        self.ax.set_xlabel("time")
        self.ax.tick_params(axis='y', colors='blue')
        self.fig.legend()

        plt.show()

    def zone_mass_distribution(self, lons, lats, file2 = None, m0 = None, pathtofile2 = None, read_file = None, write_file = None):
        plt.close()

        if read_file is None:
            import os
            curdir = os.getcwd()
            os.chdir(pathtofile2)
            o2 = Plot(file2, self.distribution)
            os.chdir(curdir)
            plt.close()
        
        heights = []
        heights2 = []

        for i, zone in enumerate([-200, -1000, "sea floor"]):
            if read_file is not None:
                if m0 is None:
                    raise ValueError("Initial mass must be provided, if mass is to be read from saved file.")
                mass = np.loadtxt(str(zone).replace(" ", "") + "_2000_" + read_file)
                mass2 = np.loadtxt(str(zone).replace(" ", "") + "_2003_" + read_file)
            else:
                bad1 = self.clean_dataset()
                bad2 = o2.clean_dataset()
                mass, m0 = self.zone_crossing_event(lons, lats, zone, bad1, bad2)
                mass2, m02 = o2.zone_crossing_event(lons, lats, zone, bad1, bad2)
                if m0 != m02:
                    raise ValueError("Initial mass of files do not match.")
            total_mass = np.sum(mass)
            total_mass2 = np.sum(mass2)
            
            if write_file is not None:
                np.savetxt(str(zone).replace(" ", "") + "_2000_" + write_file, mass)
                np.savetxt(str(zone).replace(" ", "") + "_2003_" + write_file, mass2)
            
            heights.append(total_mass / m0)
            heights2.append(total_mass2 / m0)
        
        x = np.arange(3)

        plt.bar(x-0.25, heights, width = 0.4, color = "navy", label = "2000")
        plt.bar(x+0.25, heights2, width = 0.4, color = "red", label = "2003")
        plt.legend()
        plt.xticks(x, ["-200 m", "-1000 m", "sea floor"])
        plt.xlabel("depth")
        plt.ylabel(r"$\sum m(z = \mathrm{depth}) / m_0$")
        plt.show()
    
    def mass_map(self, lons, lats, file2, pathtofile2 = None, read_file = None, write_file = None):
        """Plot the mass reached at depths [-200m, -1000m, sea_floor] over a cartopy map.

        Parameters
        -----------
        lons: 1D ndarray
            Longtidue values for grid.
        lats: 1D ndarray
            Latitude values for grid.
        file2: Name of file for comparison.
        pathotofile2: str
            Path to file2 (default None is same as original file)
        read_file: str
            Name of file with saved values, to avoid recalculating. Should be the same as the name used with write_file.
        write_file: str
            Name of file to store values.
        """

        plt.close()

        if read_file is None:
            if pathtofile2 is not None:
                import os
                curdir = os.getcwd()
                os.chdir(pathtofile2)
                o2 = Plot(file2, self.distribution)
                os.chdir(curdir)
            else:
                o2 = Plot(file2, self.distribution)
            plt.close()
        
        fig, axes = plt.subplots(3, 3, figsize=self.figsize, subplot_kw={'projection': ccrs.PlateCarree()})
        
        cmaps = ["Reds", "Blues", "Greens"]
        cmap_diffs = mpl.cm.RdBu_r
        subfig_titles = ["climatology", "2003", "difference"]

        #linear
        # vmin = [0.73, 0.25, 0.0]
        # vmax = [0.85, 0.38, 0.1]
        
        #linear_dm/m
        # vmin_diff = [-0.016, -0.028, -0.04]
        # vmax_diff = [0.02, 0.05, 0.07]
        # cbar_ticks = [[-0.01, 0, 0.01, 0.02], [-0.02, 0.00, 0.02, 0.04], [-0.03, 0, 0.03, 0.06]]

        #linear dm/m0
        # vmin_diff = [-0.016, -0.01, -0.01]
        # vmax_diff = [0.025, 0.02, 0.02]
        # cbar_ticks = [[-0.01, 0, 0.01, 0.02], [-0.01, 0.00, 0.01], [-0.01, 0, 0.01, 0.02]]

        #exp
        vmin = [0.6, 0.18, 0.0]
        vmax = [0.8, 0.3, 0.03]
        
        #exp_dm/m
        # vmin_diff = [-0.08, -0.08, -0.08]
        # vmax_diff = [0.15, 0.15, 0.15]
        # cbar_ticks = [[-0.05, 0, 0.05, 0.10], [-0.05, 0.0, 0.05, 0.10], [0, 0.1]]

        #exp dm/m0
        vmin_diff = [-0.05, -0.015, -0.02]
        vmax_diff = [0.1, 0.03, 0.04]
        cbar_ticks = [[-0.05, 0, 0.05, 0.10], [0.0, 0.01, 0.02], [-0.02, 0, 0.02, 0.04]]

        plt.subplots_adjust(wspace=0, hspace=0.2)

        #Manualy create colorbar axes.
        def cb_ax(ax1, ax2, i):
            left = (ax1.get_position().x0 + ax1.get_position().x1) / 2 - 0.2
            width = (ax2.get_position().x0 + ax2.get_position().x1) / 2 - left +0.1
            bottom = ax1.get_position().y0 - 0.005 - 0.03*i
            return fig.add_axes([left, bottom, width, 0.01])
        
        #Manualy create colorbar axes.
        def cb_ax2(ax1, i):
            left = ax1.get_position().x0 + 0.05
            width = ax1.get_position().x1 - ax1.get_position().x0
            bottom = ax1.get_position().y0 - 0.005 - 0.03*i
            return fig.add_axes([left, bottom, width, 0.01])
        
        for i, zone in enumerate([-200, -1000, "sea floor"]):
            ax = axes[i, :]
            if read_file is not None:
                mass1 = np.loadtxt(str(zone).replace(" ", "") + "_clim_" + read_file)
                mass2 = np.loadtxt(str(zone).replace(" ", "") + "_2003_" + read_file)
                mass = np.loadtxt(str(zone).replace(" ", "") + "_diff_" + read_file)
                # mass[mass > 0.8] = np.nan
                mass = mass1 - mass2
            else:
                bad1 = self.clean_dataset()
                bad2 = o2.clean_dataset()
                mass1, m0 = self.zone_crossing_event(lons, lats, zone, bad1, bad2)
                mass2, m02 = o2.zone_crossing_event(lons, lats, zone, bad1, bad2)
                mass = (np.copy(mass1) - mass2) / np.copy(mass1)
            
            if write_file is not None:
                np.savetxt(str(zone).replace(" ", "") + "_clim_" + write_file, mass1)
                np.savetxt(str(zone).replace(" ", "") + "_2003_" + write_file, mass2)
                np.savetxt(str(zone).replace(" ", "") + "_diff_" + write_file, mass)
            
            #Clip values of third column for better visualization.
            m1, M1 = np.min(mass1[np.invert(np.isnan(mass1))]), np.max(mass1[np.invert(np.isnan(mass1))])
            m2, M2 = np.min(mass2[np.invert(np.isnan(mass2))]), np.max(mass2[np.invert(np.isnan(mass2))])
            mass = np.copy(mass)
            row, col = np.where(mass < vmin_diff[i])
            mass[row, col] = vmin_diff[i]
            row, col = np.where(mass > vmax_diff[i])
            mass[row, col] = vmax_diff[i]
            m, mid, M = self.get_colormap_midpoint(mass)
            Vmin, Vmax = min(m1, m2), max(M1, M2)

            #Clip values of first two columns for better visualization. 
            mass1 = np.copy(mass1)
            mass2 = np.copy(mass2)
            row, col = np.where(mass1 < vmin[i])
            mass1[row, col] = vmin[i]
            row, col = np.where(mass1 > vmax[i])
            mass1[row, col] = vmax[i]
            row, col = np.where(mass2 < vmin[i])
            mass2[row, col] = vmin[i]
            row, col = np.where(mass2 > vmax[i])
            mass2[row, col] = vmax[i]
            
            #Add coastlines.
            for j in range(3):
                # ax[j].add_feature(cartopy.feature.LAND, facecolor="gray",edgecolor='black', zorder = 2, scale='10m')
                ax[j].coastlines(zorder = 3, resolution='10m')
                if i == 0:
                    ax[j].set_title(subfig_titles[j], fontweight= "bold", fontsize = 15)
            
            #Plot.
            ax[0].contourf(lons, lats, mass1.T, 20, transform=ccrs.PlateCarree(), cmap = cmaps[i], zorder = 1,extend='both', extendfrac='auto')
            sm2 = ax[1].contourf(lons, lats, mass2.T, 20, transform=ccrs.PlateCarree(), cmap = cmaps[i], zorder = 1,extend='both', extendfrac='auto')
            cmap = self.shiftedColorMap(cmap_diffs, midpoint = mid, name='shifted{}'.format(i))
            sm = ax[2].contourf(lons, lats, mass.T, 20, transform=ccrs.PlateCarree(), cmap = cmap, zorder = 1,extend='both', extendfrac='auto')

            #Create colorbars for first two columns.
            cb2axes = cb_ax(ax[0], ax[1], i)
            cb2 = plt.colorbar(sm2, cax = cb2axes, orientation="horizontal", shrink = 1)
            cb2.set_label(r"$m/m_0$")
            cb2.ax.set_xticks(cb2.get_ticks()[::2])
            cb2.ax.set_xticklabels(["{:.2f}".format(j) for j in cb2.get_ticks()])
            
            #Create colorbars for third column.
            cb = plt.colorbar(sm, cax = cb_ax2(ax[2], i), orientation="horizontal")
            cb.set_label(r"$\Delta m/m_0$")
            cb.ax.set_xticks(cbar_ticks[i])
            if i == 0:
                cb.ax.set_xticklabels(["{:.2f}".format(round(j, 2)) for j in cb.get_ticks()])
        
        plt.savefig("tmp.pdf", dpi = 300, bbox_inches = "tight")
        # plt.show()
    
    def zone_crossing_event(self, lons, lats, depth, bad1, bad2):
        """Calculate particle properties when crossing a certain depth."""

        z = self.o.get_property('z')
        z = np.ma.filled(z, np.nan)

        mass = self.o.get_property("mass")
        mass = np.ma.filled(mass, np.nan)

        lon = self.o.get_property("lon")
        lon = np.ma.filled(lon, np.nan)

        lat = self.o.get_property("lat")
        lat = np.ma.filled(lat, np.nan)
        status = self.o.get_property("status")
        status = np.ma.MaskedArray(status.data, status.mask, float)
        status = np.ma.filled(status, np.nan)

        # @jit(nopython = True)
        def mass_sum(lons, lats, depth, z, mass, lon, lat, bad1 = bad1, bad2 = bad2):
            mass_at_depth = np.zeros((len(lons), len(lats)))
            m0 = 0
            for i in range(mass.shape[1]):

                #Remove bad particles.
                if i in bad1 or i in bad2:
                    continue

                drifter_trajectory = z[:, i]
                drifter_mass = mass[:, i]

                trajectory_nans = np.invert(np.isnan(drifter_trajectory))

                m0 += drifter_mass[0]
                
                if type(depth) != str:
                    if np.min(drifter_trajectory[trajectory_nans])>depth:
                        continue
                    depth_idx = np.argmin(np.abs(drifter_trajectory[trajectory_nans] - depth))
                
                else:
                    #Sea_floor.
                    bounce = np.where(drifter_trajectory[:-1] < np.roll(drifter_trajectory, -1)[:-1])[0]
                    if len(bounce) > 0:
                        raise ValueError("Depth decreases.")
                    else:
                        s = status[:, i]
                        sea_floor = np.where(s == 2)[0]
                        if len(sea_floor) >0:
                            depth_idx = sea_floor[0]
                        else:
                            depth_idx = len(drifter_trajectory)-1


                if np.isnan(drifter_mass[depth_idx]) or drifter_mass[depth_idx] == 1:
                    continue

                drifter_lon = lon[depth_idx, i]
                drifter_lat = lat[depth_idx, i]

                if np.isnan(drifter_lon) or np.isnan(drifter_lat):
                    continue
                
                #Find closest position on grid.
                lon_grid = np.argmin(np.abs(lons - drifter_lon))
                lat_grid = np.argmin(np.abs(lats - drifter_lat))

                mass_at_depth[lon_grid, lat_grid] = mass_at_depth[lon_grid, lat_grid] + drifter_mass[depth_idx]
            
            return mass_at_depth, m0
        
        return mass_sum(lons, lats, depth, z, mass, lon, lat)

    def mass_z(self, file2):
        """Plot z(mass) for two files at specified locations."""

        plt.close()
        fig, ax = plt.subplots(2, 2, sharey = True, sharex = "row", figsize = self.figsize)

        o2 = Plot(file2, self.distribution)
        z2 = o2.o.get_property('z')
        z2 = np.ma.filled(z2, np.nan)

        mass2 = o2.o.get_property("mass")
        mass2 = np.ma.filled(mass2, np.nan)
        plt.close()
        
        z = self.o.get_property('z')
        z = np.ma.filled(z, np.nan)

        mass = self.o.get_property("mass")
        mass = np.ma.filled(mass, np.nan)

        lons = self.o.get_property('lon')
        lons = np.ma.filled(lons, np.nan)

        lats = self.o.get_property("lat")
        lats = np.ma.filled(lats, np.nan)

        trajectories = self.find_closest_location(lons, lats)
        bad1 = self.clean_dataset()
        bad2 = o2.clean_dataset()
        for k in range(len(ax)):
            trajectory = trajectories[k]
            if (trajectory in bad1) or (trajectory in bad2):
                raise ValueError("This trajectory has a problem. Change the location coordinates slightly to select a new trajectory.")
            
            ax[0, k].plot(mass[:, trajectory], z[:, trajectory], label = "climatology", color = "purple", linewidth = 2)
            ax[0, k].plot(mass2[:, trajectory], z2[:, trajectory], label = "2003", color = "orange", linewidth = 2)
            #dm/m
            ax[1, k].plot(np.divide(mass[:, trajectory] - mass2[:len(mass[:, trajectory]), trajectory], mass[:, trajectory]), z[:, trajectory], color = "blue", linewidth = 2)
            #dm/m_0
            # ax[1, k].plot(mass[:, trajectory] - mass2[:len(mass[:, trajectory]), trajectory], z[:, trajectory], color = "blue", linewidth = 2)
            ax[0, k].set_title("L{}".format(k + 1), fontweight="bold")
            ax[0, k].set_xlabel(r"$m(z)/m_0$")
            ax[1, k].set_xlabel(r"$\Delta m(z)/m$")
            ax[0, k].grid()
            ax[1, k].grid()
        
        ax[0, 0].legend()
        ax[0, 0].set_ylabel("depth [m]")
        ax[1, 0].set_ylabel("depth [m]")
        plt.show()

    def time_z(self, file2):
        """Plot z(t) for two files at specified locations."""

        plt.close()
        fig, ax = plt.subplots(1, 2, sharey = True, figsize = self.figsize)

        o2 = Plot(file2, self.distribution)
        z2 = o2.o.get_property('z')
        z2 = np.ma.filled(z2, np.nan)
        plt.close()
        
        z = self.o.get_property('z')
        z = np.ma.filled(z, np.nan)

        lons = self.o.get_property('lon')
        lons = np.ma.filled(lons, np.nan)

        lats = self.o.get_property("lat")
        lats = np.ma.filled(lats, np.nan)

        trajectories = self.find_closest_location(lons, lats)

        for k in range(len(ax)):
            trajectory = trajectories[k]
            zplot = z[:, trajectory]
            zplot = zplot[np.invert(np.isnan(zplot))]
            z2plot = z2[:, trajectory]
            z2plot = z2plot[np.invert(np.isnan(z2plot))]
            ax[k].plot(self.create_timedelta_array(len(zplot)), zplot, label = "climatology", color = "purple", linewidth = 2)
            ax[k].plot(self.create_timedelta_array(len(z2plot)), z2plot, label = "2003", color = "orange", linewidth = 2)
            ax[k].set_title("L{}".format(k + 1), fontweight="bold")
            ax[k].set_xlabel("sinking time [h]")
            ax[k].grid()
        
        ax[0].legend()
        ax[0].set_ylabel("depth [m]")
        plt.show()
    
    def hovmoller(self, dz = 10, maxdepth = 200, title = None, file2 = None, title2 = None):

        """Plot a hovmoller type of plot of mass distribution allong the z axis for all times."""

        plt.close()
        fig, ax = plt.subplots(1, 1 if file2 is None else 2, figsize = self.figsize)

        if title is None:
            title = "Vertical mass distribution"
            title2 = title
        
        z = self.o.get_property('z')
        z = np.ma.filled(z, np.nan)
        mass = self.o.get_property("mass")
        mass = np.ma.filled(mass, 0)

        lon = self.o.get_property("lon")
        lon = np.ma.filled(lon, np.nan)

        lat = self.o.get_property("lat")
        lat = np.ma.filled(lat, np.nan)

        if file2 is not None:
            o2 = Plot(file2, distribution=self.distribution)
            z2 = o2.o.get_property('z')
            z2 = np.ma.filled(z2, np.nan)
            mass2 = o2.o.get_property("mass")
            mass2 = np.ma.filled(mass2, 0)

            lon2 = o2.o.get_property("lon")
            lon2 = np.ma.filled(lon2, np.nan)

            lat2 = o2.o.get_property("lat")
            lat2 = np.ma.filled(lat2, np.nan)
            plt.close()

        if maxdepth is None:
            maxrange = np.min(z[np.invert(np.isnan(z))])
            if file2 is not None:
                maxrange2 = np.min(z2[np.invert(np.isnan(z2))])
        else:
            maxrange = -np.abs(maxdepth)
            if file2 is not None:
                maxrange2 = maxrange
        

        @jit(nopython = True)
        def mass_dist(mass, z, m0, maxrange, dz=dz):
            dist = []
            zmax = 0
            zmin = -dz
            while zmin > maxrange:
                dist.append(np.sum(mass[np.where((z > zmin) & (z<=zmax) & np.invert(np.isnan(z))& np.invert(np.isnan(mass)))[0]]) / m0)
                zmax = zmin
                zmin -= dz
            dist.append(np.sum(mass[np.where((z >= zmin) & (z<zmax)& np.invert(np.isnan(z))& np.invert(np.isnan(mass)))[0]]) / m0)
            return dist
        
        @jit(nopython=True)
        def precalculate(mass, z, lon, maxrange = maxrange):
            dists = []
            at_sea = np.invert(np.isnan(lon))
            for timestep in range(mass.shape[0]):
                mi = mass[timestep, at_sea]
                zi = z[timestep, at_sea]
                if timestep > 0:
                    idx = np.where((mi < 1) &np.invert(np.isnan(mi)) &(zi <0))[0]
                    mi = np.copy(mi)[idx]
                    zi = np.copy(zi)[idx]
                else:
                    m0 = np.sum(mi)
                dists.append(mass_dist(mi, zi, m0, maxrange))
                # if dists[-1][-1] > 0:
                #     break
            return dists
        
        dists = np.asarray(precalculate(mass, z, lon[0, :]))
        dists[np.where(dists == 0)[0], np.where(dists == 0)[1]] = np.nan

        if file2 is not None:
            dists2 = np.asarray(precalculate(mass2, z2, lon2[0, :], maxrange2))
            dists2[np.where(dists2 == 0)[0], np.where(dists2 == 0)[1]] = np.nan

        if file2 is None:
            ax = [ax, 0]

        xslice = 70
        yslice = 14

        im = ax[0].imshow(dists.T, cmap = "Reds", norm = mpl.colors.PowerNorm(gamma = 0.17))
        cb = plt.colorbar(im, ax = ax[0], orientation = "horizontal")
        cb.set_label(r"$\sum_i m_i/m_i(t=0)$")
        cb.set_ticks([0, 0.2, 0.4, 0.6, 1])

        ax[0].set_xlabel("time")
        ax[0].set_ylabel(r"$z$ [m]")
        ax[0].set_xticks(np.arange(0, len(dists), xslice))
        ax[0].set_xticklabels(self.time[:len(dists):xslice])
        l = ax[0].get_xticklabels()
        for i in l:
            i._text = i._text[8:]
        ax[0].set_xticklabels(l)

        yticks = np.arange(0, len(dists[0]), yslice)
        ax[0].set_yticks(yticks)
        ax[0].set_yticklabels(-yticks *dz - dz/2)
        ax[0].set_title(title, fontweight="bold")
        plt.gca().grid()
        if file2 is not None:
            im = ax[1].imshow(dists2.T, cmap = "Reds", norm = mpl.colors.PowerNorm(gamma = 0.17))
            cb = plt.colorbar(im, ax = ax[1], orientation = "horizontal")
            cb.set_label(r"$\sum_i m_i/m_i(t=0)$")
            cb.set_ticks([0, 0.2, 0.4, 0.6, 1])

            ax[1].set_xlabel("time")
            ax[1].set_ylabel(r"$z$ [m]")
            ax[1].set_xticks(np.arange(0, len(dists), xslice))
            ax[1].set_xticklabels(self.time[:len(dists):xslice])
            ax[1].set_xticklabels(self.time[:len(dists):xslice])
            l = ax[1].get_xticklabels()
            for i in l:
                i._text = i._text[8:]
            ax[1].set_xticklabels(l)

            yticks = np.arange(0, len(dists[0]), yslice)
            ax[1].set_yticks(yticks)
            ax[1].set_yticklabels(-yticks *dz - dz/2)
            ax[1].set_title(title2, fontweight="bold")
        ax[0].grid()
        plt.tight_layout()
        plt.show()
    
    def area(self, lons, lats):
        """Plot the position of specified locations on a cartopy map.

        Parameters
        ----------
        lons: 1D ndarray
            Longtidue values for grid.
        lats: 1D ndarray
            Latitude values for grid.
        """

        from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

        plt.close()
        fig, ax = plt.subplots(1, 1, subplot_kw={'projection': ccrs.PlateCarree()}, figsize = self.figsize)
        
        lons, lats = np.meshgrid(lons, lats)
        lon = np.min(lons)
        lat = np.min(lats)
        dx = np.max(lons) - np.min(lons)
        dy = np.max(lats) - np.min(lats)

        lons = self.o.get_property('lon')
        lons = np.ma.filled(lons, np.nan)

        lats = self.o.get_property("lat")
        lats = np.ma.filled(lats, np.nan)

        trajectories = self.find_closest_location(lons, lats)

        for i, t in enumerate(trajectories):
            x, y = lons[0, t], lats[0, t]
            plt.plot(x, y,"o", color = "navy", markersize = 10)
            plt.text(x+0.01, y-0.1, "L{}".format(i + 1), color = "navy", horizontalalignment="left", verticalalignment="top",fontsize = 35, fontweight = "bold")
        
        simulation_area = mpl.patches.Rectangle((lon, lat), dx, dy, zorder = 3, facecolor="none", edgecolor="red", linewidth=6)
        ax.set_extent([-5, 20, 33, 48], crs=ccrs.PlateCarree())
        ax.add_feature(cartopy.feature.LAND, facecolor="gray",edgecolor='black', zorder = 1)
        ax.coastlines(zorder = 2)
        gl = ax.gridlines(draw_labels=False, linewidth=0.5, alpha=0.4, color='k', linestyle='--')
        # gl.xlocator = mpl.ticker.FixedLocator(np.round([lon, *[lons[0, i] for i in trajectories], lon + dx], 1))
        # gl.ylocator = mpl.ticker.FixedLocator(np.round([lat, *[lats[0, i] for i in trajectories], lat + dy], 1))
        # gl.xformatter = LONGITUDE_FORMATTER
        # gl.yformatter = LATITUDE_FORMATTER
        ax.add_patch(simulation_area)

        plt.show()

    def find_closest_location(self, lons, lats):
        """Find the drifters with initial position closest to the locations specified."""

        l = lons[0, np.invert(np.isnan(lons[0, :]))]
        L = lats[0, np.invert(np.isnan(lats[0, :]))]
        trajectories = []
        for i, j in self.loc:
            lon = np.where(np.abs(lons[0, :] - i) == np.min(np.abs(l - i)))[0]
            lat = np.where(np.abs(lats[0, lon] - j) == np.min(np.abs(L - j)))
            trajectories.append(lon[lat[0][0]])
        return trajectories
    
    def create_timedelta_array(self, n):
        dt = datetime.strptime(str(self.time[1]),'%Y-%m-%d %H:%M:%S') - datetime.strptime(str(self.time[0]),'%Y-%m-%d %H:%M:%S')
        dt = dt.total_seconds() / 3600
        return np.arange(0, n * dt, dt)
    
    def m_T_correlations(self, markevery = 1):
        from matplotlib.patches import Ellipse
        import matplotlib.transforms as transforms

        z = self.o.get_property('z')
        z = np.ma.filled(z, np.nan)

        mass = self.o.get_property("mass")
        mass = np.ma.filled(mass, np.nan)
        T = self.o.get_property("sea_water_temperature")
        T = np.ma.filled(T, np.nan)
        lon = self.o.get_property("lon")
        lon = np.ma.filled(lon, np.nan)
        lat = self.o.get_property("lat")
        lat = np.ma.filled(lat, np.nan)
        bad_trajectories = self.clean_dataset()

        def mass_sum(depth, z, mass, lon, lat, T, bad = bad_trajectories):

            def grad_calc(T, z, i, j = 0):
                return (T[j] - T[i]) / (z[j] - z[i])
            
            mass_at_depth = []
            T_surface = []
            T_100 = []
            T_200 = []
            T_grad_100_200 = []
            T_grad_100 = []
            T_grad_200 = []

            for i in range(mass.shape[1]):
                if i in bad:
                    bad.remove(i)
                    continue
                drifter_trajectory = z[:, i]
                drifter_mass = mass[:, i]

                trajectory_nans = np.invert(np.isnan(drifter_trajectory))

                if np.min(drifter_trajectory[trajectory_nans])>depth:
                    continue
                
                depth_idx = np.argmin(np.abs(drifter_trajectory[trajectory_nans] - depth))

                idx_100 = np.argmin(np.abs(drifter_trajectory[trajectory_nans] + 15))

                if np.isnan(drifter_mass[depth_idx]) or drifter_mass[depth_idx] == 1:
                    continue

                drifter_lon = lon[depth_idx, i]
                drifter_lat = lat[depth_idx, i]

                if np.isnan(drifter_lon) or np.isnan(drifter_lat):
                    continue

                if np.isnan(drifter_mass[depth_idx]) or drifter_mass[depth_idx] == 1:
                    continue

                mass_at_depth.append(drifter_mass[depth_idx])

                T_surface.append(T[0, i])
                T_100.append(T[idx_100, i])
                T_200.append(T[depth_idx, i])
                
                T_grad_100_200.append(grad_calc(T[:, i], drifter_trajectory, depth_idx, idx_100))
                T_grad_100.append(grad_calc(T[:, i], drifter_trajectory, idx_100))
                T_grad_200.append(grad_calc(T[:, i], drifter_trajectory, depth_idx))

            return mass_at_depth, [T_surface, T_100, T_200, T_grad_100, T_grad_200, T_grad_100_200]
    
        mass_at_depth, temperature = mass_sum(-200, z, mass, lon, lat, T)

        def correlation_elipse(x, y):
            cov = np.cov(x, y)
            cov_coef = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
            print(cov_coef)
            radius_x = np.sqrt(1 + cov_coef)
            radius_y = np.sqrt(1 - cov_coef)
            return np.mean(x), np.mean(y), np.sqrt(cov[0, 0]), np.sqrt(cov[1, 1]), radius_x, radius_y
        

        plt.close()
        fig, ax = plt.subplots(2, 3, sharey=True, figsize = self.figsize)
        ax = ax.flatten()
        xlabels = ["T(z = 0 m)", "T(z = -15 m)", "T(z = -200 m)", r"$\frac{T_0-T_{15}}{0 + 15}$", r"$\frac{T_0-T_{200}}{0 + 200}$", r"$\frac{T_{15}-T_{200}}{-15 + 200}$"]
        for i in range(len(temperature)):
            cov = correlation_elipse(temperature[i], mass_at_depth)
            ellipse = Ellipse((0, 0), width= cov[4] * 2, height = cov[5] * 2, facecolor = "none", edgecolor = "red", zorder = 2, linewidth = 3)
            transf = transforms.Affine2D() \
                    .rotate_deg(45) \
                    .scale(cov[2], cov[3]) \
                    .translate(cov[0], cov[1])

            ellipse.set_transform(transf + ax[i].transData)
            ax[i].add_patch(ellipse)
            ax[i].plot(temperature[i][::markevery], mass_at_depth[::markevery], ".", zorder = 1, color = "black")
            ax[i].set_xlabel(xlabels[i])
            print(max(temperature[i]))
        ax[0].set_ylabel(r"$m_i/m_0$")
        ax[3].set_ylabel(r"$m_i/m_0$")
        plt.show()

    def clean_dataset(self):
        """Get indicies of trajectories which are on land / do not interact (mass stays one forever) or had a problem with reading the temperature."""

        from global_land_mask import globe

        z = self.o.get_property('z')
        z = np.ma.filled(z, np.nan)

        mass = self.o.get_property("mass")
        mass = np.ma.filled(mass, np.nan)

        T = self.o.get_property("sea_water_temperature")
        T = np.ma.filled(T, np.nan)

        lon = self.o.get_property("lon")
        lon = np.ma.filled(lon, np.nan)

        lat = self.o.get_property("lat")
        lat = np.ma.filled(lat, np.nan)

        bad_trajectories = []

        for i in range(mass.shape[1]):
            drifter_trajectory = z[:, i]
            drifter_mass = mass[:, i]

            if np.isnan(lat[0, i]) or globe.is_land(lat[0, i], lon[0, i]):
                bad_trajectories.append(i)
                continue

            trajectory_nans = np.invert(np.isnan(drifter_trajectory))
            
            if (drifter_trajectory[trajectory_nans]>0).any() or np.any(drifter_mass[1:] == 1):
                bad_trajectories.append(i)
                continue

            if np.min(drifter_trajectory[trajectory_nans])>-100:
                continue

            idx_100 = np.argmin(np.abs(drifter_trajectory[trajectory_nans] + 100))
            
            if T[idx_100, i] > 16:
                bad_trajectories.append(i)
        
        return bad_trajectories
    
    def exp_lin_diff(self, file2, file3, file4):
        def diff(object1, object2, self):
            z2 = object2.get_property('z')
            z2 = np.ma.filled(z2, np.nan)

            mass2 = object2.get_property("mass")
            mass2 = np.ma.filled(mass2, np.nan)
            
            z = object1.get_property('z')
            z = np.ma.filled(z, np.nan)

            mass = object1.get_property("mass")
            mass = np.ma.filled(mass, np.nan)

            lons = object1.get_property('lon')
            lons = np.ma.filled(lons, np.nan)

            lats = object1.get_property("lat")
            lats = np.ma.filled(lats, np.nan)

            trajectory = self.find_closest_location(lons, lats)[0]
            T = object1.get_property("sea_water_temperature")
            T = np.ma.filled(T, np.nan)
            T2 = object2.get_property("sea_water_temperature")
            T2 = np.ma.filled(T2, np.nan)
            return np.divide(mass[:, trajectory] - mass2[:len(mass[:, trajectory]), trajectory], mass[:, trajectory]), z[:, trajectory], T[:, trajectory], T2[:, trajectory]
        
        
        o2 = Plot(file2, self.distribution)
        plt.close()
        # o3 = Plot(file3, self.distribution, locations=self.loc)
        # plt.close()
        # o4 = Plot(file4, self.distribution)
        # plt.close()

        lin, z, T1, T12 = diff(self.o, o2.o, self)
        # exp, z2, T2, T22 = diff(o3.o, o4.o, o3)
        DT1 = T12[:len(T1)] - T1
        # DT2 = T22[:len(T2)] - T2
        # if len(z)> len(z2):
        #     Z = z2
        # else:
        #     Z = z
        # plt.plot(np.divide(exp[:len(Z)], lin[:len(Z)]), Z)
        a, c, d = 0.064, 0.14, 0.145
        # y = c * d / a * np.exp(d * T2[:len(Z)]) * DT2[:len(Z)] / DT1[:len(Z)]
        atT = a * np.arange(0, len(DT1) * 0.5, 0.5) / 24 * T1
        atDT = a * np.arange(0, len(DT1) * 0.5, 0.5) / 24 * DT1
        lin_analytic = (atDT-atDT**2/2) * (np.cosh(atT) - np.sinh(atT)) * np.exp(atT)
        lin_analytic2 = (atDT) * (np.cosh(atT) - np.sinh(atT)) * np.exp(atT)

        # plt.plot(y[:len(Z)], Z)
        plt.plot(lin_analytic, z)
        plt.plot(lin_analytic2, z)
        plt.plot(lin, z)
        plt.show()
    
    def mass_at_depth_time(self, file2):

        plt.close()
        o2 = Plot(file2, self.distribution)
        plt.close()

        z2 = o2.o.get_property('z')
        z2 = np.ma.filled(z2, np.nan)

        mass2 = o2.o.get_property("mass")
        mass2 = np.ma.filled(mass2, np.nan)
        
        z = self.o.get_property('z')
        z = np.ma.filled(z, np.nan)

        mass = self.o.get_property("mass")
        mass = np.ma.filled(mass, np.nan)
        bad1 = self.clean_dataset()
        bad2 = o2.clean_dataset()

        BAD = np.zeros(z.shape[1], dtype = bool)
        BAD[bad1] = True

        crossed_particles = np.zeros(len(BAD), dtype = bool)
        m1 = [0]

        for i in range(0, mass.shape[0]):
            new_crossed_particles = np.logical_and(np.logical_and(np.logical_and(np.logical_and(z[i, :] <= -1000, np.invert(np.isnan(z[i, :]))), np.invert(BAD)), np.invert(crossed_particles)), np.invert(np.isnan(mass[i, :])))
            m1.append(m1[-1] + np.sum(mass[i, new_crossed_particles]))
            crossed_particles[new_crossed_particles] = True
        
        
        m2 = [0]
        BAD = np.zeros(z2.shape[1], dtype = bool)
        BAD[bad2] = True
        crossed_particles = np.zeros(len(BAD), dtype = bool)

        for i in range(0, mass2.shape[0]):
            new_crossed_particles = np.logical_and(np.logical_and(np.logical_and(np.logical_and(z2[i, :] <= -1000, np.invert(np.isnan(z2[i, :]))), np.invert(BAD)), np.invert(crossed_particles)), np.invert(np.isnan(mass2[i, :])))
            m2.append(m2[-1] + np.sum(mass2[i, new_crossed_particles]))
            crossed_particles[new_crossed_particles] = True
        
        plt.plot(self.create_timedelta_array(len(m1)), m1)
        plt.plot(self.create_timedelta_array(len(m2)), m2)
        plt.show()

        

