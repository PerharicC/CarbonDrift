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
    """The alternative for opening .nc files, as I have found that opendrift has problems with understanding masked datasets."""

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

    def __init__(self, filename, distribution, fig_size = (20, 20), fontsize = 20):

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
        # self.o = cdo.open(filename, distribution = distribution) #Doesnt work with masked arrays.
        self.o = Open(filename)
        self.time = self.o.get_time_array()
        self.distribution = distribution

        plt.rcParams.update({'font.size': fontsize})

        fig, ax = plt.subplots(1, 1, figsize = fig_size)
        self.fig = fig
        self.ax = ax
        self.figsize = fig_size

        self.steps_output = self.o.steps_output
    
    @staticmethod
    def get_cmap(x, map):

        """Create a scalar mappable from array."""

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
        
        from mpl_toolkits.axes_grid1 import make_axes_locatable

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

        sm = ax.scatter(lon[0, :], lat[0, :], z[0, :], s=mass[0, :] * k, c = mass[0, :], cmap = "plasma_r", norm = mpl.colors.PowerNorm(vmin = 0, vmax = 1, gamma = 0.2))
        cb = fig.colorbar(sm, ax = ax)
        cb.set_label(r"$m / m_0$")
        cb_ticks = [0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 1]
        cb.set_ticks(cb_ticks)

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
        ax.set_xticks(ax.get_xticks()[::2])
        ax.set_yticks(ax.get_yticks()[::2])

        norm = plt.Normalize(vmin=0, vmax=1)

        def update(i):
            ax.cla()
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
            ax.set_zlim(*zlim)
            ax.set_xlabel("lon")
            ax.set_ylabel("lat")
            ax.set_zlabel("z")
            ax.yaxis.labelpad=30
            ax.xaxis.labelpad=30
            ax.zaxis.labelpad=30
            ax.set_xticks(ax.get_xticks()[::2])
            ax.set_yticks(ax.get_yticks()[::2])
            idx = np.invert((z[1, :] >= 0) & (mass[1, :] == 1))
            sm = ax.scatter(lon[i, idx], lat[i, idx], z[i, idx], s=mass[i, idx] * k, c = norm(mass[i, idx]), cmap = "plasma_r", norm = mpl.colors.PowerNorm(vmin = 0, vmax = 1, gamma = 0.2))
            ax.set_title(self.time[i])
            ax.view_init(elev=10, azim=45, roll=0)

        frames = len(self.time)

        anim=animation.FuncAnimation(fig, update, blit=False, frames = frames, interval=100)

        if outfile is not None:
            writer = animation.PillowWriter(fps=10)
            anim.save(outfile, writer = writer)

        plt.show()

    def vertical_particle_distribution(self, maxdepth=None, bins=None, maxnum=None):
        """
        Plot the vertical distribution of particles at each time step.
        
        Parameters
        ----------
        maxdepth: float
            Maxdepth of plotting. If None, the maxdepth reached by the particle is used.
        bins: int
            Number of bins for plotting. If None, the number is calculated from length of interval, 
            which, if not configured, will be set to 1 m.
        maxnum: float
            The high limit of bars. If None, the max number in bin is used.
        """

        plt.close()
        fig = plt.figure(figsize = self.figsize)
        
        mainplot = fig.add_axes([.15, .3, .8, .5])
        sliderax = fig.add_axes([.15, .08, .75, .05])

        tslider = Slider(sliderax, 'Timestep', 0, self.steps_output-1,
                         valinit = 0, valfmt='%0.0f', color = "black",
                         valstep = 1)
        
        z = self.o.get_property('z')
        z = np.ma.filled(z, np.nan)

        mass = self.o.get_property("mass")
        mass = np.ma.filled(mass, np.nan)

        if maxdepth is not None:
            maxrange = -np.abs(maxdepth)
        else:
            maxrange = np.min(z)
        
        if bins is not None:
            dz = -maxrange/bins
        else:
            try:
                dz = self.o.get_config('vertical_mixing:verticalresolution')
            except:
                dz = 1.
        
        if maxnum is None:
            hist_series = np.zeros((int(-maxrange/dz), self.steps_output-1))
            bin_series = np.zeros((int(-maxrange/dz)+1, self.steps_output-1))
            for i in range(self.steps_output-1):
                hist_series[:,i], bin_series[:,i] = np.histogram(z[i,:][np.where(mass[i, :] > 0)][np.isfinite(z[i,:][np.where(mass[i, :] > 0)])], bins=int(-maxrange/dz), range=[maxrange, 0])
            maxnum = hist_series.max()
        

        def update(val):
            tindex = int(tslider.val)
            mainplot.cla()
            mainplot.grid()
            mainplot.hist(z[tindex, :][np.where(mass[tindex, :] > 0)], bins=int(-maxrange/dz),
                          range=[maxrange, 0], orientation='horizontal', color = "blue")
            mainplot.set_ylim([maxrange, 0])
            mainplot.set_xlim([0, maxnum])
            mainplot.set_xlabel('number of particles')
            mainplot.set_ylabel('depth [m]')
            fig.canvas.draw_idle()

        update(0)
        tslider.on_changed(update)
        plt.show()

        return fig, mainplot, sliderax, tslider
    
    def vertical_mass_distribution_with_thresholds(self, file_names, thresholds, label = "threshold", maxdepth = None, bins = None, maxnum = None):
        """
        Plot the mass sum in each specified interval of z for given threshold values.

        Parameters
        ----------
        file_names: list of strings
            List of file names to be plotted.
        thresholds: list or ndarray
            List of threshold values, corresponding to the file names.
        label: str
            Label of threshold slider (default "threshold").
        maxdepth: float
            Maxdepth of plotting. If None, the maxdepth reached by the particle is used.
        bins: int
            Number of bins for plotting. If None, the number is calculated from length of interval, 
            which, if not configured, will be set to 1 m.
        maxnum: float
            The high limit of bars. If None, it is calculated as the sum of the initial masses.
        """

        plt.close()
        fig = plt.figure(figsize = self.figsize)

        mainplot = fig.add_axes([.15, .3, .8, .5])
        tsliderax = fig.add_axes([.15, .08, .75, .05])
        thsliderax = fig.add_axes([.15, .01, .75, .05])

        tslider = Slider(tsliderax, 'Timestep', 0, self.steps_output-1,
                         valinit=0, valfmt='%0.0f', valstep = 1, color = "black")
        
        thslider = Slider(thsliderax, label, thresholds[0], thresholds[-1], thresholds[0],
                          valstep = thresholds[1] - thresholds[0], color = "black")

        @jit(nopython = True)
        def mass_dist(mass, z, dz, maxrange):
            dist = []
            zmax = 0
            zmin = -dz
            while zmin >= maxrange:
                if zmin == maxrange:
                    dist.append(np.sum(mass[np.where((z >= zmin) & (z<=zmax))[0]]))
                    break
                dist.append(np.sum(mass[np.where((z > zmin) & (z<=zmax))[0]]))
                zmax = zmin
                zmin -= dz
            return dist[::-1]
        
        @jit(nopython = True)
        def precalculate(mass, z, dz, maxrange):
            dists = []
            for timestep in range(mass.shape[0]):
                dists.append(mass_dist(mass[timestep, :], z[timestep, :], dz, maxrange))
            return dists

        mn = 0
        mr = 0
        ms = []
        zs = []

        for file in file_names:
            o = Open(file)

            z = o.get_property('z')
            z = np.ma.filled(z, np.nan)

            mass = o.get_property("mass")

            mass = np.ma.filled(mass, np.nan)
            mr = min(mr, np.min(z))
            mn = max(mn, np.sum(mass[0,:]))

            ms.append(mass)
            zs.append(z)

        if maxdepth is not None:
            maxrange = -np.abs(maxdepth)
        else:
            maxrange = mr
        
        if maxnum is None:
            maxnum = mn
        
        if bins is not None:
            dz = -maxrange/bins
        else:
            try:
                dz = o.get_config('vertical_mixing:verticalresolution')
            except:
                dz = 1.
        
        dists = []

        for i in range(len(ms)):
            dists.append(precalculate(ms[i], zs[i], dz, maxrange))
          
        yaxis = np.linspace(maxrange+ dz/2, -dz/2, len(dists[0][0]))

        def update(val):
            tindex = int(tslider.val)
            thindex = np.where(thresholds == thslider.val)[0][0]
            mainplot.cla()
            mainplot.grid()
            mainplot.barh(yaxis, dists[thindex][tindex], color = "blue")
            mainplot.set_ylim([maxrange, 0])
            mainplot.set_xlim([0, maxnum])
            mainplot.set_xlabel('mass')
            mainplot.set_ylabel("z")
            mainplot.set_title("Vertical mass distribution", fontweight = "bold")
            fig.canvas.draw_idle()
        
        update(0)

        tslider.on_changed(update)
        thslider.on_changed(update)
        
        plt.show()
        
        return fig, mainplot, tsliderax, tslider, thsliderax, thslider
    
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

    def reset_figure(self):
        fig, ax = plt.subplots(1, 1, figsize = self.figsize)
        self.fig = fig
        self.ax = ax

    def zone_mass_distribution(self, lons, lats, title = None, file2 = None, m0 = None, pathtofile2 = None, read_file = None, write_file = None, labels = ["", ""]):
        plt.close()

        if (file2 is not None) and (read_file is None):
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
                mass = np.loadtxt(str(zone).replace(" ", "") +read_file[0])
                if file2 is not None:
                    mass2 = np.loadtxt(str(zone).replace(" ", "") +read_file[1])
            else:
                mass, m0 = self.zone_crossing_event(lons, lats, zone)
            total_mass = np.sum(mass)

            if file2 is not None and read_file is None:
                mass2, m02 = o2.zone_crossing_event(lons, lats, zone)
                if m0 != m02:
                    raise ValueError("Initial mass of files do not match.")
            if file2 is not None:
                total_mass2 = np.sum(mass2)
            if write_file is not None:
                np.savetxt(str(zone).replace(" ", "") + write_file[0], mass)
                if file2 is not None:
                    np.savetxt(str(zone).replace(" ", "") + write_file[1], mass2)
            
            heights.append(total_mass / m0)
            if file2 is not None:
                heights2.append(total_mass2 / m0)
        
        x = np.arange(3)

        plt.bar(x-0.25, heights, width = 0.4, color = "navy", label = labels[0])

        if file2 is not None:
            plt.bar(x+0.25, heights2, width = 0.4, color = "red", label = labels[1])
        if len(labels) > 0:
            plt.legend()
        if title is not None:
            plt.title(title, fontweight="bold")
        plt.xticks(x, [-200, -1000, "sea floor"])
        plt.xlabel("depth")
        plt.ylabel(r"$\sum m(z = \mathrm{depth}) / m_0$")
        plt.show()
    
    def mass_map(self, lons, lats, title = None, file2 = None, pathtofile2 = None, total_mass_reached = False, temp_file = None, read_file = None, write_file = None):
        plt.close()

        if (file2 is not None) and (read_file is None):
            import os
            curdir = os.getcwd()
            os.chdir(pathtofile2)
            o2 = Plot(file2, self.distribution)
            os.chdir(curdir)
            plt.close()
        
        # fig, ax = plt.subplots(2, 3, subplot_kw={'projection': ccrs.PlateCarree()}, figsize = self.figsize)
        fig = plt.figure(layout="constrained", figsize=self.figsize)
        gs = GridSpec(1, 3, figure=fig)
        minmap, maxmap = 0, 0
        mass_buff = []
        # cmap_buff = []

        for i, zone in enumerate([-200, -1000, "sea floor"]):
            if read_file is not None:
                mass = np.loadtxt(str(zone).replace(" ", "") + read_file)
            else:
                mass, m0 = self.zone_crossing_event(lons, lats, zone)
                total_mass = np.sum(mass)

            if file2 is not None and read_file is None:
                mass2, m02 = o2.zone_crossing_event(lons, lats, zone)
                total_mass2 = np.sum(mass2)
                mass = (np.copy(mass) - mass2) / np.copy(mass)
            if write_file is not None:
                np.savetxt(str(zone).replace(" ", "") + write_file, mass)
            
            m, mid, M = self.get_colormap_midpoint(mass)
            if total_mass_reached:
                print("---------------------\nDepth:", zone)
                print("Total mass 2000: {}\nTotal mass 2003: {}\nTotal initial mass: {}".format(round(total_mass, 2), round(total_mass2, 2), m0))
            # cmap = mpl.cm.RdBu
            # cmap = self.shiftedColorMap(cmap, midpoint = mid, name='shifted{}'.format(i))
            mass_buff.append(mass)
            minmap = min(minmap, m)
            maxmap = max(maxmap, M)
            # cmap_buff.append(cmap)
        cmap = mpl.cm.RdBu
        norm = MidpointNormalize(minmap, maxmap, 0)
        for i, zone in enumerate([-200, -1000, "sea floor"]):
            ax = fig.add_subplot(gs[i], projection= ccrs.PlateCarree())
            ax.add_feature(cartopy.feature.LAND, facecolor="gray",edgecolor='black', zorder = 2)

            sm = ax.contourf(lons, lats, mass_buff[i].T, 60, transform=ccrs.PlateCarree(), cmap = cmap, norm = norm, zorder = 1)
            cb = plt.colorbar(sm, ax = ax, orientation="horizontal")
            cb.set_label(r"$\Delta m/m$")
            cb.ax.set_xticks(cb.get_ticks())
            cb.ax.set_xticklabels(["{:.2f}".format(i) for i in cb.get_ticks()])
            ax.coastlines(zorder = 3)

            if type(zone) == int:
                ax.set_title("depth: " + str(zone) + " m", fontweight = "bold", fontsize=17)
            else:
                ax.set_title("depth: " + zone, fontweight = "bold", fontsize=17)
            
            if total_mass_reached and read_file is None:
                print(m0== m02)
                ax.text(0.5,-0.2, "Total mass 2000: {}\nTotal mass 2003: {}\nTotal initial mass: {}".format(round(total_mass, 2), round(total_mass2, 2), m0), size=15, ha="center", transform=ax.transAxes)
        
        # if temp_file is not None:
        #     from netCDF4 import Dataset
        #     temp = Dataset(temp_file)
        #     ax = fig.add_subplot(gs[1, :3], projection= ccrs.PlateCarree())
        #     ax.add_feature(cartopy.feature.LAND, facecolor="gray",edgecolor='black', zorder = 2)
        #     sm = ax.contourf(temp["lon"][:], temp["lat"][:], temp["thetao"][0, 0, :, :], 60, transform=ccrs.PlateCarree(), cmap = "RdBu", zorder = 1)
        #     cb = plt.colorbar(sm, ax = ax)
        #     cb.set_label(r"$\Delta T\,[\mathrm{^\circ C}]$")
        #     ax.coastlines(zorder = 3)
        #     ax.set_title("Sea surface temperature", fontweight = "bold")
        
        # ax = fig.add_subplot(gs[1, 3:], projection= ccrs.PlateCarree())
        # ax.add_feature(cartopy.feature.LAND, facecolor="gray",edgecolor='black', zorder = 2)
        # sm = ax.contourf(lons, lats, np.ma.filled(self.o.get_property("depth")[0], np.nan)[0, :].reshape(len(lats), len(lons)), 60, transform=ccrs.PlateCarree(), cmap = "RdBu", zorder = 1)
        # cb = plt.colorbar(sm, ax = ax)
        # cb.set_label(r"$z$ [m]")
        # ax.coastlines(zorder = 3)
        # ax.set_title("depth", fontweight = "bold")

        if title is not None:
            fig.suptitle(title, fontweight = "bold")

        plt.show()
    
    def zone_crossing_event(self, lons, lats, depth):
        from global_land_mask import globe

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
        def mass_sum(lons, lats, depth, z, mass, lon, lat):
            mass_at_depth = np.zeros((len(lons), len(lats)))
            m0 = 0
            for i in range(mass.shape[1]):
                drifter_trajectory = z[:, i]
                drifter_mass = mass[:, i]
                if np.isnan(lat[0, i]):
                    continue
                if globe.is_land(lat[0, i], lon[0, i]):
                    continue

                trajectory_nans = np.invert(np.isnan(drifter_trajectory))
                
                if (drifter_trajectory[trajectory_nans]>0).any():
                    continue
                
                if np.any(drifter_mass[1:] == 1):
                    continue

                m0 += drifter_mass[0]
                
                if type(depth) != str:
                    if np.min(drifter_trajectory[trajectory_nans])>depth:
                        continue
                    depth_idx = np.argmin(np.abs(drifter_trajectory[trajectory_nans] - depth))
                
                else:
                    bounce = np.where(drifter_trajectory[:-1] < np.roll(drifter_trajectory, -1)[:-1])[0]
                    if len(bounce) > 0:
                        raise ValueError("Depth decreases.")
                        # depth_idx = bounce[0]
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
                
                lon_grid = np.argmin(np.abs(lons - drifter_lon))
                lat_grid = np.argmin(np.abs(lats - drifter_lat))
                
                mass_at_depth[lon_grid, lat_grid] = mass_at_depth[lon_grid, lat_grid] + drifter_mass[depth_idx]
            return mass_at_depth, m0
        
        return mass_sum(lons, lats, depth, z, mass, lon, lat)

    def mass_z(self, file2 = None):
        plt.close()

        if file2 is not None:
            import os
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
        np.random.seed()
        i = np.random.randint(0, z.shape[1])
        while np.min(z[:, i]) > -1000:
            i = np.random.randint(0, z.shape[1])
        plt.plot(mass[:, i], z[:, i], label = "file1")

        if file2 is not None:
            plt.plot(mass2[:, i], z2[:, i], label = "file2")
        plt.ylim(-200, 0)
        plt.legend()
        plt.show()

    def trajectory_mass_time(self, file2 = None):
        plt.close()

        if file2 is not None:
            import os
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

        moving = self.o.get_property("moving")
        moving = np.ma.filled(moving, np.nan)

        i = 0
        j = 0
        j_sum = 0
        for i in range(moving.shape[1]):
            if np.all(moving[:, i] == 1):
                break
            else:
                if np.sum(moving[:, i])> j_sum:
                    j_sum = np.sum(moving[:, i])
                    j = i
        if not np.all(moving[:, i] == 1):
            i = j
        
        z = z[:, i]
        mass = mass[:, i]

        if file2 is not None:
            mass2 = mass2[:, i]
            z2 = z2[:, i]
        

        fig, [ax1, ax2] = plt.subplots(1, 2)
        
        ax1.plot(mass, z, "-", color = "orange", label = "file1")
        ax2.plot([i.total_seconds() / 3600 for i in self.time], z, "-", color = "orange")

        if file2 is not None:
            ax1.plot(mass2, z2, "-", color = "red", label = "file2")
            ax2.plot([i.total_seconds() / 3600 for i in self.time], z2, "-", color = "red")
        
        ax1.set_xlabel(r"$\frac{m}{m_0}$")
        ax1.set_ylabel("z [m]")
        ax2.set_xlabel("t [h]")

        ax1.legend()
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
    
