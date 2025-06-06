import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib import animation
from matplotlib import gridspec
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
from matplotlib.lines import Line2D
import cartopy
from cartopy import config
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

import carbondrift
from carbondrift.models.massdecay.carbondrift import *
from carbondrift.models.logger import Logger

log = Logger("CarbonDrift.model.plots")
logger = log.LOGGER

from netCDF4 import Dataset, num2date

from numba import jit
from datetime import datetime, timedelta
import os
import json

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

    def get_properties_from_list(self, properties):
        p = []
        for property in properties:
            if property in ["status", "origin_marker"]:
                prop = self.get_property(property)
                prop = np.ma.MaskedArray(prop.data, prop.mask, float)
                prop = np.ma.filled(prop, np.nan)
            else:
                prop = self.get_property(property)
                prop = np.ma.filled(prop, np.nan)
            p.append(prop)
        return p

def get_status_info(data):
    """Assign index values to status info."""
    status = data.variables["status"]
    meanings = status.flag_meanings.split()
    values = status.flag_values
    # ice = "Ice"
    seafloor = "Reached_Sea_Floor"
    stranded = "stranded"
    decayed = "Fully_Decayed"
    # ice_idx = meanings.index(ice)
    if seafloor in meanings:
        seafloor_idx = meanings.index(seafloor)
        sf = values[seafloor_idx]
    else:
        sf = -10
    if stranded in meanings:
        stranded_idx = meanings.index(stranded)
        st = values[stranded_idx]
    else:
        st = -11
    if decayed in meanings:
        decayed_idx = meanings.index(decayed)
        sd = values[decayed_idx]
    else:
        sd = -11
    return sf, st, sd

class Plot:
    """CarbonDrift plotting module. Contains many methods for simulation visualization."""

    def __init__(self, *files, cmap = None, color = None, linestyle = None,
                 lons = None, lats = None, figsize = (20, 20), diffidx = None,
                 fontsize = 17, title = None, depth = -200, add = False, suptitle = None,
                 diff = True, absolute = False, fontweight = "normal", martincurve = None,
                 outfile = None, shrink = 1, clip = None, locations = None, bins = None,
                 prop1 = None, prop2 = None, colorbarlabel = None,
                 xlabel = None, ylabel = None, xlim = None, ylim = None, linewidth = 2,
                 legend = None, areagridpath = f"{os.path.dirname(os.path.dirname(carbondrift.__file__))}/supplementary_data/area_grid.npy",
                 biomegridpath = f"{os.path.dirname(os.path.dirname(carbondrift.__file__))}/supplementary_data/biome_grid.npy", group = None,
                 mhwintesitypath = None, dpi = 300):
        """
        Parameters
        -----------
        files: str
            Unlimited number of string simulation paths. At least one is required!
        cmap: str or None
            Matplolib cmap name (default None)
        color: list or None
            List of matplolib color names. E.g ["red","blue","orange"] (default None)
        linestyle: list or None
            List of linestyles. E.g. ["solid","dashed","dotted"] (default None)
        lons: numpy array or None
            Specific longitudes to use in plotting methods. If None the np.array(-180,180,1) is asummed.
        lons: numpy array or None
            Specific latitudes to use in plotting methods. If None the np.array(-90,90,1) is asummed.
        figsize: tuple
            Figure size (default(20, 20))
        diffidx: int or None
            Index of files which separates masses that are added with a plus sign and masses that are
            added with a minus sign. Only use if diff and add are both True. If None len(files) // 2 is asummed.
        fontsize: float
            Fontsize (default 17)
        title: str or None
            Figure title (default None)
        depth: float
            Depth in meters at which a calculation takes place. If positive a depth*(-1) preprocess is performed.
            If depth <= -5000, calculations at seafloor are asummed (default -200)
        add: bool
            Add file outputs (deafult False)
        suptitle: list or None
            List of suptitles (default None)
        diff: bool
            Take difference between two or more files if specified (default True)
        absolute: bool
            Show absolute value of mass/mass flux... (default False)
        fontweight: str
            Fontweight (default "normal")
        martincurve: float or None
            Martin curve coefficient. If specified the curve is plotted (default None)
        outfile: str or None
            Outfile path, if None plt.show() is called (default None)
        shrink: float
            Coorbar shrink (default 1)
        clip: list of length 2
            Len 2 list with min and max value for clipping (default None)
        locations: list of tuples or None
            Locations for dynamics analysis. Example [(longitude1, latitude1), (longitude2, latitude2),...]
            (default None)
        bins: int or None
            Number of bins for histogram (default None)
        prop1: str or None
            Property to plot on x axes. Must be the short name in the .nc file (default None)
        prop2: str or None
            Property to plot on y axes. Must be the short name in the .nc file (default None)
        xlabel: str or None
            xlabel (default None)
        ylabel: str or None
            ylabel (default None)
        xlim: list of length 2
            Len 2 list with min and max value for xlim (default None)
        ylim: list of length 2
            Len 2 list with min and max value for ylim (default None)
        linewidth: int
            linewidth (default 2)
        legend: list or None
            List of legend labels (default None)
        areagridpath: str
            Path to npy file containg area of each 1 by 1 grid cell (default f"./supplementary_data/area_grid.npy")
        biomegridpath: str
            Path to npy file containing biome int labels for each 1 by 1 grid cell (default f"./supplementary_data/biomegrid.npy")
        group: str or None
            Grouping of tracers. Choose between biome, lonmean, latmean or none (default None)
        mhwintensitypath: TODO EXPLAIN
        dpi: int
            Output image resolution (default 300)
        """

        logger.debug("Setting up figure.")
        plt.rcParams.update({'font.size': fontsize})
        fig, ax = plt.subplots(1, 1, figsize = figsize)
        self.fig = fig
        self.ax = ax
        
        if len(files) == 1 and (diff or add):
            logger.warning("Only one filepath given, changing diff and add to False.")
            diff = False
            add = False
        elif len(files) > 1 and not (diff or add):
            logger.warning("Multiple files are given but difference and add are set to False.")
        if diff and add and diffidx is None:
            logger.warning("Both diff and add are set to True, but diffidx is None. Setting diffidx to half the number of files.")
            diffidx = len(files) // 2
        
        logger.debug("Importing data.")

        self.objects =  []
        self.decay_type = []
        self.initial_velocities = []
        self.fragmentation = False
        for k in range(len(files)):
            obj = Open(files[k])
            setattr(self, "obj" + str(k + 1), obj)
            self.objects.append(getattr(self, "obj" + str(k + 1)))
            param_file = files[k].replace(".nc", ".json")
            if os.path.exists(param_file):
                with open(param_file, "r") as f:
                    param = json.load(f)
                    self.decay_type.append(param["decaytype"])
                    self.initial_velocities.append(param["initialvelocity"])
                    self.fragmentation = not param["fragmentation"]
                f.close()
        
        logger.debug("Adding attributes.")

        if depth > 0:
            depth*=-1
        if depth <= -5000:
            logger.debug("Setting depth to sea_floor.")
            depth = "sea_floor"
        
        self.depth = depth
        self.cmap = cmap
        self.bins = bins
        self.color = color
        if linestyle is None:
            self.linestyle = ["solid", "dashed", "dotted", "dashdot", (0, (3, 5, 1, 5, 1, 5))]
        else:
            self.linestyle = linestyle
        
        if lons is None:
            logger.warning("Lons is not given - asumming 1 degree resolution over all space.")
            lons = np.arange(-180, 180, 1)
        if lats is None:
            logger.warning("Lats is not given - asumming 1 degree resolution over all space.")
            lats = np.arange(-90, 90, 1)
        self.lons = lons
        self.lats = lats
        self.figsize = figsize
        self.shrink = shrink
        
        self.title = title
        self.suptitle = suptitle
        self.diff = diff
        self.add = add
        self.diffidx = diffidx
        self.abs = absolute
        self.outfile = outfile
        self.fontweight = fontweight

        logger.debug("Creating time array.")
        self.time = []
        for i in range(len(self.objects)):
            self.time.append(self.objects[i].get_time_array())

        logger.debug("Decrypting status numberings.")
        seafloor, stranded, decayed = [], [], []
        for i in range(len(self.objects)):
            stinfo = get_status_info(self.objects[i].data)
            seafloor.append(stinfo[0])
            stranded.append(stinfo[1])
            decayed.append(stinfo[2])
        # self.ice_idx = ice
        self.seafloor_idx = seafloor
        self.stranded_idx = stranded
        self.decayed_idx = decayed

        logger.debug("Setting up clipping.")
        if clip is None:
            self.clip = False
        else:
            self.clip = True
            self.Vmin, self.Vmax = clip
        
        self.areagridpath = areagridpath
        self.biomegridpath = biomegridpath
        self.mhwIpath = mhwintesitypath
        self.group = group
        self.loc = locations

        self.prop1 = prop1
        self.prop2 = prop2
        if martincurve is not None:
            self.plot_martin_curve = True
            self.martin_curve_coef = float(martincurve)
        else:self.plot_martin_curve = False
        self.cb_units = colorbarlabel
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.xlim = xlim
        self.ylim = ylim

        self.lw = linewidth

        if legend is None:
            self.legend = False
        else:
            self.labels = legend
            self.legend = True

        self.dpi = dpi
    
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
        # mpl.cm.register(cmap=newcmap)

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
    
    def clip_array(self, array):
        m, M = np.min(array[np.invert(np.isnan(array))]), np.max(array[np.invert(np.isnan(array))])
        row, col = np.where(array < self.Vmin)
        array[row, col] = self.Vmin
        row, col = np.where(array > self.Vmax)
        array[row, col] = self.Vmax
        return array

    def mass_map(self):
        """Plot the mass at depth over a cartopy map."""

        plt.close("all")
        fig, ax = plt.subplots(1, 1, figsize=self.figsize, subplot_kw={'projection': ccrs.PlateCarree()})
        
        logger.debug("Initialize colormap.")
        if self.cmap is None:
            if self.diff:
                cmap = mpl.cm.RdBu_r
            else:
                cmap = mpl.cm.Reds
        else:
            try:
                cmap = getattr(mpl.cm, self.cmap)
            except AttributeError:
                logger.debug("Specified colormap doesn't exist, changing to Reds.")
                cmap = "Reds"
        
        h = self.depth
        
        logger.debug("Searching for bad trajectories.")
        bad = []
        for i in range(len(self.objects)):
            bad.append(self.clean_dataset(self.objects[i]))
        bad = np.ravel(bad)

        logger.debug("Start calculating mass at given depth.")
        mass1 = self.zone_crossing_event(self.obj1, self.lons, self.lats, h, bad,self.seafloor_idx[0])

        if self.clip and not self.diff and not self.add:
            mass1 = self.clip_array(np.copy(mass1))
        
        logger.debug("Finished calculating mass at given depth.")

        if self.diff:
            if not self.add:
                logger.debug("Start calculating mass at given depth for file2.")
                mass2 = self.zone_crossing_event(self.obj2, self.lons, self.lats, h, bad,self.seafloor_idx[1])
                logger.debug("Finished calculating mass at given depth for file2.")

                if self.abs:
                    mass = np.copy(mass1) - mass2
                else:
                    mass = (np.copy(mass1) - mass2) / np.copy(mass1)
                if self.clip:
                    mass = self.clip_array(np.copy(mass))
                m, mid, M = self.get_colormap_midpoint(mass)
            else:
                mass2 = np.zeros(mass1.shape)
                for i in range(1, self.diffidx):
                    mass1 += self.zone_crossing_event(self.objects[i], self.lons, self.lats, h, bad,self.seafloor_idx[i])
                for i in range(self.diffidx, len(self.objects)):
                    mass2 += self.zone_crossing_event(self.objects[i], self.lons, self.lats, h, bad,self.seafloor_idx[i])
                if self.abs:
                    mass = np.copy(mass1) - mass2
                else:
                    mass = (np.copy(mass1) - mass2) / np.copy(mass1)
                if self.clip:
                    mass = self.clip_array(np.copy(mass))
                m, mid, M = self.get_colormap_midpoint(mass)
        elif self.add:
            logger.debug("Start calculating mass at given depth for other files.")
            for i, obj in enumerate(self.objects[1:]):
                mass1 += self.zone_crossing_event(obj, self.lons, self.lats, h, bad,self.seafloor_idx[i])
            logger.debug("Finished calculating mass at given depth for other files.")
            mass = np.copy(mass1)
            if self.clip:
                mass = self.clip_array(np.copy(mass))
        if self.diff or self.add:
            NODATA = self.find_grid_cells_with_no_data(self.obj1)
            mass[np.isnan(NODATA)] = np.nan
        else:
            NODATA = self.find_grid_cells_with_no_data(self.obj1)
            mass1[np.isnan(NODATA)] = np.nan
        logger.debug("Start plotting")
        ax.coastlines(zorder = 3, resolution='10m')

        extend = 'both' if self.clip else None

        if not self.diff and not self.add:
            sm = ax.contourf(self.lons, self.lats, mass1.T, 20, transform=ccrs.PlateCarree(), cmap = cmap, zorder = 1,extend=extend, extendfrac='auto')
            cb = plt.colorbar(sm, ax = ax, orientation="horizontal", shrink = self.shrink)
        else:
            if self.diff:
                cmap = self.shiftedColorMap(cmap, midpoint = mid, name='shifted')
            sm = ax.contourf(self.lons, self.lats, mass.T, 20, transform=ccrs.PlateCarree(), cmap = cmap, zorder = 1,extend = extend, extendfrac='auto')
            cb = plt.colorbar(sm, ax = ax, orientation="horizontal", shrink = self.shrink)
        if self.cb_units is not None:
            cb.set_label(f"{self.cb_units}")
        self.axis_setup(ax)

        ax.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k', facecolor = "beige")
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='k', alpha=0.8, linestyle='--')

        plt.tight_layout()

        if self.outfile is not None:
            logger.debug("Saving output file.")
            plt.savefig(self.outfile, dpi = self.dpi, bbox_inches = "tight")
        else:
            plt.show()

    def find_grid_cells_with_no_data(self, obj:Open):
        lon, lat = obj.get_properties_from_list(["lon", "lat"])
        lon, lat = lon[0, :], lat[0, :]
        X = np.zeros((len(self.lons), len(self.lats)))
        for i in range(len(self.lats)):
            for j in range(len(self.lons)):
                if not np.any((lon == self.lons[j]) & (lat == self.lats[i])):
                    X[j, i] = np.nan
        return X
    
    @staticmethod        
    def interpolate(x, y, x0):
        if np.any(x == x0):
            idx2 = np.where(x == x0)[0][0]
            return y[idx2], idx2
        idx2 = np.where(x < x0)[0][0]
        y1, y2 = y[idx2 - 1], y[idx2]
        x1, x2 = x[idx2 - 1], x[idx2]
        k = (y2 - y1) / (x2 - x1)
        n = y2 - k * x2
        return k * x0 + n, idx2
    
    def zone_crossing_event(self, obj:Open, lons, lats, depth, bad, seafloor_idx):
        
        """Calculate particle properties when crossing a certain depth."""

        logger.debug("Reading required simulation properties.")
        z, mass, lon, lat, status = obj.get_properties_from_list([
                                        "z", "mass", "lon", "lat", "status"
                                        ])

        return self.mass_sum(lons, lats, depth, z, mass, lon, lat, status, seafloor_idx, bad)

    # @jit(nopython = True)
    def mass_sum(self, lons, lats, depth, z, mass, lon, lat, status, sfidx, bad):
        '''
        Return mass at each grid cell at given depth.
        '''
        logger.debug("Start summing masses.")
        mass_at_depth = np.zeros((len(lons), len(lats)))
        for i in range(mass.shape[1]):

            #Remove bad particles.
            if i in bad:
                continue

            drifter_trajectory = z[:, i]
            drifter_mass = mass[:, i]
            m0 = drifter_mass[0]

            trajectory_nans = np.invert(np.isnan(drifter_trajectory))
            
            if type(depth) != str:
                if np.min(drifter_trajectory[trajectory_nans])>depth:
                    continue

                m, depth_idx = self.interpolate(drifter_trajectory[trajectory_nans], drifter_mass[trajectory_nans], depth)
            
            else:
                #Sea_floor.
                # bounce = np.where(drifter_trajectory[:-1] < np.roll(drifter_trajectory, -1)[:-1])[0]
                # if len(bounce) > 0:
                #     raise ValueError("Depth decreases.")
                
                s = status[:, i]
                sea_floor = np.where(s == sfidx)[0]
                if len(sea_floor) >0:
                    depth_idx = sea_floor[0]
                else:
                    depth_idx = len(drifter_trajectory)-1
                m = drifter_mass[depth_idx]


            if np.isnan(m):
                continue

            drifter_lon = lon[depth_idx, i]
            drifter_lat = lat[depth_idx, i]

            if np.isnan(drifter_lon) or np.isnan(drifter_lat):
                continue
            
            #Find closest position on grid.
            lon_grid = np.argmin(np.abs(lons - drifter_lon))
            lat_grid = np.argmin(np.abs(lats - drifter_lat))
            if not self.diff and not self.abs and not self.add:
                m/=m0
            mass_at_depth[lon_grid, lat_grid] = mass_at_depth[lon_grid, lat_grid] + m
        logger.debug("Finish summing masses.")
        return mass_at_depth
    
    def create_timedelta_array(self, n, objidx):
        dt = datetime.strptime(str(self.time[objidx][1]),'%Y-%m-%d %H:%M:%S') - datetime.strptime(str(self.time[objidx][0]),'%Y-%m-%d %H:%M:%S')
        dt = dt.total_seconds() / 3600
        return np.arange(0, n * dt, dt)

    def clean_dataset(self, obj:Open):
        """Get indicies of trajectories which are on land / have any other problems."""

        z, mass, T, lon, lat, status = obj.get_properties_from_list([
            "z", "mass", "sea_water_temperature", "lon", "lat", "status"
        ])

        bad_trajectories = []

        for i in range(mass.shape[1]):
            drifter_trajectory = z[:, i]
            drifter_mass = mass[:, i]

            if np.isnan(lat[0, i]):
                bad_trajectories.append(i)
                continue

            trajectory_nans = np.invert(np.isnan(drifter_trajectory))
            
            if (drifter_trajectory[trajectory_nans]>drifter_trajectory[0]).any() or np.any(drifter_mass[1:] >= drifter_mass[0]):
                bad_trajectories.append(i)
        return bad_trajectories
    
    def current_strength(self):
        """Plot horizontal distance traveled between intial position and position at some depth, for ecah grid cell."""
        plt.close()
        fig, ax = plt.subplots(1, 1, figsize=self.figsize, subplot_kw={'projection': ccrs.PlateCarree()})
        
        logger.debug("Initialize colormap.")
        if self.cmap is None:
            cmap = "Reds"
        else:
            try:
                cmap = getattr(mpl.cm, self.cmap)
            except AttributeError:
                logger.debug("Specified colormap doesn't exist, changing to Reds.")
                cmap = "Reds"
        
        logger.debug("Searching for bad trajectories.")
        bad = self.clean_dataset(self.obj1)

        
        logger.debug("Start calculating distance.")
        distance, in_cell, is_at_seafloor, is_stranded = self.calculate_horizontal_distance(self.obj1, self.lons, self.lats, bad)
        
        D = np.copy(distance)
        rows, cols = np.where(D == 0) #Land
        D[rows, cols] = np.nan
        # rows2, cols2 = np.where(in_cell == 0)
        # D[rows2, cols2] = np.nan
        logger.debug("Start plotting")
        ax.add_feature(cartopy.feature.LAND, zorder=2, edgecolor='k', facecolor = "beige")

        if self.clip:
            D = self.clip_array(np.copy(D))
        
        sm = ax.contourf(self.lons, self.lats, D.T, 20, transform=ccrs.PlateCarree(), cmap = cmap, zorder = 1, extend = "max")
        # ax.contour(self.lons, self.lats, in_cell.T, levels = 0, transform = ccrs.PlateCarree(), colors = "red")
        ax.contour(self.lons, self.lats, is_at_seafloor.T, levels = 1, transform = ccrs.PlateCarree(), colors = "red", label = "seafloor")
        # ax.contour(self.lons, self.lats, is_stranded.T, levels = 1, transform = ccrs.PlateCarree(), colors = "orange", label = "stranded")
        
        # incell_label = Line2D([0], [0], color='red', lw=2)
        seafloor_label = Line2D([0], [0], color='red', lw=2)
        # stranded_label = Line2D([0], [0], color='orange', lw=2)

        # ax.legend([incell_label, seafloor_label, stranded_label], ['moved out of cell', "seafloor", "stranded"], loc='upper right')
        ax.legend([seafloor_label], ["seafloor"], loc="upper right")

        cb = plt.colorbar(sm, ax = ax, orientation="horizontal", shrink = self.shrink)
        if self.cb_units is not None:
            cb.set_label(f"{self.cb_units}")
        self.axis_setup(ax)

        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='k', alpha=0.8, linestyle='--', zorder = 100)

        plt.tight_layout()

        if self.outfile is not None:
            logger.debug("Saving output file.")
            plt.savefig(self.outfile, dpi = self.dpi, bbox_inches = "tight")
        else:
            plt.show()
    
    def calculate_horizontal_distance(self, obj: Open, lons, lats, bad):
        '''
        Calculate horizontal distanace between initial position and position at given depth for each grid cell.
        
        Returns
        -------
        distance traveled
            ndarray of size (lon,lat) with distances
        is_in_cell
            boolean ndarray of size (lon, lat) of drifters which are still in the same cell as in the beegining
        is_at_seafloor
            boolean ndarray of size (lon, lat) of drifters which have reached the seafloor
        is_stranded
            boolean ndarray of size (lon, lat) of drifters which were stranded on the castline'''

        from haversine import haversine, Unit
        from tqdm import trange

        logger.debug("Reading required simulation properties.")

        lon, lat, z, status = obj.get_properties_from_list([
            "lon", "lat", "z", "status"
        ])

        distance_traveled = np.zeros((len(lons), len(lats)))
        is_in_cell = np.ones((len(lons), len(lats)))
        is_stranded = np.zeros((len(lons), len(lats)))
        is_at_seafloor = np.zeros((len(lons), len(lats)))
        for i in trange(lon.shape[1]):
            stranded = False
            reached_sea_floor = False
            #Remove bad particles.
            if i in bad:
                # progress.update(1)
                continue

            drifter_lon = lon[:, i]
            drifter_lat = lat[:, i]
            drifter_z = z[:, i]
            drifter_status = status[:, i]
            nans = np.invert(np.isnan(drifter_lat))
            if np.min(drifter_z[nans]) > self.depth:
                if len(np.where(drifter_status == self.stranded_idx[0])[0]) > 0:
                    idx = np.where(drifter_status == self.stranded_idx[0])[0][0]
                    stranded = True
                elif len(np.where(drifter_status == self.decayed_idx[0])[0]) > 0:
                    idx = np.where(drifter_status == self.decayed_idx[0])[0][0]
                    stranded = True
                else:
                    idx = np.where(drifter_status == self.seafloor_idx[0])[0]
                    if len(idx) == 0:
                        logger.warning("Drifter at initial position ({}, {}) has not been deactivated/passed the specified depth.".format(drifter_lat[0], drifter_lon[0]))
                        logger.warning("Taking the last position.")
                        idx = len(drifter_z[nans]) - 1
                    else:
                        idx = idx[0]
                        reached_sea_floor = True
                phi = drifter_lat[idx]
                lam = drifter_lon[idx]
            else:
                phi, idx = self.interpolate(drifter_z, drifter_lat, self.depth)
                lam, _ = self.interpolate(drifter_z, drifter_lon, self.depth)
            
            pos1 = (drifter_lat[0], drifter_lon[0])
            pos2 = (phi, lam)

            lon_grid = np.argmin(np.abs(lons - pos1[1]))
            lat_grid = np.argmin(np.abs(lats - pos1[0]))
            if np.abs(pos2[0] - pos1[0]) > (lats[1] - lats[0]) / 2 or np.abs(pos1[1] - pos2[1]) > (lons[1] - lons[0]) / 2:
                is_in_cell[lon_grid, lat_grid] = 0
                d = haversine(pos1, pos2, unit = Unit.KILOMETERS)
                distance_traveled[lon_grid, lat_grid] = d
            else:
                d = haversine(pos1, pos2, unit = Unit.KILOMETERS)
                distance_traveled[lon_grid, lat_grid] = d
            if stranded:
                is_stranded[lon_grid, lat_grid] = 1
            if reached_sea_floor:
                is_at_seafloor[lon_grid, lat_grid] = 1
        return distance_traveled, is_in_cell, is_at_seafloor, is_stranded

    def drifter_locations(self):
        '''Plot locations of specified drifters on a world map.'''
        from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
        import matplotlib.style as style

        style.use('seaborn-v0_8-bright')

        if self.loc is None:
            raise AttributeError("Locations are not specified.")
        
        plt.close()
        fig, ax = plt.subplots(1, 1, subplot_kw={'projection': ccrs.PlateCarree()}, figsize = self.figsize)

        logger.debug("Searching for bad trajectories.")
        bad = self.clean_dataset(self.obj1)
        logger.debug("Reading required simulation properties.")

        lon, lat = self.obj1.get_properties_from_list(["lon", "lat"])

        k = 0
        markers = "s^vox"

        for i, j in self.loc:
            lon_idx = np.where(lon[0, :] == i)
            lat_idx = np.where(lat[0, :] == j)
            idx = np.intersect1d(lon_idx, lat_idx)
            if len(idx) == 0:
                logger.info("Could not find drifter with starting position ({}, {}).".format(i, j))
                continue
            else:
                idx = idx[0]
            
            if idx in bad:
                logger.info("({}, {}) is a bad trajectory and will not be included.".format(i, j))
                continue
            if self.color is not None:
                if len(self.color) >1:
                    if k >= len(self.color):
                        color = None
                    else:
                        color = self.color[k]
                else:
                    color = self.color[0]
            else:
                color = None

            ax.plot(lon[0, idx], lat[0, idx], markers[k%len(markers)], color = color, markersize = 10, label = self.labels[k] if self.legend else None, zorder = 100)
            k += 1
        
        ax.set_global()
        ax.add_feature(cartopy.feature.LAND, facecolor="beige",edgecolor='black', zorder = 1)
        ax.coastlines(zorder = 2)
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='k', alpha=0.8, linestyle='--')
        self.axis_setup(ax)
        plt.tight_layout()

        if self.outfile is not None:
            logger.debug("Saving output file.")
            plt.savefig(self.outfile, dpi = self.dpi, bbox_inches = "tight")
        else:
            plt.show()

    def martin_curve(self, z, flux_100 = None):
        if flux_100 is None:
            return (np.abs(z) / 100) ** (-self.martin_curve_coef)
        else: return flux_100 * (np.abs(z) / 100) ** (-self.martin_curve_coef)

    def drifter_properties(self):
        '''Plot drifter prop1 vs prop2 at a specified location.'''
        import matplotlib.style as style
        
        style.use('tableau-colorblind10')
        
        plt.close()

        if self.loc is None:
            raise AttributeError("Locations are not specified.")
        
        self.fig, self.ax = plt.subplots(1, len(self.loc), figsize = self.figsize)
        try:
            self.ax[0]
        except TypeError:
            self.ax = [self.ax]
        logger.debug("Searching for bad trajectories.")
        bad = []
        for i in range(len(self.objects)):
            bad.append(self.clean_dataset(self.objects[i]))
        logger.debug("Reading required simulation properties.")

        if self.prop1 is None:
            logger.info("Property 1 not specified. Setting x to time.")
            self.prop1 = "time"
        if self.prop2 is None:
            logger.info("Property 2 not specified. Setting y to mass.")
            self.prop2 = "mass"
                
        lines = []

        m = 0
        for i, j in self.loc:
            for k, obj in enumerate(self.objects):
                lon, lat = obj.get_properties_from_list(["lon", "lat"])
                lon_idx = np.where(lon[0, :] == i)
                lat_idx = np.where(lat[0, :] == j)
                idx = np.intersect1d(lon_idx, lat_idx)

                if len(idx) == 0:
                    logger.info("Could not find drifter with starting position ({}, {}).".format(i, j))
                    continue
                else:
                    idx = idx[0]
            
                if idx in bad[k]:
                    logger.info("({}, {}) is a bad trajectory and will not be included.".format(i, j))
                    continue
                
                if self.prop1 == "total_horizontal_velocity":
                    u = obj.get_property("x_sea_water_velocity")
                    v = obj.get_property("y_sea_water_velocity")
                    x = np.sqrt(u ** 2 + v ** 2)
                    x = np.ma.filled(x, np.nan)
                elif self.prop1 == "tau":
                    T = obj.get_property("sea_water_temperature")
                    T = np.ma.filled(T, np.nan)
                    T[T< 0]=0
                    decay = self.decay_type[k]
                    if decay == "linear":
                        decayrate = 0.064 * T + 0.02
                    else:
                        decayrate = 0.140 * np.exp(0.145 * T)
                    x = 1 / decayrate
                    x *= 24
                elif self.prop1 != "time":
                    x = obj.get_property(self.prop1)
                    x = np.ma.filled(x, np.nan)
                else:
                    x = self.create_timedelta_array(lon.shape[0], k)
                
                if self.prop2 == "total_horizontal_velocity":
                    u = obj.get_property("x_sea_water_velocity")
                    v = obj.get_property("y_sea_water_velocity")
                    y = np.sqrt(u ** 2 + v ** 2)
                    Y = np.ma.filled(x, np.nan)
                elif self.prop2 == "tau":
                    T = obj.get_property("sea_water_temperature")
                    T = np.ma.filled(T, np.nan)
                    T[T< 0]=0
                    decay = self.decay_type[k]
                    if decay == "linear":
                        decayrate = 0.064 * T + 0.02
                    else:
                        decayrate = 0.140 * np.exp(0.145 * T)
                    y = 1 / decayrate
                    Y *= 24
                elif self.prop2 != "time":
                    y = obj.get_property(self.prop2)
                    y = np.ma.filled(y, np.nan)
                else:
                    y = self.create_timedelta_array(lon.shape[0], k)

                if self.prop1 == "time":
                    X = x
                else:
                    X = x[:, idx]
                    if not self.abs and self.prop1 == "mass":
                        X /= X[0]
            
                if self.prop2 == "time":
                    Y = y
                else:
                    Y = y[:, idx]
                    if not self.abs and self.prop2 == "mass":
                        Y /= Y[0]
                
                if self.color is not None:
                    try:
                        color = self.color[k % len(self.color)]
                    except IndexError:
                        logger.warning(f"Setting color to None for object number {k}", exc_info=IndexError)
                        color = None
                else:
                    color = None
                linestyle = self.linestyle[k % len(self.linestyle)]
                line, = self.ax[m].plot(X, Y, lw = self.lw, linestyle = linestyle, color = color)
                lines.append(line)
        
            if self.plot_martin_curve:
                if self.abs:
                    #ADD FLUXES
                    pass
                else:
                    line, = self.ax[m].plot(self.martin_curve(Y), Y, lw = self.lw, linestyle = "dotted")
                if m ==0 and self.legend:
                    lines.append(line)
                    self.labels.append("Martin Curve")

            self.axis_setup(self.ax[m], False, True)
            if self.suptitle is not None:
                self.ax[m].set_title(self.suptitle[m], fontweight = self.fontweight)

            # self.ax[m].tick_params(labelbottom=False,labeltop=True)
            self.ax[m].xaxis.set_ticks_position('top')
            self.ax[m].xaxis.set_label_position('top')
            m+=1
        if self.title is not None:
            self.fig.suptitle(self.title, fontweight = self.fontweight)
        if self.legend: self.ax[0].legend(lines, [f"{label}" for label in self.labels])
        plt.tight_layout()

        if self.outfile is not None:
            logger.debug("Saving output file.")
            plt.savefig(self.outfile, dpi = self.dpi, bbox_inches = "tight")
        else:
            plt.show()

    def get_mass_sum_at_depth(self):

        '''Return sum of mass over all grid points at given depth.'''

        M = self.zone_crossing_event(self.obj1, self.lons, self.lats, self.depth, self.clean_dataset(self.obj1))
        if not self.abs:
            x = len(np.where(M > 0)[0])
        else:
            x = 1
        return np.sum(M[np.invert(np.isnan(M))]) / x
    
    def get_biome_weighted_mass_at_depth(self, obj = None):
        '''Return sum of mass over all grid points in the 4 biomes at given depth.'''
        logger.debug("Reading required simulation properties.")
        if obj is None:
            obj = self.obj1
        z, mass, lon, lat, biome, status = obj.get_properties_from_list([
            "z", "mass", "lon", "lat", "origin_marker", "status"
        ])

        logger.debug("Start summing masses.")

        depth = self.depth
        mass_by_biomes = np.zeros(4)
        for i in range(mass.shape[1]):

            drifter_trajectory = z[:, i]
            drifter_mass = mass[:, i]
            b = int(biome[0, i])

            trajectory_nans = np.invert(np.isnan(drifter_trajectory))
            
            if type(depth) != str:
                if np.min(drifter_trajectory[trajectory_nans])>depth:
                    continue

                m, depth_idx = self.interpolate(drifter_trajectory[trajectory_nans], drifter_mass[trajectory_nans], depth)
            
            else:
                
                s = status[:, i]
                sea_floor = np.where(s == self.seafloor_idx[0])[0]
                if len(sea_floor) >0:
                    depth_idx = sea_floor[0]
                else:
                    depth_idx = len(drifter_trajectory)-1
                m = drifter_mass[depth_idx]


            if np.isnan(m):
                continue
            
            mass_by_biomes[b] += m
        
        logger.debug("Finish summing masses.")
        return mass_by_biomes#, np.sum(mass_by_biomes) * 10 ** (-15)
    
    def get_mean_biome_mass_flux_at_depth(self, obj = None):
        '''Return sum of mass flux (per area) over all grid points in the 4 biomes at a given depth.'''
        logger.debug("Reading required simulation properties.")
        if obj is None:
            obj = self.obj1
        z, mass, lon, lat, biome, status = obj.get_properties_from_list([
            "z", "mass", "lon", "lat", "origin_marker", "status"
        ])

        area = self.load_area()

        logger.debug("Start summing masses.")

        depth = self.depth
        mass_by_biomes = [[], [], [], []]
        for i in range(mass.shape[1]):

            drifter_trajectory = z[:, i]
            drifter_mass = mass[:, i]
            if drifter_mass[0] <=0:continue
            b = int(biome[0, i])

            trajectory_nans = np.invert(np.isnan(drifter_trajectory))
            
            if type(depth) != str:
                if np.min(drifter_trajectory[trajectory_nans])>depth:
                    continue

                m, depth_idx = self.interpolate(drifter_trajectory[trajectory_nans], drifter_mass[trajectory_nans], depth)
            
            else:
                #Sea_floor.
                # bounce = np.where(drifter_trajectory[:-1] < np.roll(drifter_trajectory, -1)[:-1])[0]
                # if len(bounce) > 0:
                #     raise ValueError("Depth decreases.")
                
                s = status[:, i]
                sea_floor = np.where(s == self.seafloor_idx[0])[0]
                if len(sea_floor) >0:
                    depth_idx = sea_floor[0]
                else:
                    depth_idx = len(drifter_trajectory)-1
                m = drifter_mass[depth_idx]


            if np.isnan(m):
                continue
            A = area[np.where(self.lats == lat[0, i])[0][0], np.where(self.lons == lon[0, i])[0][0]]
            if A == 0:
                continue
            mass_by_biomes[b].append(m / A)
        
        logger.debug("Finish summing masses.")
        return np.asarray([np.mean(i) for i in mass_by_biomes])#, np.sum(mass_by_biomes) * 10 ** (-15)
    
    def load_area(self):
        return np.load(self.areagridpath)

    def mean_export_biome_flux(self):
        '''
        Bar plot of mean fluxes of multiple groups (e.g. Cnidaria egestion, Ctenophora mortality...) for each biome.
        '''
        logger.debug("Setting up figure.")
        plt.close()
        fig = plt.figure(figsize = self.figsize)
        outer_grid = gridspec.GridSpec(1, 4, wspace=0.05)
        ax = []
        biome_titles = ["HCSS", "LC", "HCPS", "COAST"]
        for i in range(4):
            ax_outer = fig.add_subplot(outer_grid[i])
            ax_outer.axis('off')
            ax_outer.set_title(biome_titles[i], y = 1.05)
            innergrid = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer_grid[i], wspace=0, hspace=0)
            for j in range(2):
                if not (i == 0 and j == 0):
                    ax.append(plt.Subplot(fig, innergrid[j], sharey=ax[0]))
                else:
                    ax.append(plt.Subplot(fig, innergrid[j]))
        logger.debug("Figure set.")
        bottom = np.zeros((len(ax), len(self.objects) // 8+1))
        if self.color is None:
            colors = ["blue", "orange", "green"]
        else:
            colors = self.color
        if self.xlabel is not None:
            xticks = self.xlabel.split(",")
        else:
            xticks = np.arange(len(self.objects) // 8)
        j = 0
        for i, obj in enumerate(self.objects):
            flux = self.get_mean_biome_mass_flux_at_depth(obj)

            for k in range(4):
                y = np.zeros(len(xticks))
                y[j] = flux[k]
                ax[i%2 + 2*k].bar(xticks, y, bottom = bottom[i%2 + 2*k], color = colors[(i%6) // 2], alpha = 0.7)
                bottom[i%2 + 2*k]+=y
            if i % 6 == 5:
                j+=1
        if self.suptitle is None:
            suptitle = ["Carcasses", "Fecal Matter"]
        else:
            suptitle = self.suptitle
        for i in range(len(ax)):
            if i > 0:
                ax[i].yaxis.set_tick_params(labelleft=False)
            if i % 2 == 0:
                ax[i].set_title(suptitle[0], fontsize=15, fontweight=self.fontweight)
            else:
                ax[i].set_title(suptitle[1], fontsize=15, fontweight=self.fontweight)
            ax[i].tick_params(axis='x', labelrotation=90, labelsize = 15)
            
            fig.add_subplot(ax[i])
        
        if self.ylabel is not None:
            ax[0].set_ylabel(f"{self.ylabel}")
        
        if self.legend:
            lines = [mpl.patches.Patch(facecolor=colors[0]),
                     mpl.patches.Patch(facecolor=colors[1]),
                     mpl.patches.Patch(facecolor=colors[2])]
            fig.legend(lines, [f"{label}" for label in self.labels], loc = "center right")
        if self.xlim is not None: ax[0].set_xlim(*self.xlim)
        if self.ylim is not None: ax[0].set_ylim(*self.ylim)
        plt.tight_layout()        
        
        if self.outfile is not None:
            logger.debug("Saving output file.")
            plt.savefig(self.outfile, dpi = self.dpi, bbox_inches = "tight")
        else:
            plt.show()     
    
    def mass_flux_map(self):
        """Plot the mass flux reached at given depth on world map.
        """

        plt.close("all")
        fig, ax = plt.subplots(1, 1, figsize=self.figsize, subplot_kw={'projection': ccrs.PlateCarree()})
        
        logger.debug("Initialize colormap.")
        if self.cmap is None:
            if self.diff:
                cmap = mpl.cm.RdBu_r
            else:
                cmap = mpl.cm.Reds
        else:
            try:
                cmap = getattr(mpl.cm, self.cmap)
            except AttributeError:
                logger.debug("Specified colormap doesn't exist, changing to Reds.")
                cmap = "Reds"
        
        h = self.depth
        area = np.load(self.areagridpath).T
        i = np.where(np.arange(-180, 180, 1) == np.min(self.lons))[0][0]
        j = np.where(np.arange(-180, 180, 1) == np.max(self.lons))[0][0]
        k = np.where(np.arange(-90, 90, 1) == np.min(self.lats))[0][0]
        l = np.where(np.arange(-90, 90, 1) == np.max(self.lats))[0][0]
        area = area[i:j+1, k:l+1]
        logger.debug("Searching for bad trajectories.")
        bad = []
        for i in range(len(self.objects)):
            bad.append(self.clean_dataset(self.objects[i]))
            break

        logger.debug("Start calculating mass at given depth.")
        flux1 = self.zone_crossing_event(self.obj1, self.lons, self.lats, h, bad, self.seafloor_idx[0]) / area

        if self.clip and not self.diff and not self.add:
            flux1 = self.clip_array(np.copy(flux1))
        
        logger.debug("Finished calculating mass at given depth.")

        if self.diff:
            if not self.add:
                logger.debug("Start calculating mass at given depth for file2.")
                flux2 = self.zone_crossing_event(self.obj2, self.lons, self.lats, h, bad,self.seafloor_idx[1]) / area
                logger.debug("Finished calculating mass at given depth for file2.")
                
                if self.abs:
                    flux = np.copy(flux1) - flux2
                else:
                    flux= (np.copy(flux1) - flux2) / np.copy(flux1)
                if self.clip:
                    flux = self.clip_array(np.copy(flux))
                m, mid, M = self.get_colormap_midpoint(flux)
            else:
                flux2 = np.zeros(flux1.shape)
                for i in range(1, self.diffidx):
                    flux1 += self.zone_crossing_event(self.objects[i], self.lons, self.lats, h, bad, self.seafloor_idx[i]) / area
                for i in range(self.diffidx, len(self.objects)):
                    flux2 += self.zone_crossing_event(self.objects[i], self.lons, self.lats, h, bad, self.seafloor_idx[i]) / area
                if self.abs:
                    flux = np.copy(flux1) - flux2
                else:
                    flux= (np.copy(flux1) - flux2) / np.copy(flux1)
                if self.clip:
                    flux = self.clip_array(np.copy(flux))
                m, mid, M = self.get_colormap_midpoint(flux)
        elif self.add:
            logger.debug("Start calculating mass at given depth for other files.")
            for i, obj in enumerate(self.objects[1:]):
                flux1 += self.zone_crossing_event(obj, self.lons, self.lats, h, bad, self.seafloor_idx[i]) / area
            logger.debug("Finished calculating mass at given depth for other files.")
            flux = np.copy(flux1)
            if self.clip:
                flux = self.clip_array(np.copy(flux))
        if self.diff or self.add:
            NODATA = self.find_grid_cells_with_no_data(self.obj1)
            flux[np.isnan(NODATA)] = np.nan
        else:
            NODATA = self.find_grid_cells_with_no_data(self.obj1)
            flux1[np.isnan(NODATA)] = np.nan
        logger.debug("Start plotting")
        ax.coastlines(zorder = 3, resolution='110m')

        extend = 'both' if self.clip else None

        if not self.diff and not self.add:
            sm = ax.contourf(self.lons, self.lats, flux1.T, 20, transform=ccrs.PlateCarree(), cmap = cmap, zorder = 1,extend=extend, extendfrac='auto')
            cb = plt.colorbar(sm, ax = ax, orientation="horizontal", shrink = self.shrink)
        else:
            if self.diff:
                cmap = self.shiftedColorMap(cmap, midpoint = mid, name='shifted')
            sm = ax.contourf(self.lons, self.lats, flux.T, 20, transform=ccrs.PlateCarree(), cmap = cmap, zorder = 1,extend = extend, extendfrac='auto')
            cb = plt.colorbar(sm, ax = ax, orientation="horizontal", shrink = self.shrink)
        if self.cb_units is not None:
            cb.set_label(f"{self.cb_units}")
        self.axis_setup(ax)

        ax.add_feature(cartopy.feature.LAND, zorder=2, edgecolor='k', facecolor = "beige")
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='k', alpha=0.8, linestyle='--', zorder = 100)

        plt.tight_layout()

        if self.outfile is not None:
            logger.debug("Saving output file.")
            plt.savefig(self.outfile, dpi = self.dpi, bbox_inches = "tight")
        else:
            plt.show()

    def animate_3D(self):
        """A 3d animation of the simulation."""
        
        plt.close()

        fig = plt.figure(figsize=self.figsize)
        ax = fig.add_subplot(projection='3d')
        
        logger.debug("Reading required simulation properties.")
        lon, lat, z, mass = self.obj1.get_properties_from_list([
            "lon", "lat", "z", "mass"
        ])
        
        logger.debug("Initialize colormap.")
        if self.cmap is None:
            logger.warning("Colormap hasn't been provided, using jet.")
            cmap = mpl.cm.jet
        else:
            try:
                cmap = getattr(mpl.cm, self.cmap)
            except AttributeError:
                logger.debug("Specified colormap doesn't exist, changing to jet.")
                cmap = "jet"

        if self.fragmentation:
            k = 100
        else:
            k = 1 / mass[0, :]
        dt = (self.time[0][1]-self.time[0][0]).total_seconds()
        idx = np.where(np.invert(np.logical_or(np.isnan(mass[1, :]), np.isnan(z[1, :]))))[0][0]
        try:
            w0 = abs(self.initial_velocities[0]) * 24 * 3600
        except Exception:
            w0 = abs((z[1, idx] - z[0, idx]) * (mass[0, idx]/mass[1, idx]) ** (1/6) / dt)*24*3600
        sc = ax.scatter(lon[0, :], lat[0, :], z[0, :], s=mass[0, :] * k, c = w0 * np.ones(len(mass[0, :])), cmap = cmap, vmin=0, vmax = w0)

        
        logger.debug("Setting up axes properties.")
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
        ax.set_title(self.time[0][0])
        ax.view_init(elev=10, azim=45, roll=0)

        cb = plt.colorbar(sc, shrink = self.shrink)
        cb.set_label(r"$|w(t)|\,\mathrm{[m\,day^{-1}]}$")

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
            idx = np.invert(np.isnan(lon[i, :]))
            w = abs((z[i - 1, idx] - z[i, idx]) / dt)*24*3600
            if self.fragmentation:
                q = k
            else:
                q = k[idx]
            sc = ax.scatter(lon[i, idx], lat[i, idx], z[i, idx], s=mass[i, idx] * q, c = w, cmap = cmap, vmin=0, vmax = w0)
            cb.update_normal(sc)
            ax.set_title(self.time[0][i])
            ax.view_init(elev=10, azim=45, roll=0)

        frames = len(self.time[0])

        logger.debug("Initializing animation.")
        anim=animation.FuncAnimation(fig, update, blit=False, frames = frames, interval=100)

        if self.outfile is not None:
            logger.debug("Adding PillowWriter to animation.")
            writer = animation.PillowWriter(fps=10)
            logger.debug("Saving animation.")
            anim.save(self.outfile, writer = writer)
        else:
            plt.show()

    def mean_lon_mass_flux(self):
        """Plot latitude vs longitude averaged mass fluxes"""
        area = np.load(self.areagridpath).T
        logger.debug("Searching for bad trajectories.")
        bad = self.clean_dataset(self.objects[0])
        bad = np.ravel(bad)

        logger.debug("Start calculating mass at given depth.")
        if self.legend is not None:
            model_num = len(self.labels)
        else:
            model_num = 1
        sim_num = len(self.objects) // model_num
        for mnum in range(model_num):
            flux1 = self.zone_crossing_event(self.objects[mnum*sim_num], self.lons, self.lats, self.depth, bad,self.seafloor_idx[mnum*sim_num]) / area
            if self.add:
                logger.debug("Start calculating mass at given depth for other files.")
                for i, obj in enumerate(self.objects[mnum*sim_num + 1:mnum*sim_num+sim_num]):
                    flux1 += self.zone_crossing_event(obj, self.lons, self.lats, self.depth, bad, self.seafloor_idx[mnum*sim_num + 1 + i]) / area
                logger.debug("Finished calculating mass at given depth for other files.")
            flux = np.copy(flux1)
            flux[np.isnan(flux)] = 0
            flux_mask = flux == 0
            flux = np.ma.array(data = flux, mask = flux_mask)
            mean_lon_flux = np.ma.mean(flux, axis = 0)
            # print(mean_lon_flux.shape)
            if self.color is not None:
                try:
                    color = self.color[mnum]
                except IndexError:
                    color = "black"
                    logger.warning(f"Setting color to None.", exc_info=IndexError)
            else:
                color = "black"
            if self.linestyle is not None:
                try:
                    ls = self.linestyle[mnum]
                except IndexError:
                    ls = "solid"
                    logger.warning(f"Setting linestyle to solid.", exc_info=IndexError)
            else:
                ls = "solid"
            if self.legend is not None:
                label = self.labels[mnum]
            else:
                label = None
            plt.plot(self.lats, mean_lon_flux, color = color, linewidth = self.lw, label = label, ls = ls)
        if self.xlabel is not None: plt.xlabel(f"{self.xlabel}")
        if self.ylabel is not None: plt.ylabel(f"{self.ylabel}")
        if self.xlim is not None:
            plt.xlim(*self.xlim)
        else:
            plt.xlim((-90, 90))
            plt.xticks([-90, -45, 0, 45, 90], [r"$-90^\circ$", r"$-45^\circ$",r"$0^\circ$", r"$45^\circ$", r"$90^\circ$"])
        if self.ylim is not None: plt.ylim(*self.ylim)
        if self.legend is not None: plt.legend()
        plt.grid()
        if self.title is not None:
            plt.title(self.title)
        plt.tight_layout()

        if self.outfile is not None:
            logger.debug("Saving output file.")
            plt.savefig(self.outfile, dpi = self.dpi, bbox_inches = "tight")
        else:
            plt.show()
    
    def mass_flux_distribution(self):
        """Mass flux distributions at some depth for different grouping of tarcers."""
        area = np.load(self.areagridpath).T
        logger.debug("Searching for bad trajectories.")
        bad = []
        for i in range(len(self.objects)):
            bad.append(self.clean_dataset(self.objects[i]))
        bad = np.ravel(bad)

        logger.debug("Start calculating mass at given depth.")
        flux1 = self.zone_crossing_event(self.obj1, self.lons, self.lats, self.depth, bad,self.seafloor_idx[0]) / area
        if self.add:
            logger.debug("Start calculating mass at given depth for other files.")
            for i, obj in enumerate(self.objects[1:]):
                flux1 += self.zone_crossing_event(obj, self.lons, self.lats, self.depth, bad,self.seafloor_idx[i]) / area
            logger.debug("Finished calculating mass at given depth for other files.")
        flux = np.copy(flux1)
        flux[np.isnan(flux)] = 0
        flux_mask = flux == 0
        if self.group is None or self.group == "none":
            if self.color is not None:
                color = self.color[0]
            else:
                color = None
            self.ax.hist(flux[np.invert(flux_mask)], bins=self.bins, color = color, density = True)
        elif self.group == "biome":
            biome_names = ["HCSS", "LC", "HCPS", "COAST"]
            biomes = np.load(self.biomegridpath).T
            for i in range(4):
                if self.color is not None:
                    color = self.color[i % len(self.color)]
                else:
                    color = None
                mask = np.invert(flux_mask) & (biomes == i)
                self.ax.hist(flux[mask], bins = self.bins, color = color, alpha = 0.8, label = biome_names[i], density=True)
        elif self.group in ["lonmean", "latmean"]:
            if self.color is not None:
                color = self.color[0]
            else:
                color = None
            flux = np.ma.array(data = flux, mask = flux_mask)
            self.ax.hist(np.ma.mean(flux, axis = 0 if self.group == "lonmean" else 1), bins = self.bins, color = color, density=True)
        self.ax.axvline(np.ma.mean(flux[mask]), color = "black", lw = self.lw, linestyle = "dashed", label = f"povprečje; {round(np.ma.mean(flux[mask]), 2)}")
        self.ax.legend()
        self.axis_setup(self.ax)

        plt.tight_layout()

        if self.outfile is not None:
            logger.debug("Saving output file.")
            plt.savefig(self.outfile, dpi = self.dpi, bbox_inches = "tight")
        else:
            plt.show()
            
    def mhw_correlations(self):
        logger.debug("Searching for bad trajectories.")
        bad = []
        for i in range(len(self.objects)):
            bad.append(self.clean_dataset(self.objects[i]))
        logger.debug("Reading required simulation properties.")
        h = self.depth
        area = np.load(self.areagridpath).T
        I = np.load(self.mhwIpath).T
        flux1 = np.zeros((360, 180))
        flux2 = np.zeros((360, 180))
        for i in range(self.diffidx):
            flux1 += self.zone_crossing_event(self.objects[i], self.lons, self.lats, h, bad,self.seafloor_idx[i]) / area
        for i in range(self.diffidx, len(self.objects)):
            flux2 += self.zone_crossing_event(self.objects[i], self.lons, self.lats, h, bad,self.seafloor_idx[i]) / area
        flux = (np.copy(flux1) - flux2) / np.copy(flux1)
        self.ax.scatter(I.ravel(), flux.ravel(), color = "k")
        mask = np.invert((flux1 == 0) | (flux2==0) | np.isnan(flux1) | np.isnan(flux2) | np.isnan(I))
        if self.color is not None:
            color = self.color[0]
        else:
            color = "k"
        self.ax.scatter(I.ravel()[mask.ravel()], flux.ravel()[mask.ravel()], color = color)
        cov = np.cov(I.ravel()[mask.ravel()], flux.ravel()[mask.ravel()])
        pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
        ell_radius_x = np.sqrt(1 + pearson)
        ell_radius_y = np.sqrt(1 - pearson)
        ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2, edgecolor = "red", facecolor = 'none', linewidth = 3)
        scale_x = np.sqrt(cov[0, 0]) * 5
        mean_x = np.mean(I.ravel()[mask.ravel()])

        # calculating the standard deviation of y ...
        scale_y = np.sqrt(cov[1, 1]) * 5
        mean_y = np.mean(flux.ravel()[mask.ravel()])

        transf = transforms.Affine2D()
        transf.rotate_deg(45)
        transf.scale(scale_x, scale_y)
        transf.translate(mean_x, mean_y)

        ellipse.set_transform(transf + self.ax.transData)
        self.ax.add_patch(ellipse)
        self.axis_setup(self.ax, None, True)
        plt.tight_layout()

        if self.outfile is not None:
            logger.debug("Saving output file.")
            plt.savefig(self.outfile, dpi = self.dpi, bbox_inches = "tight")
        else:
            plt.show()

    def get_mass_over_mhw_area(self, obj = None):
        '''Return sum of mass over all grid points in the 4 biomes at given depth.'''
        logger.debug("Reading required simulation properties.")
        if obj is None:
            obj = self.obj1
        loc = self.loc
        if loc is None:
            raise AttributeError("Missing mhw Area. Locations must be added.")
        lonmin = min([i[0] for i in loc])
        lonmax = max([i[0] for i in loc])
        latmin = min([i[1] for i in loc])
        latmax = max([i[1] for i in loc])
        z, mass, lon, lat, biome, status = obj.get_properties_from_list([
            "z", "mass", "lon", "lat", "origin_marker", "status"
        ])

        logger.debug("Start summing masses.")

        depth = self.depth
        mass_by_biomes = 0
        for i in range(mass.shape[1]):

            drifter_trajectory = z[:, i]
            drifter_mass = mass[:, i]
            if lon[0, i]<lonmin or lon[0, i]>lonmax or lat[0, i]<latmin or lat[0, i]>latmax:
                continue

            trajectory_nans = np.invert(np.isnan(drifter_trajectory))
            
            if type(depth) != str:
                if np.min(drifter_trajectory[trajectory_nans])>depth:
                    continue

                m, depth_idx = self.interpolate(drifter_trajectory[trajectory_nans], drifter_mass[trajectory_nans], depth)
            
            else:
                s = status[:, i]
                sea_floor = np.where(s == self.seafloor_idx)[0]
                if len(sea_floor) >0:
                    depth_idx = sea_floor[0]
                else:
                    depth_idx = len(drifter_trajectory)-1
                m = drifter_mass[depth_idx]


            if np.isnan(m):
                continue
            
            mass_by_biomes += m
        
        logger.debug("Finish summing masses.")
        return mass_by_biomes * 10 ** (-15)

    def animate_current_3D(self):
        """A 3d and 2d current animation of the simulation.
        Sizes of particles are relative to their initial mass."""
        
        plt.close()

        fig = plt.figure(figsize=self.figsize)
        ax = fig.add_subplot(121, projection='3d')
        ax2 = fig.add_subplot(122, projection = ccrs.PlateCarree())
        
        logger.debug("Reading required simulation properties.")
        z, mass, lon, lat, status = self.obj1.get_properties_from_list([
            "z", "mass", "lon", "lat", "status"
        ])

        seafloor_tracker = {"lon":[], "lat":[], "z":[], "m":[]}
        u = self.obj1.get_property("x_sea_water_velocity")
        v = self.obj1.get_property("y_sea_water_velocity")
        U = np.sqrt(u ** 2 + v ** 2)
        U = np.ma.filled(U, np.nan)
        
        logger.debug("Initialize colormap.")
        if self.cmap is None:
            logger.warning("Colormap hasn't been provided, using jet.")
            cmap = mpl.cm.jet
        else:
            try:
                cmap = getattr(mpl.cm, self.cmap)
            except AttributeError:
                logger.debug("Specified colormap doesn't exist, changing to jet.")
                cmap = "jet"

        dt = (self.time[0][1]-self.time[0][0]).total_seconds()
        idx = np.where(np.invert(np.logical_or(np.isnan(mass[1, :]), np.isnan(z[1, :]))))[0][0]
        try:
            w0 = abs(self.initial_velocities[0]) * 24 * 3600
        except Exception:
            w0 = abs((z[1, idx] - z[0, idx]) * (mass[0, idx]/mass[1, idx]) ** (1/6) / dt)*24*3600
        
        lonlim = np.where((lon[0,:]>=self.xlim[0]) & (lon[0,:]<=self.xlim[1]))[0]
        latlim = np.where((lat[0,:]>=self.ylim[0]) & (lat[0,:]<=self.ylim[1]))[0]
        intersectlim = np.intersect1d(lonlim, latlim)
        xlim = [np.min(lon[:, intersectlim][np.invert(np.isnan(lon[:, intersectlim]))]),
                np.max(lon[:, intersectlim][np.invert(np.isnan(lon[:, intersectlim]))])]
        ylim = [np.min(lat[:, intersectlim][np.invert(np.isnan(lat[:, intersectlim]))]),
                np.max(lat[:, intersectlim][np.invert(np.isnan(lat[:, intersectlim]))])]
        zlim = [np.min(z[:, intersectlim][np.invert(np.isnan(z[:, intersectlim]))]),
                np.max(z[:, intersectlim][np.invert(np.isnan(z[:, intersectlim]))])]
        U0 = np.max(U[:, intersectlim][np.invert(np.isnan(U[:, intersectlim]))])
        zmin = zlim[0]
        k = 1 / mass[0, intersectlim]*100
        sc = ax.scatter(lon[0, intersectlim], lat[0, intersectlim], z[0, intersectlim],
                        s=mass[0, intersectlim] * k, c = U[0, intersectlim], cmap = cmap, vmin=0, vmax = U0)
        sc2 = ax2.scatter(lon[0, intersectlim], lat[0, intersectlim], s=mass[0, intersectlim] * k,
                           c=z[0, intersectlim], cmap = cmap, vmin=zmin, vmax = 0)
        ax2.quiver(lon[0, :], lat[0, :], u[0, :], v[0, :], color='b', units='xy', scale=1/5)
        # ax.scatter(lon[0, intersectlim], lat[0, intersectlim], z[0, intersectlim], s=mass[0, intersectlim] * k, c = U[0, intersectlim], cmap = cmap, vmin=0, vmax = U0)
        logger.debug("Setting up axes properties.")
        ax.yaxis.labelpad=30
        ax.xaxis.labelpad=30
        ax.zaxis.labelpad=30
        ax.set_xlabel("lon")
        ax.set_ylabel("lat")
        ax.set_zlabel("z")
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_zlim(*zlim)
        fig.suptitle(self.time[0][0])
        ax.view_init(elev = 8, azim=45, roll=0)
        
        cb = plt.colorbar(sc, shrink = self.shrink, orientation = "horizontal")
        cb.set_label(r"$|\vec{u}_H(t)|\,\mathrm{[m\,s^{-1}]}$")

        ax2.add_feature(cartopy.feature.LAND, zorder=2, edgecolor='k', facecolor = "beige")
        gl = ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=2, color='k', alpha=0.8, linestyle='--', zorder = 100)
        ax2.set_xlabel("lon")
        ax2.set_ylabel("lat")
        ax2.set_xlim(*xlim)
        ax2.set_ylim(*ylim)
        ax2.grid()

        cb2 = plt.colorbar(sc2, shrink = self.shrink, orientation = "horizontal")
        cb2.set_label(r"$z\,\mathrm{[m]}$")

        def update(i):
            ax.cla()
            ax2.cla()
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
            ax.set_zlim(*zlim)
            ax2.set_xlim(*xlim)
            ax2.set_ylim(*ylim)
            ax.yaxis.labelpad=30
            ax.xaxis.labelpad=30
            ax.zaxis.labelpad=30
            ax.set_xlabel("lon")
            ax.set_ylabel("lat")
            ax.set_zlabel("z")
            ax2.set_xlabel("lon")
            ax2.set_ylabel("lat")
            idx = np.where(np.invert(np.isnan(lon[i, :])))[0]
            idx = np.intersect1d(idx, intersectlim)
            sf = np.where(status[i, :] == self.seafloor_idx)[0]
            sf = np.intersect1d(sf, intersectlim)
            for da_tracer in sf:
                seafloor_tracker["lon"].append(lon[i, da_tracer])
                seafloor_tracker["lat"].append(lat[i, da_tracer])
                seafloor_tracker["z"].append(z[i, da_tracer])
                seafloor_tracker["m"].append(mass[i, da_tracer]/mass[0, da_tracer]*100)
            k = 1 / mass[0, idx]*100
            sc = ax.scatter(lon[i, idx], lat[i, idx], z[i, idx], s=mass[i, idx] * k, c = U[i, idx], cmap = cmap, vmin=0, vmax = U0)
            sc2 = ax2.scatter(lon[i, idx], lat[i, idx], s=mass[i, idx] * k, c=z[i, idx], cmap = cmap, vmin=zmin, vmax = 0)
            ax2.quiver(lon[i, idx], lat[i, idx], u[i, idx], v[i, idx], color='b', units='xy', scale=1/5)
            ax.scatter(seafloor_tracker["lon"], seafloor_tracker["lat"], seafloor_tracker["z"], color = "black", s=seafloor_tracker["m"])
            ax2.scatter(seafloor_tracker["lon"], seafloor_tracker["lat"], color = "black", s=seafloor_tracker["m"])
            ax2.grid()
            ax2.add_feature(cartopy.feature.LAND, zorder=2, edgecolor='k', facecolor = "beige")
            gl = ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=2, color='k', alpha=0.8, linestyle='--', zorder = 100)
            cb.update_normal(sc)
            cb2.update_normal(sc2)
            fig.suptitle(self.time[0][i])
            ax.view_init(elev=8, azim=45, roll=0)

        frames = len(self.time[0])

        logger.debug("Initializing animation.")
        anim=animation.FuncAnimation(fig, update, blit=False, frames = frames, interval=100)

        if self.outfile is not None:
            logger.debug("Adding PillowWriter to animation.")
            writer = animation.PillowWriter(fps=10)
            logger.debug("Saving animation.")
            anim.save(self.outfile, writer = writer)
        else:
            plt.show()

    def vertical_particle_distribution(self):
        """Create an interactive slider vertical distribution of particle count or save it as an animation."""

        plt.close()
        fig = plt.figure(figsize = self.figsize)
        
        if self.outfile is not None:
            logger.warning(
                "Outfile is given, the plot will be saved as an animation. If you want slider control, pelase remove the outfile argument."
                )
            mainplot = fig.add_subplot(111)
        else:
            logger.warning(
                "Outfile is not given. If you want to save the plot as an animation, please add the outfile argument."
                )
            mainplot = fig.add_axes([.15, .3, .8, .5])
            sliderax = fig.add_axes([.15, .08, .75, .05])

        maxdepth = np.min(self.obj1.get_property("depth")) if self.depth is None else self.depth

        z_dist = np.linspace(maxdepth, 0, self.bins)

        logger.debug("Retrieving time information.")
        maxtime = 0
        time_steps = []
        for obj in self.objects:
            t = obj.data["time"].shape[0]
            maxtime = max(t, maxtime)
            time_step = datetime.strptime(obj.data.time_step_output, "%H:%M:%S")
            time_steps.append(timedelta(hours = time_step.hour, minutes=time_step.minute, seconds=time_step.second).total_seconds() / 3600)
        time_steps = np.asarray(np.asarray(time_steps) // np.min(time_steps), dtype=np.int64)
        
        if self.outfile is None:
            logger.debug("Setting up slider.")
            tslider = Slider(sliderax, 'time', 0, maxtime-1,
                            valinit = 0, valfmt='%0.0f', color = "black",
                            valstep = 1)
        
        logger.debug("Precalculating distribution.")
        particle_count = np.zeros((maxtime, len(self.objects), self.bins))
        sea_floor_particles = np.zeros((maxtime, self.bins))
        
        for i, obj in enumerate(self.objects):
            z, status = obj.get_properties_from_list([
                "z", "status"
            ])
            time_step = time_steps[i]
            for t in range(obj.data["time"].shape[0]):
                active = np.invert(np.isnan(z[t, :]))
                hist = np.histogram(z[t, active], bins = self.bins,
                                                    range=[maxdepth, 0])
                particle_count[t*time_step:(t+1)*time_step, i, :] = hist[0]
                sea_floor_particles = self.update_sea_floor_particle_count(sea_floor_particles, self.bins,z, active,
                                                                           status, hist[1], t, self.seafloor_idx[i], time_step)
        
        if not self.abs:
            totnum = np.sum(particle_count[0, :, :])
            particle_count /= totnum
            sea_floor_particles /= totnum
        
        if self.outfile is not None:
            mainplot.barh(z_dist, sea_floor_particles[0, :], color = "black", height = maxdepth/self.bins, label = "seafloor")
            for i in range(particle_count.shape[1]):
                color = self.color[i] if self.color is not None else None
                label = self.labels[i] if self.legend else None
                mainplot.barh(z_dist, particle_count[0, i, :],
                             left = sea_floor_particles[0, :] + np.sum(particle_count[0, :i, :], axis=0),
                             height = maxdepth/self.bins, color = color, label = label)
            self.axis_setup(mainplot, "lower right", True)
        
        def update(val):
            tindex = int(tslider.val) if self.outfile is None else val
            mainplot.cla()
            mainplot.barh(z_dist, sea_floor_particles[tindex, :], color = "black", height = maxdepth/self.bins, label = "seafloor")
            for i in range(particle_count.shape[1]):
                color = self.color[i] if self.color is not None else None
                label = self.labels[i] if self.legend else None
                mainplot.barh(z_dist, particle_count[tindex, i, :],
                             left = sea_floor_particles[tindex, :] + np.sum(particle_count[tindex, :i, :], axis=0),
                             height = maxdepth/self.bins, color = color, label = label)
            if self.outfile is None: fig.canvas.draw_idle()
            self.axis_setup(mainplot, "lower right", True)

        if self.outfile is not None:
            logger.debug("Initializing animation.")
            anim=animation.FuncAnimation(fig, update, blit=False, frames = maxtime, interval=100)
            logger.debug("Adding PillowWriter to animation.")
            writer = animation.PillowWriter(fps=20)
            logger.debug("Saving animation.")
            anim.save(self.outfile, writer = writer)
        else:
            update(0)
            tslider.on_changed(update)
            plt.show()
    
    #Fast computation of seafloor vertical particle distribution.
    @staticmethod
    @jit(nopython=True)
    def update_sea_floor_particle_count(sea_floor_particles, bins, z, active, status, dist, t, sf, time_step):
        for j in range(bins):
            if j == 0:
                zj = z[t, active]<= dist[j]
            else:
                zj = np.logical_and(z[t, active]<= dist[j], z[t, active]> dist[j-1])
            sea_floor_particles[t*time_step:, j] += len(np.where(status[t, active][zj] == sf)[0])
        return sea_floor_particles
    
    def axis_setup(self, ax, legendloc = None, grid=False):
        if legendloc!=False:
            if self.legend:ax.legend(loc = legendloc)
        if self.xlim is not None:ax.set_xlim(*self.xlim)
        if self.ylim is not None: ax.set_ylim(*self.ylim)
        if self.xlabel is not None: ax.set_xlabel(f"{self.xlabel}")
        if self.ylabel is not None: ax.set_ylabel(f"{self.ylabel}")
        if self.title is not None: ax.set_title(self.title, fontweight = self.fontweight)
        if grid:ax.grid()