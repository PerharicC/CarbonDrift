import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, CheckButtons
from matplotlib import animation
from matplotlib import gridspec

import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
import cartopy
from cartopy import config
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from model.massdecay.carbondrift import *
from model.logger import Logger

log = Logger("CarbonDrift.model.plots")
logger = log.LOGGER

from netCDF4 import Dataset, num2date

from numba import jit
from datetime import datetime, timedelta
import cftime

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

def get_status_info(data):
    status = data.variables["status"]
    meanings = status.flag_meanings.split()
    values = status.flag_values
    # ice = "Ice"
    seafloor = "Reached_Sea_Floor"
    stranded = "stranded"
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
    return sf, st

class Plot:

    def __init__(self, *files, cmap = None, color = None, linestyle = None,
                 lons = None, lats = None, figsize = (20, 20),
                 fontsize = 17, title = None, depth = -200, add = False,
                 diff = True, absolute = False, fontweight = "normal",
                 outfile = None, shrink = 1, clip = None, locations = None, bins = None,
                 loclines = None, prop1 = None, prop2 = None, colorbarlabel = None,
                 xlabel = None, ylabel = None, xlim = None, ylim = None, linewidth = 2,
                 legend = None, areagridpath = f"./supplementary_data/area_grid.npy",
                 biomegridpath = f"./supplementary_data/biomegrid.npy", group = None):
        
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
        if diff and add:
            logger.warning("Both diff and add are set to True, setting add to False.")
            add = False
        
        logger.debug("Importing data.")

        self.objects =  []
        for k in range(len(files)):
            obj = Open(files[k])
            setattr(self, "obj" + str(k + 1), obj)
            self.objects.append(getattr(self, "obj" + str(k + 1)))
        
        logger.debug("Adding attributes.")

        if depth > 0:
            depth*=-1
        if depth <= -5000:
            logger.debug("Setting depth to sea_floor.")
            depth = "sea_floor"
        
        self.depth = depth
        self.cmap = cmap
        self.bins = bins
        if color is None:
            self.color = None
        else:
            self.color = color.split(",")
        if linestyle is None:
            self.linestyle = ["solid", "dashed", "dotted", "dashdot", (0, (3, 5, 1, 5, 1, 5))]
        else:
            self.linestyle = linestyle.split(",")
        
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
        self.diff = diff
        self.add = add
        self.abs = absolute
        self.outfile = outfile
        self.fontweight = fontweight

        logger.debug("Creating time array.")
        self.time = self.obj1.get_time_array()

        logger.debug("Decrypting status numberings.")

        seafloor, stranded = get_status_info(self.obj1.data)
        # self.ice_idx = ice
        self.seafloor_idx = seafloor
        self.stranded_idx = stranded

        logger.debug("Setting up clipping.")
        if clip is None:
            self.clip = False
        else:
            self.clip = True
            self.Vmin, self.Vmax = [float(i.replace("m", "-")) for i in clip.split(":")]

        # if loclines is None:
        #     self.loclines = False
        # else:
        #     self.loclines = loclines
        
        self.areagridpath = areagridpath
        self.biomegridpath = biomegridpath
        self.group = group
        self.loc = locations

        self.prop1 = prop1
        self.prop2 = prop2
        self.cb_units = colorbarlabel
        self.xlabel = xlabel
        self.ylabel = ylabel
        if xlim is not None:
            x1, x2 = xlim.split(":")
            if "m" in x1:
                x1 = -float(x1.strip("m"))
            else:
                x1 = float(x1)
            if "m" in x2:
                x2 = -float(x2.strip("m"))
            else:
                x2 = float(x2)
            self.xlim = [x1, x2]
        else:
            self.xlim = None
        
        if ylim is not None:
            y1, y2 = ylim.split(":")
            if "m" in y1:
                y1 = -float(y1.strip("m"))
            else:
                y1 = float(y1)
            if "m" in y2:
                y2 = -float(y2.strip("m"))
            else:
                y2 = float(y2)
            self.ylim = [y1, y2]
        else:
            self.ylim = None

        self.lw = linewidth

        if legend is None:
            self.legend = False
        else:
            self.labels = legend.split(",")
            self.legend = True

    
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
        """Plot the mass reached at depths [-200m, -1000m, sea_floor] over a cartopy map.

        Parameters
        -----------
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
        
        logger.debug("Searching for bad trajectories.")
        bad = []
        for i in range(len(self.objects)):
            bad.append(self.clean_dataset(self.objects[i]))
        bad = np.ravel(bad)

        logger.debug("Start calculating mass at given depth.")
        mass1 = self.zone_crossing_event(self.obj1, self.lons, self.lats, h, bad)

        if self.clip and not self.diff:
            mass1 = self.clip_array(np.copy(mass1))
        
        logger.debug("Finished calculating mass at given depth.")

        if self.diff:
            logger.debug("Start calculating mass at given depth for file2.")
            mass2 = self.zone_crossing_event(self.obj2, self.lons, self.lats, h, bad)
            logger.debug("Finished calculating mass at given depth for file2.")

            if self.abs:
                mass = np.copy(mass1) - mass2
            else:
                mass = (np.copy(mass1) - mass2) / np.copy(mass1)
            if self.clip:
                mass = self.clip_array(np.copy(mass))
            m, mid, M = self.get_colormap_midpoint(mass)
        elif self.add:
            logger.debug("Start calculating mass at given depth for other files.")
            for obj in self.objects[1:]:
                mass1 += self.zone_crossing_event(obj, self.lons, self.lats, h, bad)
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
        if self.title is not None:
            ax.set_title(self.title, fontweight = self.fontweight)

        ax.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k', facecolor = "beige")
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='k', alpha=0.8, linestyle='--')
        # gl.xlocator = mticker.FixedLocator([-120, -60, 0, 60, 120])
        # gl.ylocator = mticker.FixedLocator([-60, -30, 0, 30, 60])
        # gl.xformatter = LONGITUDE_FORMATTER
        # gl.yformatter = LATITUDE_FORMATTER
        # ax.set_xticks([-120, -60, 0, 60, 120])
        # ax.set_xticklabels([r"$-120^\circ$", r"$-60^\circ$", r"$0^\circ$", r"$60^\circ$", r"$120^\circ$"])
        # ax.set_yticks([-60, -30, 0, 30, 60])
        # ax.set_yticklabels([r"$-60^\circ$", r"$-30^\circ$", r"$0^\circ$", r"$30^\circ$", r"$60^\circ$"])

        plt.tight_layout()

        if self.outfile is not None:
            logger.debug("Saving output file.")
            plt.savefig(self.outfile, dpi = 300, bbox_inches = "tight")
        else:
            plt.show()

    def find_grid_cells_with_no_data(self, obj):
        lon = obj.get_property("lon")
        lon = np.ma.filled(lon, np.nan)[0, :]

        lat = obj.get_property("lat")
        lat = np.ma.filled(lat, np.nan)[0, :]
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
    
    def zone_crossing_event(self, obj:Open, lons, lats, depth, bad):
        
        """Calculate particle properties when crossing a certain depth."""

        logger.debug("Reading required simulation properties.")
        z = obj.get_property('z')
        z = np.ma.filled(z, np.nan)

        mass = obj.get_property("mass")
        mass = np.ma.filled(mass, np.nan)

        lon = obj.get_property("lon")
        lon = np.ma.filled(lon, np.nan)

        lat = obj.get_property("lat")
        lat = np.ma.filled(lat, np.nan)
        status = obj.get_property("status")
        status = np.ma.MaskedArray(status.data, status.mask, float)
        status = np.ma.filled(status, np.nan)

        return self.mass_sum(lons, lats, depth, z, mass, lon, lat, status, self.seafloor_idx, bad)

    # @jit(nopython = True)
    def mass_sum(self, lons, lats, depth, z, mass, lon, lat, status, sfidx, bad):
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
    
    def create_timedelta_array(self, n):
        dt = datetime.strptime(str(self.time[1]),'%Y-%m-%d %H:%M:%S') - datetime.strptime(str(self.time[0]),'%Y-%m-%d %H:%M:%S')
        dt = dt.total_seconds() / 3600
        return np.arange(0, n * dt, dt)

    def clean_dataset(self, obj):
        """Get indicies of trajectories which are on land / do not interact (mass stays one forever) or had a problem with reading the temperature."""

        from global_land_mask import globe

        z = obj.get_property('z')
        z = np.ma.filled(z, np.nan)

        mass = obj.get_property("mass")
        mass = np.ma.filled(mass, np.nan)

        T = obj.get_property("sea_water_temperature")
        T = np.ma.filled(T, np.nan)

        lon = obj.get_property("lon")
        lon = np.ma.filled(lon, np.nan)

        lat = obj.get_property("lat")
        lat = np.ma.filled(lat, np.nan)

        status = obj.get_property("status")
        status = np.ma.MaskedArray(status.data, status.mask, float)
        status = np.ma.filled(status, np.nan)

        bad_trajectories = []

        for i in range(mass.shape[1]):
            drifter_trajectory = z[:, i]
            drifter_mass = mass[:, i]

            if np.isnan(lat[0, i]) or globe.is_land(lat[0, i], lon[0, i]):
                bad_trajectories.append(i)
                continue

            trajectory_nans = np.invert(np.isnan(drifter_trajectory))
            
            if (drifter_trajectory[trajectory_nans]>drifter_trajectory[0]).any() or np.any(drifter_mass[1:] >= drifter_mass[0]):
                bad_trajectories.append(i)
            
            # elif np.any(status[:, i] == self.ice_idx): #Ice
            #     bad_trajectories.append(i)
        return bad_trajectories
    
    def current_strength(self):
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
        rows, cols = np.where(D == 0)
        D[rows, cols] = np.nan
        # rows2, cols2 = np.where(in_cell == 0)
        # D[rows2, cols2] = np.nan
        logger.debug("Start plotting")
        ax.add_feature(cartopy.feature.LAND, zorder=3, edgecolor='k', facecolor = "beige")

        if self.clip:
            D = self.clip_array(np.copy(D))
        
        sm = ax.contourf(self.lons, self.lats, D.T, 20, transform=ccrs.PlateCarree(), cmap = cmap, zorder = 1, extend = "max")
        # ax.contour(self.lons, self.lats, in_cell.T, levels = 0, transform = ccrs.PlateCarree(), colors = "red")
        ax.contour(self.lons, self.lats, is_at_seafloor.T, levels = 1, transform = ccrs.PlateCarree(), colors = "black", label = "seafloor")
        # ax.contour(self.lons, self.lats, is_stranded.T, levels = 1, transform = ccrs.PlateCarree(), colors = "orange", label = "stranded")
        
        # incell_label = Line2D([0], [0], color='red', lw=2)
        seafloor_label = Line2D([0], [0], color='black', lw=2)
        # stranded_label = Line2D([0], [0], color='orange', lw=2)

        # ax.legend([incell_label, seafloor_label, stranded_label], ['moved out of cell', "seafloor", "stranded"], loc='upper right')
        ax.legend([seafloor_label], ["seafloor"], loc="upper right")

        cb = plt.colorbar(sm, ax = ax, orientation="horizontal", shrink = self.shrink)
        if self.cb_units is not None:
            cb.set_label(f"{self.cb_units}")
        if self.title is not None:
            ax.set_title(self.title, fontweight = self.fontweight)

        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='k', alpha=0.8, linestyle='--')
        # gl.xlocator = mticker.FixedLocator([-180, -90, 0, 90, 180])
        # gl.ylocator = mticker.FixedLocator([-90, -45, 0, 45, 90])
        # gl.xformatter = LONGITUDE_FORMATTER
        # gl.yformatter = LATITUDE_FORMATTER
        # ax.set_xticks([-180, -90, 0, 90, 180])
        # ax.set_xticklabels([r"$-180^\circ$", r"$-90^\circ$", r"$0^\circ$", r"$90^\circ$", r"$180^\circ$"])
        # ax.set_yticks([-90, -45, 0, 45, 90])
        # ax.set_yticklabels([r"$-90^\circ$", r"$-45^\circ$", r"$0^\circ$", r"$45^\circ$", r"$90^\circ$"])

        plt.tight_layout()

        if self.outfile is not None:
            logger.debug("Saving output file.")
            plt.savefig(self.outfile, dpi = 300, bbox_inches = "tight")
        else:
            plt.show()
    
    def calculate_horizontal_distance(self, obj: Open, lons, lats, bad):

        from haversine import haversine, Unit
        from tqdm import trange

        logger.debug("Reading required simulation properties.")

        lon = obj.get_property("lon")
        lon = np.ma.filled(lon, np.nan)

        lat = obj.get_property("lat")
        lat = np.ma.filled(lat, np.nan)
        z = obj.get_property("z")
        z = np.ma.filled(z, np.nan)
        status = obj.get_property("status")
        status = np.ma.MaskedArray(status.data, status.mask, float)
        status = np.ma.filled(status, np.nan)
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
                if len(np.where(drifter_status == self.stranded_idx)[0]) > 0:
                    idx = np.where(drifter_status == self.stranded_idx)[0][0]
                    stranded = True
                else:
                    idx = np.where(drifter_status == self.seafloor_idx)[0]
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

        lon = self.obj1.get_property("lon")
        lon = np.ma.filled(lon, np.nan)

        lat = self.obj1.get_property("lat")
        lat = np.ma.filled(lat, np.nan)

        # xticks = []
        # yticks = []
        k = 1

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
            
            ax.plot(lon[0, idx], lat[0, idx], "o", markersize = 10, label = self.labels[k-1] if self.legend else None)
            # if self.loclines:
            #     ax.hlines(lat[0, idx], -180, lon[0, idx], linestyles="dashed", color = "black")
            #     ax.vlines(lon[0, idx], -90, lat[0, idx], linestyle = "dashed", color = "black")
            # xticks.append(lon[0, idx])
            # yticks.append(lat[0, idx])
            k += 1
        
        ax.set_global()
        ax.add_feature(cartopy.feature.LAND, facecolor="gray",edgecolor='black', zorder = 1)
        ax.coastlines(zorder = 2)
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='k', alpha=0.8, linestyle='--')
        # ax.set_xticks(xticks)
        # ax.set_xticklabels([r"${}^\circ $".format(i) for i in xticks] , rotation = 30)
        # ax.set_yticks(yticks)
        # ax.set_yticklabels(r"${}^\circ $".format(i) for i in yticks)
        if self.legend:
            ax.legend()

        plt.tight_layout()

        if self.outfile is not None:
            logger.debug("Saving output file.")
            plt.savefig(self.outfile, dpi = 300, bbox_inches = "tight")
        else:
            plt.show()

    def drifter_properties(self):
        '''Plot drifter prop1 vs prop2 at a specified location.'''
        import matplotlib.style as style
        
        style.use('tableau-colorblind10')
        
        plt.close()

        self.fig, self.ax = plt.subplots(1, 1, figsize = self.figsize)

        if self.loc is None:
            raise AttributeError("Locations are not specified.")
        

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

        for i, j in self.loc:
            for k, obj in enumerate(self.objects):
                lon = obj.get_property("lon")
                lon = np.ma.filled(lon, np.nan)

                lat = obj.get_property("lat")
                lat = np.ma.filled(lat, np.nan)
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
                
                if self.prop1 != "time":
                    x = obj.get_property(self.prop1)
                    x = np.ma.filled(x, np.nan)
                    unitsx = "[" + obj.data[self.prop1].units + "]"
                else:
                    x = self.create_timedelta_array(lon.shape[0])
                    unitsx = "[h]"
                
                if self.prop2 != "time":
                    y = obj.get_property(self.prop2)
                    y = np.ma.filled(y, np.nan)
                    unitsy = "[" + self.obj1.data[self.prop2].units + "]"
                else:
                    y = self.create_timedelta_array(lon.shape[0])
                    unitsy = "[h]"

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
                    color = self.color[k % len(self.color)]
                else:
                    color = None
                linestyle = self.linestyle[k % len(self.linestyle)]
                line, = self.ax.plot(X, Y, lw = self.lw, linestyle = linestyle, color = color)
                lines.append(line)
                
        
        if self.xlabel is None:
            self.ax.set_xlabel(self.prop1 + " " + unitsx)
        else:
            self.ax.set_xlabel(f"{self.xlabel}")
        if self.ylabel is None:
            self.ax.set_ylabel(self.prop2 + " " + unitsy)
        else:
            self.ax.set_ylabel(f"{self.ylabel}")

        if self.title is not None:
            self.ax.set_title(self.title)
        
        if self.legend: self.ax.legend(lines, [f"{label}" for label in self.labels])
        if self.xlim is not None: self.ax.set_xlim(*self.xlim)
        if self.ylim is not None: self.ax.set_ylim(*self.ylim)
        plt.grid()
        plt.tight_layout()

        if self.outfile is not None:
            logger.debug("Saving output file.")
            plt.savefig(self.outfile, dpi = 300, bbox_inches = "tight")
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
        z = obj.get_property('z')
        z = np.ma.filled(z, np.nan)

        mass = obj.get_property("mass")
        mass = np.ma.filled(mass, np.nan)

        lon = obj.get_property("lon")
        lon = np.ma.filled(lon, np.nan)

        lat = obj.get_property("lat")
        lat = np.ma.filled(lat, np.nan)

        biome = obj.get_property("origin_marker")
        biome = np.ma.MaskedArray(biome.data, biome.mask, float)
        biome = np.ma.filled(biome, np.nan)

        status = obj.get_property("status")
        status = np.ma.MaskedArray(status.data, status.mask, float)
        status = np.ma.filled(status, np.nan)

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
                #Sea_floor.
                # bounce = np.where(drifter_trajectory[:-1] < np.roll(drifter_trajectory, -1)[:-1])[0]
                # if len(bounce) > 0:
                #     raise ValueError("Depth decreases.")
                
                s = status[:, i]
                sea_floor = np.where(s == self.seafloor_idx)[0]
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
        z = obj.get_property('z')
        z = np.ma.filled(z, np.nan)

        mass = obj.get_property("mass")
        mass = np.ma.filled(mass, np.nan)

        lon = obj.get_property("lon")
        lon = np.ma.filled(lon, np.nan)

        lat = obj.get_property("lat")
        lat = np.ma.filled(lat, np.nan)

        biome = obj.get_property("origin_marker")
        biome = np.ma.MaskedArray(biome.data, biome.mask, float)
        biome = np.ma.filled(biome, np.nan)

        status = obj.get_property("status")
        status = np.ma.MaskedArray(status.data, status.mask, float)
        status = np.ma.filled(status, np.nan)

        area = self.load_area()

        logger.debug("Start summing masses.")

        depth = self.depth
        mass_by_biomes = [[], [], [], []]
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
                #Sea_floor.
                # bounce = np.where(drifter_trajectory[:-1] < np.roll(drifter_trajectory, -1)[:-1])[0]
                # if len(bounce) > 0:
                #     raise ValueError("Depth decreases.")
                
                s = status[:, i]
                sea_floor = np.where(s == self.seafloor_idx)[0]
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
        colors = ["blue", "orange", "green"]
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
                ax[i%2 + 2*k].bar(xticks, y, bottom = bottom[i%2 + 2*k], color = colors[(i%6) // 2])
                bottom[i%2 + 2*k]+=y
            if i % 6 == 5:
                j+=1
        for i in range(len(ax)):
            if i > 0:
                ax[i].yaxis.set_tick_params(labelleft=False)
            if i % 2 == 0:
                ax[i].set_title("Carcasses", fontsize=15)
            else:
                ax[i].set_title("Fecal Matter", fontsize=15)
            ax[i].tick_params(axis='x', labelrotation=90, labelsize = 15)
            
            fig.add_subplot(ax[i])
        
        if self.ylabel is not None:
            ax[0].set_ylabel(f"{self.ylabel}")
        
        if self.legend:
            lines = [mpl.patches.Patch(facecolor='blue'),
                     mpl.patches.Patch(facecolor='orange'),
                     mpl.patches.Patch(facecolor='green')]
            fig.legend(lines, [f"{label}" for label in self.labels], loc = "center right")
        if self.xlim is not None: ax[0].set_xlim(*self.xlim)
        if self.ylim is not None: ax[0].set_ylim(*self.ylim)
        plt.tight_layout()        
        
        if self.outfile is not None:
            logger.debug("Saving output file.")
            plt.savefig(self.outfile, dpi = 300, bbox_inches = "tight")
        else:
            plt.show()     
    
    def mass_flux_map(self):
        """Plot the mass flux reached at given depth on world map.

        Parameters
        -----------
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
        logger.debug("Searching for bad trajectories.")
        bad = []
        for i in range(len(self.objects)):
            bad.append(self.clean_dataset(self.objects[i]))
        bad = np.ravel(bad)

        logger.debug("Start calculating mass at given depth.")
        flux1 = self.zone_crossing_event(self.obj1, self.lons, self.lats, h, bad) / area

        if self.clip and not self.diff:
            flux1 = self.clip_array(np.copy(flux1))
        
        logger.debug("Finished calculating mass at given depth.")

        if self.diff:
            logger.debug("Start calculating mass at given depth for file2.")
            flux2 = self.zone_crossing_event(self.obj2, self.lons, self.lats, h, bad) / area
            logger.debug("Finished calculating mass at given depth for file2.")

            if self.abs:
                flux = np.copy(flux1) - flux2
            else:
                flux= (np.copy(flux1) - flux2) / np.copy(flux1)
            if self.clip:
                flux = self.clip_array(np.copy(flux))
            m, mid, M = self.get_colormap_midpoint(flux)
        elif self.add:
            logger.debug("Start calculating mass at given depth for other files.")
            for obj in self.objects[1:]:
                flux1 += self.zone_crossing_event(obj, self.lons, self.lats, h, bad) / area
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
        ax.coastlines(zorder = 3, resolution='10m')

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
        if self.title is not None:
            ax.set_title(self.title, fontweight = self.fontweight)

        ax.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k', facecolor = "beige")
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='k', alpha=0.8, linestyle='--')
        # gl.xlocator = mticker.FixedLocator([-120, -60, 0, 60, 120])
        # gl.ylocator = mticker.FixedLocator([-60, -30, 0, 30, 60])
        # gl.xformatter = LONGITUDE_FORMATTER
        # gl.yformatter = LATITUDE_FORMATTER
        # ax.set_xticks([-120, -60, 0, 60, 120])
        # ax.set_xticklabels([r"$-120^\circ$", r"$-60^\circ$", r"$0^\circ$", r"$60^\circ$", r"$120^\circ$"])
        # ax.set_yticks([-60, -30, 0, 30, 60])
        # ax.set_yticklabels([r"$-60^\circ$", r"$-30^\circ$", r"$0^\circ$", r"$30^\circ$", r"$60^\circ$"])

        plt.tight_layout()

        if self.outfile is not None:
            logger.debug("Saving output file.")
            plt.savefig(self.outfile, dpi = 300, bbox_inches = "tight")
        else:
            plt.show()

    def animate_3D(self):
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
        
        logger.debug("Reading required simulation properties.")
        lon = self.obj1.get_property('lon')
        lon = np.ma.filled(lon, np.nan)

        lat = self.obj1.get_property("lat")
        lat = np.ma.filled(lat, np.nan)

        z = self.obj1.get_property('z')
        z = np.ma.filled(z, np.nan)

        mass = self.obj1.get_property("mass")
        mass = np.ma.filled(mass, np.nan)
        
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

        k = 1 / mass[0, :]
        dt = (self.time[1]-self.time[0]).total_seconds()
        idx = np.where(np.invert(np.isnan(mass[1, :])))[0][0]
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
        ax.set_title(self.time[0])
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
            sc = ax.scatter(lon[i, idx], lat[i, idx], z[i, idx], s=mass[i, idx] * k[idx], c = w, cmap = cmap, vmin=0, vmax = w0)
            cb.update_normal(sc)
            ax.set_title(self.time[i])
            ax.view_init(elev=10, azim=45, roll=0)

        frames = len(self.time)

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
        area = np.load(self.areagridpath).T
        logger.debug("Searching for bad trajectories.")
        bad = []
        for i in range(len(self.objects)):
            bad.append(self.clean_dataset(self.objects[i]))
        bad = np.ravel(bad)

        logger.debug("Start calculating mass at given depth.")
        flux1 = self.zone_crossing_event(self.obj1, self.lons, self.lats, self.depth, bad) / area
        if self.add:
            logger.debug("Start calculating mass at given depth for other files.")
            for obj in self.objects[1:]:
                flux1 += self.zone_crossing_event(obj, self.lons, self.lats, self.depth, bad) / area
            logger.debug("Finished calculating mass at given depth for other files.")
        flux = np.copy(flux1)
        flux[np.isnan(flux)] = 0
        flux_mask = flux == 0
        # flux = np.ma.array(data = flux, mask = flux_mask)
        mean_lon_flux = np.ma.mean(flux, axis = 0)
        # print(mean_lon_flux.shape)
        plt.plot(mean_lon_flux, self.lats, color = "black", linewidth = self.lw)
        plt.xlabel(self.xlabel, fontsize = 17)
        plt.ylabel(self.ylabel, fontsize = 17)
        
        plt.tight_layout()

        if self.outfile is not None:
            logger.debug("Saving output file.")
            plt.savefig(self.outfile, dpi = 300, bbox_inches = "tight")
        else:
            plt.show()
    
    def mass_flux_distribution(self):
        area = np.load(self.areagridpath).T
        logger.debug("Searching for bad trajectories.")
        bad = []
        for i in range(len(self.objects)):
            bad.append(self.clean_dataset(self.objects[i]))
        bad = np.ravel(bad)

        logger.debug("Start calculating mass at given depth.")
        flux1 = self.zone_crossing_event(self.obj1, self.lons, self.lats, self.depth, bad) / area
        if self.add:
            logger.debug("Start calculating mass at given depth for other files.")
            for obj in self.objects[1:]:
                flux1 += self.zone_crossing_event(obj, self.lons, self.lats, self.depth, bad) / area
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
            self.ax.legend()
        elif self.group in ["lonmean", "latmean"]:
            if self.color is not None:
                color = self.color[0]
            else:
                color = None
            flux = np.ma.array(data = flux, mask = flux_mask)
            self.ax.hist(np.ma.mean(flux, axis = 0 if self.group == "lonmean" else 1), bins = self.bins, color = color, density=True)
        
        if self.xlabel is None:
            self.ax.set_xlabel("C flux [g C Y^-1 m^-2]")
        else:
            self.ax.set_xlabel(f"{self.xlabel}")
        if self.ylabel is None:
            self.ax.set_ylabel("probability density")
        else:
            self.ax.set_ylabel(f"{self.ylabel}")

        if self.title is not None:
            self.ax.set_title(self.title)

        plt.tight_layout()

        if self.outfile is not None:
            logger.debug("Saving output file.")
            plt.savefig(self.outfile, dpi = 300, bbox_inches = "tight")
        else:
            plt.show()
            

