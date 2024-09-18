import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, CheckButtons
from matplotlib import animation
from matplotlib.gridspec import GridSpec

import matplotlib.ticker as mticker

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
    # ice_idx = meanings.index(ice)
    seafloor_idx = meanings.index(seafloor)
    return values[seafloor_idx]#, values[ice_idx] 

class Plot:

    def __init__(self, file1, file2 = None, file3 = None, file4 = None, cmap = None,
                 lons = None, lats = None, figsize = (20, 20),
                 fontsize = 17, title = None, depth = -200,
                 diff = True, absolute = False, fontweight = "normal",
                 outfile = None, shrink = 1, clip = None, locations = None,
                 loclines = None, prop1 = None, prop2 = None, colorbarlabel = None,
                 xlabel = None, ylabel = None, xlim = None, ylim = None, linewidth = 2,
                 legend = None):

        logger.debug("Setting up figure.")
        fig, ax = plt.subplots(1, 1, figsize = figsize)
        self.fig = fig
        self.ax = ax
        plt.rcParams.update({'font.size': fontsize})
        
        if file2 is None and diff:
            logger.warning("Only one filepath given, changing diff to False.")
            diff = False
        elif file2 is not None and not diff:
            logger.warning("File2 is given but difference is set to False. Plotting file1.")
        
        logger.debug("Importing data.")

        self.obj = Open(file1)

        self.objects = [self.obj]
        if file2 is not None:
            self.obj2 = Open(file2)
            self.objects.append(self.obj2)
        if file3 is not None:
            self.obj3 = Open(file3)
            self.objects.append(self.obj3)
        if file4 is not None:
            self.obj4 = Open(file4)
            self.objects.append(self.obj4)
        
        logger.debug("Adding attributes.")

        if depth > 0:
            depth*=-1
        if depth <= -5000:
            logger.debug("Setting depth to sea_floor.")
            depth = "sea_floor"
        
        self.depth = depth
        self.cmap = cmap
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
        self.abs = absolute
        self.outfile = outfile
        self.fontweight = fontweight

        logger.debug("Creating time array.")
        self.time = self.obj.get_time_array()

        logger.debug("Decrypting status numberings.")

        seafloor = get_status_info(self.obj.data)
        # self.ice_idx = ice
        self.seafloor_idx = seafloor

        logger.debug("Setting up clipping.")
        if clip is None:
            self.clip = False
        else:
            self.clip = True
            self.Vmin, self.Vmax = [float(i.replace("m", "-")) for i in clip.split(":")]

        if loclines is None:
            self.loclines = False
        else:
            self.loclines = loclines
        
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
        bad1 = self.clean_dataset(self.obj)
        bad2 = []
        if self.diff: bad2 = self.clean_dataset(self.obj2)
        
        logger.debug("Start calculating mass at given depth.")
        mass1 = self.zone_crossing_event(self.obj, self.lons, self.lats, h, bad1, bad2)

        if self.clip and not self.diff:
            mass1 = self.clip_array(np.copy(mass1))
        
        logger.debug("Finished calculating mass at given depth.")

        if self.diff:
            logger.debug("Start calculating mass at given depth for file2.")
            mass2 = self.zone_crossing_event(self.obj2, self.lons, self.lats, h, bad1, bad2)
            logger.debug("Finished calculating mass at given depth for file2.")

            if self.abs:
                mass = np.copy(mass1) - mass2
            else:
                mass = (np.copy(mass1) - mass2) / np.copy(mass1)
            if self.clip:
                mass = self.clip_array(np.copy(mass))
            m, mid, M = self.get_colormap_midpoint(mass)
        
        logger.debug("Start plotting")
        ax.coastlines(zorder = 3, resolution='10m')

        if not self.diff:
            sm = ax.contourf(self.lons, self.lats, mass1.T, 20, transform=ccrs.PlateCarree(), cmap = cmap, zorder = 1,extend='both', extendfrac='auto')
            cb = plt.colorbar(sm, ax = ax, orientation="horizontal", shrink = self.shrink)
        else:
            cmap = self.shiftedColorMap(cmap, midpoint = mid, name='shifted')
            sm = ax.contourf(self.lons, self.lats, mass.T, 20, transform=ccrs.PlateCarree(), cmap = cmap, zorder = 1,extend='both', extendfrac='auto')
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
    
    def zone_crossing_event(self, obj:Open, lons, lats, depth, bad1, bad2):
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

        return self.mass_sum(lons, lats, depth, z, mass, lon, lat, status, self.seafloor_idx, bad1, bad2)

    # @jit(nopython = True)
    def mass_sum(self, lons, lats, depth, z, mass, lon, lat, status, sfidx, bad1, bad2):
        logger.debug("Start summing masses.")
        mass_at_depth = np.zeros((len(lons), len(lats)))
        for i in range(mass.shape[1]):

            #Remove bad particles.
            if i in bad1 or i in bad2:
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
            if not self.diff and not self.abs:
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
        bad = self.clean_dataset(self.obj)

        
        logger.debug("Start calculating distance.")
        distance, in_cell = self.calculate_horizontal_distance(self.obj, self.lons, self.lats, bad)
        
        D = np.copy(distance)
        rows, cols = np.where(D == 0)
        D[rows, cols] = np.nan
        rows, cols = np.where(in_cell == 0)
        D[rows, cols] = np.nan
        logger.debug("Start plotting")
        ax.coastlines(zorder = 3, resolution='10m')

        if self.clip:
            D = self.clip_array(np.copy(D))
        
        sm = ax.contourf(self.lons, self.lats, D.T, 20, transform=ccrs.PlateCarree(), cmap = cmap, zorder = 1, extend = "max")
        # ax.contourf(self.lons, self.lats, np.ma.masked_where(distance != 0, distance).T,colors='none', hatches=['xx'], extend='both')
        ax.scatter(self.lons[np.where(in_cell==0)[0]], self.lats[np.where(in_cell==0)[1]], color = "k", label = "moved out of cell")
        ax.legend()

        cb = plt.colorbar(sm, ax = ax, orientation="horizontal", shrink = self.shrink)
        cb.set_label(r"$|\vec{r}(t_F) - \vec{r}(t_0)|$ [km]")
        
        if self.title is not None:
            ax.set_title(self.title, fontweight = self.fontweight)

        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                  linewidth=2, color='k', alpha=0.8, linestyle='--')
        gl.xlocator = mticker.FixedLocator([-180, -90, 0, 90, 180])
        gl.ylocator = mticker.FixedLocator([-90, -45, 0, 45, 90])
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        ax.set_xticks([-180, -90, 0, 90, 180])
        ax.set_xticklabels([r"$-180^\circ$", r"$-90^\circ$", r"$0^\circ$", r"$90^\circ$", r"$180^\circ$"])
        ax.set_yticks([-90, -45, 0, 45, 90])
        ax.set_yticklabels([r"$-90^\circ$", r"$-45^\circ$", r"$0^\circ$", r"$45^\circ$", r"$90^\circ$"])

        plt.tight_layout()

        if self.outfile is not None:
            logger.debug("Saving output file.")
            plt.savefig(self.outfile, dpi = 300, bbox_inches = "tight")
        else:
            plt.show()
    
    def calculate_horizontal_distance(self, obj: Open, lons, lats, bad):

        from haversine import haversine, Unit
        from numba_progress import ProgressBar
        from tqdm import trange

        logger.debug("Reading required simulation properties.")

        lon = obj.get_property("lon")
        lon = np.ma.filled(lon, np.nan)

        lat = obj.get_property("lat")
        lat = np.ma.filled(lat, np.nan)
        status = obj.get_property("status")
        status = np.ma.MaskedArray(status.data, status.mask, float)
        status = np.ma.filled(status, np.nan)
        moving = obj.get_property("moving")
        moving = np.ma.MaskedArray(moving.data, moving.mask, float)
        moving = np.ma.filled(moving, np.nan)

        # @jit(nogil=True)
        def fast_distance_calculator(lons, lats, lon, lat, moving, bad):
            distance_traveled = np.zeros((len(lons), len(lats)))
            is_in_cell = np.ones((len(lons), len(lats)))

            for i in trange(lon.shape[1]):

                #Remove bad particles.
                if i in bad:
                    # progress.update(1)
                    continue

                drifter_lon = lon[:, i]
                drifter_lat = lat[:, i]

                movement_stopped = np.where(moving[:, i] == 0)[0]

                if len(movement_stopped) > 0:
                    movement_stopped = movement_stopped[0]
                else:
                    movement_stopped = len(drifter_lat) - 1
                
                pos1 = (drifter_lat[0], drifter_lon[0])
                pos2 = (drifter_lat[movement_stopped], drifter_lon[movement_stopped])

                lon_grid = np.argmin(np.abs(lons - pos1[1]))
                lat_grid = np.argmin(np.abs(lats - pos1[0]))
                if np.abs(pos2[0] - pos1[0]) > (lats[1] - lats[0]) / 2 or np.abs(pos1[1] - pos2[1]) > (lons[1] - lons[0]) / 2:
                    is_in_cell[lon_grid, lat_grid] = 0
                    d = haversine(pos1, pos2, unit = Unit.KILOMETERS)
                    distance_traveled[lon_grid, lat_grid] = d
                else:
                    # d = 2 * 6371 * np.arcsin(np.sqrt(np.sin((pos2[0] - pos1[0]) / 2) ** 2 + 
                    #                                  np.cos(pos2[0]) * np.cos(pos1[0]) * np.sin((pos2[1] - pos1[1]) / 2) ** 2))
                    
                    d = haversine(pos1, pos2, unit = Unit.KILOMETERS)
                    # print(d)
                    distance_traveled[lon_grid, lat_grid] = d
                # progress.update(1)
            return distance_traveled, is_in_cell
        
        # with ProgressBar(total = lon.shape[1]) as progress:
        #     d, in_cell = fast_distance_calculator(lons, lats, lon, lat, moving, bad, progress)

        d, in_cell = fast_distance_calculator(lons, lats, lon, lat, moving, bad)
        return d, in_cell

    def drifter_locations(self):
        from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
        import matplotlib.style as style

        style.use('seaborn-v0_8-bright')

        if self.loc is None:
            raise AttributeError("Locations are not specified.")
        
        plt.close()
        fig, ax = plt.subplots(1, 1, subplot_kw={'projection': ccrs.PlateCarree()}, figsize = self.figsize)

        logger.debug("Searching for bad trajectories.")
        bad = self.clean_dataset(self.obj)
        logger.debug("Reading required simulation properties.")

        lon = self.obj.get_property("lon")
        lon = np.ma.filled(lon, np.nan)

        lat = self.obj.get_property("lat")
        lat = np.ma.filled(lat, np.nan)

        xticks = []
        yticks = []
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
            
            ax.plot(lon[0, idx], lat[0, idx], "o", markersize = 10, label = "L" + str(k))
            if self.loclines:
                ax.hlines(lat[0, idx], -180, lon[0, idx], linestyles="dashed", color = "black")
                ax.vlines(lon[0, idx], -90, lat[0, idx], linestyle = "dashed", color = "black")
            xticks.append(lon[0, idx])
            yticks.append(lat[0, idx])
            k += 1
        
        ax.set_global()
        ax.add_feature(cartopy.feature.LAND, facecolor="gray",edgecolor='black', zorder = 1)
        ax.coastlines(zorder = 2)
        ax.set_xticks(xticks)
        ax.set_xticklabels([r"${}^\circ $".format(i) for i in xticks] , rotation = 30)
        ax.set_yticks(yticks)
        ax.set_yticklabels(r"${}^\circ $".format(i) for i in yticks)
        
        ax.legend()

        plt.tight_layout()

        if self.outfile is not None:
            logger.debug("Saving output file.")
            plt.savefig(self.outfile, dpi = 300, bbox_inches = "tight")
        else:
            plt.show()

    def drifter_properties(self):
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
        linestyles = ["solid", "dashed", "dotted", "dashdot"]
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
                    unitsy = "[" + self.obj.data[self.prop2].units + "]"
                else:
                    y = self.create_timedelta_array(lon.shape[0])
                    unitsy = "[h]"

                if self.prop1 == "time":
                    X = x
                else:
                    X = x[:, idx]
            
                if self.prop2 == "time":
                    Y = y
                else:
                    Y = y[:, idx]
                # color = next(self.ax._get_lines.prop_cycler)['color']
                line, = self.ax.plot(X, Y, lw = self.lw, linestyle = linestyles[k % 4])
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
        M = self.zone_crossing_event(self.obj, self.lons, self.lats, self.depth, self.clean_dataset(self.obj), [])
        if not self.abs:
            x = len(np.where(M > 0)[0])
        else:
            x = 1
        return np.sum(M[np.invert(np.isnan(M))]) / x
    
    def get_biome_weighted_mass_at_depth(self):
        logger.debug("Reading required simulation properties.")
        obj = self.obj
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
        return mass_by_biomes, np.sum(mass_by_biomes) * 10 ** (-15)