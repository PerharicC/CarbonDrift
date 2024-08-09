import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, CheckButtons
from matplotlib import animation
from matplotlib.gridspec import GridSpec

import cartopy
from cartopy import config
import cartopy.crs as ccrs

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
    ice = "Ice"
    seafloor = "Reached_Sea_Floor"
    ice_idx = meanings.index(ice)
    seafloor_idx = meanings.index(seafloor)
    return values[ice_idx], values[seafloor_idx]

class Plot:

    def __init__(self, file1, file2 = None, cmap = None,
                 lons = None, lats = None, figsize = (20, 20),
                 fontsize = 17, title = None, depth = -200,
                 diff = True, absolute = False, fontweight = "normal",
                 outfile = None, shrink = 1, clip = None):

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

        if file2 is not None:
            self.obj2 = Open(file2)
        
        logger.debug("Adding attributes.")

        if depth > 0:
            depth*=-1
        if depth <= -5000:
            logger.debug("Setting depth to sea_floor.")
            depth = "sea_floor"
        
        self.depth = depth
        self.cmap = cmap
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

        logger.debug("Decrypting status numeberings.")

        ice, seafloor = get_status_info(self.obj.data)
        self.ice_idx = ice
        self.seafloor_idx = seafloor

        logger.debug("Setting up clipping.")
        if clip is None:
            self.clip = False
        else:
            self.clip = True
            self.Vmin, self.Vmax = [float(i.replace("m", "-")) for i in clip.split(":")]

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
    
    def clip_mass(self, mass):
        m, M = np.min(mass[np.invert(np.isnan(mass))]), np.max(mass[np.invert(np.isnan(mass))])
        row, col = np.where(mass < self.Vmin)
        mass[row, col] = self.Vmin
        row, col = np.where(mass > self.Vmax)
        mass[row, col] = self.Vmax
        return mass

    def mass_map(self):
        """Plot the mass reached at depths [-200m, -1000m, sea_floor] over a cartopy map.

        Parameters
        -----------
        """

        plt.close()
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
        mass1, m0 = self.zone_crossing_event(self.obj, self.lons, self.lats, h, bad1, bad2)

        if self.clip and not self.diff:
            mass1 = self.clip_mass(np.copy(mass1))
        
        logger.debug("Finished calculating mass at given depth.")

        if self.diff:
            logger.debug("Start calculating mass at given depth for file2.")
            mass2, m02 = self.zone_crossing_event(self.obj2, self.lons, self.lats, h, bad1, bad2)
            logger.debug("Finished calculating mass at given depth for file2.")

            if self.abs:
                mass = np.copy(mass1) - mass2
            else:
                mass = (np.copy(mass1) - mass2) / np.copy(mass1)
            if self.clip:
                mass = self.clip_mass(np.copy(mass))
            m, mid, M = self.get_colormap_midpoint(mass)
        
        logger.debug("Start plotting")
        ax.coastlines(zorder = 3, resolution='10m')

        if not self.diff:
            sm = ax.contourf(self.lons, self.lats, mass1.T, 20, transform=ccrs.PlateCarree(), cmap = cmap, zorder = 1,extend='both', extendfrac='auto')
            cb = plt.colorbar(sm, ax = ax, orientation="horizontal", shrink = self.shrink)
            cb.set_label(r"$m/m_0$")
        else:
            cmap = self.shiftedColorMap(cmap, midpoint = mid, name='shifted')
            sm = ax.contourf(self.lons, self.lats, mass.T, 20, transform=ccrs.PlateCarree(), cmap = cmap, zorder = 1,extend='both', extendfrac='auto')
            cb = plt.colorbar(sm, ax = ax, orientation="horizontal", shrink = self.shrink)
            if self.abs:
                cb.set_label(r"$\Delta m$")
            else:
                cb.set_label(r"$\Delta m/m$")
        
        if self.title is not None:
            ax.set_title(self.title, fontweight = self.fontweight)

        plt.tight_layout()

        if self.outfile is not None:
            logger.debug("Saving output file.")
            plt.savefig(self.outfile, dpi = 300, bbox_inches = "tight")
        else:
            plt.show()
    
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

        # @jit(nopython = True)
        def mass_sum(lons, lats, depth, z, mass, lon, lat, sfidx, bad1 = bad1, bad2 = bad2):
            logger.debug("Start summing masses.")
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
                        sea_floor = np.where(s == sfidx)[0]
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
            logger.debug("Finish summing masses.")
            return mass_at_depth, m0
        
        return mass_sum(lons, lats, depth, z, mass, lon, lat, self.seafloor_idx)
    
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
            
            if (drifter_trajectory[trajectory_nans]>0).any() or np.any(drifter_mass[1:] >= 1):
                bad_trajectories.append(i)
            
            elif np.any(status[:, i] == self.ice_idx): #Ice
                bad_trajectories.append(i)
        return bad_trajectories