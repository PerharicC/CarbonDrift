import numpy as np
from numpy.random import random
import scipy
from opendrift.models.oceandrift import OceanDrift, Lagrangian3DArray
from datetime import datetime, timedelta
# import logging
# logging.captureWarnings(True)
# logger = logging.getLogger(__name__)
from model.logger import Logger

log = Logger("CarbonDrift.model.areadecay.carbondrift")
logger = log.LOGGER

import opendrift
import geojson

from opendrift.readers import reader_global_landmask
import sys
import traceback

from model.areadecay.decaydistributions import Decay_Functions

import psutil
from numba import jit



class Carbon(Lagrangian3DArray):

    variables = Lagrangian3DArray.add_variables([
    ('mass', {'dtype': np.float32,
                      'units': 'g',
                      'default': 1})])



class CarbonDrift(OceanDrift):

    """An OceanDrift inherited class for tracking carbon particles, as they decay through area microbial decay and fragmentation.
    
    The area microbial decay is goverened by the two coupled equations:

        dm/dt = -k * S,

        w/w_0 = (m / m_0)^(1/6).
    
    The fragmentation is implemented as a random decay process, corresponding to some user specified distribution.

    """

    ElementType = Carbon

    required_variables = {
        "x_sea_water_velocity": {'fallback': 0},
        'y_sea_water_velocity': {'fallback': 0},
        'land_binary_mask': {'fallback': None},
        "sea_water_temperature": {"fallback": None},
        "depth": {"fallback": None},
    }

    def __init__(self, *args, **kwargs):

        """
        Parameters
        ----------
        mass_threshold: float
            The multiplication constant in the mass distribution. Should be between (0, 1] (default 1)
        z_threshold: float
            The multiplication constant in the z distribution. Should be between (0, 1] (default 1)
        threshold: float
            The multiplication constant in the combined mass/z distribution. Should be between (0, 1] (default 1)
        decay_stop: float
            The fraction of the initial mass at which the particles deactivate (default 0.00)
        initial_velocity: float
            The initial velocity of the particles (default -1)
        mass: float, ndarray
            The initial masses of all particles in grams (default 1 g)
        depth_norm: float
            Normalizing depth factor for z_linear distribution (default 1000)
        distribution: str
            Name of the method used as the decay distribution
        plot_distribution: bool
            Plot the distribution (default False)
        max_ram: float
            The ram threshold in GiB, which when exceeded splits the simulation into two subtasks 
            (only possible using the DivideAndConquerCarbonDrift module) (default None)
        max_num: int
            The max number of particles a simulation can contain, before it splits into two subtasks
            (only possible using the DivideAndConquerCarbonDrift module) (default None)
        decay_type: str
            Type of decay linear/exponential (default "linear")
        +All other arguments that OceanDrift can accept.
        """
        
        self.mass_threshold = kwargs.pop("mass_threshold", 1)
        self.z_threshold = kwargs.pop("z_threshold", 1)
        self.threshold = kwargs.pop("threshold", 1)
        self.decay_stop = kwargs.pop("decay_stop", 0.)
        self.w0 = kwargs.pop("initial_velocity", -0.01)
        m0 = kwargs.pop("mass", None)
        self.r0 = kwargs.pop("r0", 0.1)
        depth_norm = kwargs.pop("depth_norm", 1000)
        distribution = kwargs.pop("distribution", "mass_sqrt")
        plot_dist = kwargs.pop("plot_distribution", False)
        self.max_ram = kwargs.pop("max_ram", None)
        self.max_num = kwargs.pop("max_num", None)
        self.decay_type = kwargs.pop("decay_type", "linear")
        
        super(CarbonDrift, self).__init__(*args, **kwargs)
        
        if m0 is not None:
            self.m0 = np.max(m0)
        else:
            self.m0 = self.elements.variables["mass"]["default"]
        

        self.mass_idx = list(vars(self.elements).keys()).index("mass") + 2 #Probably dont need this.

        #Initialize fragmentation decay distribution.
        dist = Decay_Functions(self.m0, depth_norm, self.mass_threshold, self.z_threshold, self.threshold)
        self.func = dist.choose_distribution(distribution, plot_dist)(dist)

        #Initialize DivideAndConquer algorithm splitting attributes.
        self.pause = False
        self.resume = True

        #Allow fragmentation.
        self.no_decay = False

        #Dont allow deactivation.
        self.deactivate = False

        #Allow horizontal advection.
        self.horizontal_advection = False
    
    def update(self):
        self.deactivate_elements(self.environment.sea_water_temperature < 0, reason ="Ice")
        if self.horizontal_advection and (self.steps_calculation > 1 or self.task_count == 0):
            self.advect_ocean_current()

        self.elements.mass = self.microbial_decay()
        
        dz = self.time_step.total_seconds() * self.vertical_velocity()

        self.elements.z = self.elements.z + dz

        #Check if particles reached sea floor.
        depths = self.environment.depth
        bottom = self.elements.z <= depths
        self.elements.z[bottom] = depths[bottom]

        self.deactivate_elements(bottom, reason = "Reached_Sea_Floor")
        

        if not self.no_decay and (self.expected_steps_calculation - 1) * self.time_step + self.start_time != self.time:
            self.create_new_particles(*self.fragmentation())
        
        #For DivideAndConquer.
        if self.max_ram is not None:
            if self.ram_check():
                self.pause = True
        
        #For DivideAndConquer.
        if self.max_num is not None:
            if len(self.elements.mass) > self.max_num:
                self.pause = True
    
    def deactivate_fragmentation(self):
        self.no_decay = True
    
    def enable_deactivation(self):
        self.deactivate = True
    
    def deactivate_horizontal_advection(self):
        self.horizontal_advection = False
    
    def activate_horizontal_advection(self):
        self.horizontal_advection = True

    def microbial_decay(self):
        constants = 4 * np.pi * self.time_step.total_seconds() * self.r0 ** 2 / self.m0 ** (2 / 3)
        return self.elements.mass ** (2/3) * (self.elements.mass ** (1 / 3) - constants * self.decay_coef()) 
    
    def decay_coef(self):

        """Compute the microbial coeficient from current sea water temeperatures."""

        if self.decay_type == "linear":
            return (self.environment.sea_water_temperature * 0.064 + 0.02) / (24 * 3600) #Return linear decay.
        return 0.140 * np.exp(0.145 * self.environment.sea_water_temperature) / (24 * 3600) #Return exponential decay.
    
    def vertical_velocity(self):
        return self.w0 * (self.elements.mass / self.m0) ** (1 / 6)

    def decay_check(self):
        
        """Check which particles will decay in this time step.

        Returns
        -------
        ndarray:
            Boolean ndarray of size equal to current particle number, where True means that the particle will decay.
        """

        mass = self.elements.mass
        M = np.copy(mass)
        size_check = M > self.decay_stop * self.m0
        #Insure randomness.
        np.random.seed()

        return (self.threshold * self.func(mass, self.elements.z) > random(size = len(mass))) & size_check
    
    def fragmentation(self):
        mass = self.elements.mass
        M = np.copy(mass)
        decay = self.decay_check()
        decay_particles_mass = M[decay]
        return decay, decay_particles_mass, decay_particles_mass * random(size = len(decay_particles_mass))
    
    def create_new_particles(self, decay, m, new_m):

        """Creates two new particles from the original particle. The mass is conserved."""

        #Calculate second particle mass and change the original particle mass with the first particle mass.
        new_m2 = m - new_m
        self.elements.mass[decay] = new_m
        decay_IDS = self.elements.ID[decay]

        #Create new IDs.
        current_particle_number = self.history.shape[0]
        ID = np.arange(current_particle_number + 1, current_particle_number + 1 + len(m), 1)
        
        #Add the second particle to the elements attribute, inheriting all
        #original particle properties, except for the particle's ID and mass.
        for attribute, value in vars(self.elements).items():

            if attribute == "ID":
                self.elements.ID = np.append(self.elements.ID, ID)
                
            elif attribute == "mass":
                self.elements.mass = np.append(self.elements.mass, new_m2)

            elif attribute == "dtype":
                continue

            else:
                new_vals = value[decay]
                setattr(self.elements, attribute, np.append(value, new_vals))
        
        if len(np.where(decay == True)[0]) > 0:
            orig_history = np.ma.MaskedArray.copy(self.history)
            new_history, new_mask = self.update_history(orig_history.data, orig_history.mask, decay, ID, decay_IDS)
            new_history = np.ma.array(new_history, mask = new_mask)
            self.history = np.ma.vstack([self.history, new_history])
    
    @staticmethod
    @jit(nopython = True, parallel = True)
    def update_history(orig_history_data, orig_history_mask, decay, ID, decay_IDS):

        """Updates the history array with newly formed particles, so that the shape of the array is not changed through time.

        The mass of new particles is set to 0 for all times before the decay."""
        new_history = []
        new_mask = []
        for k, i in enumerate(np.where(decay == True)[0]):

            history = np.copy(orig_history_data[decay_IDS[k] - 1])
            history_mask = np.copy(orig_history_mask[decay_IDS[k] - 1])
            j = 0
            for time_step in history:
                time_step[0] = ID[k]

                time_step["mass"] = 0

                history[j] = time_step
                j+=1
            new_history.append(history)
            new_mask.append(history_mask)

        return new_history, new_mask
    
    #The same as the original run() function in basemodel, with a slight order change in the outputfile initialisation and export buffer size.
    #The pause and resume attributes are also changed according to the simulation output.

    def run(self,
            time_step=None,
            steps=None,
            time_step_output=None,
            duration=None,
            end_time=None,
            outfile=None,
            export_variables=None,
            export_buffer_length=100,
            stop_on_error=False):
        """Start a trajectory simulation, after initial configuration.

        Performs the main loop:
            - Obtain environment data for positions of all particles.
            - Call method 'update' to update (incl advect) particle properties.
        until one of the following conditions are met:
            - Maximum number of steps are reached
            - A needed variable can not be obtained by any reader
                (outside spatial/temporal domain) and has no fallback
                (default) value.
            - All particles have been deactivated (e.g. by stranding)
            - Occurance of any error, whose trace will be output to terminal.

        Before starting a model run, readers must be added for all
        required variables, unless fallback values have been specified.
        Some particles/elements must have been scheduled for seeding, and the
        run will start at the time when the first element has been scheduled..

        Arguments:
            time_step: interval between particles updates, in seconds or as
                timedelta. Default: 3600 seconds (1 hour)
            time_step_output: Time step at which element properties are stored
                and eventually written to file.
                Timedelta object or seconds.
                Default: same as time_step, meaning that all steps are stored
            The length of the simulation is specified by defining one
                (and only one) of the following parameters:
                - steps: integer, maximum number of steps. End of simulation
                will be self.start_time + steps*self.time_step
                - duration: timedelta defining the length of the simulation
                - end_time: datetime object defining the end of the simulation
            export_variables: list of variables and parameter names to be
                saved to file. Default is None (all variables are saved)
        """

        # Exporting software and hardware specification, for possible debugging
        logger.debug(opendrift.versions())

        self.timer_end('configuration')
        self.timer_start('preparing main loop')

        if self.num_elements_scheduled() == 0:
            raise ValueError('Please seed elements before starting a run.')
        self.elements = self.ElementType()

        # Export seed_geojson as FeatureCollection string
        self.add_metadata('seed_geojson',
                          geojson.FeatureCollection(self.seed_geojson))

        if outfile is None and export_buffer_length is not None:
            logger.debug('No output file is specified, '
                         'neglecting export_buffer_length')
            export_buffer_length = None

        # Some cleanup needed if starting from imported state
        if self.steps_calculation >= 1:
            self.steps_calculation = 0
        if self.history is not None:
            # Delete history matrix before new run
            self.history = None
            # Renumbering elements from 0 to num_elements, necessary fix when
            # importing from file, where elements may have been deactivated
            # TODO: should start from 1?
            self.elements.ID = np.arange(0, self.num_elements_active())

        ########################
        # Simulation time step
        ########################
        if time_step is None:
            time_step = timedelta(
                minutes=self.get_config('general:time_step_minutes'))
        if type(time_step) is not timedelta:
            # Time step may be given in seconds, as alternative to timedelta
            time_step = timedelta(seconds=time_step)
        self.time_step = time_step
        if time_step_output is None:
            time_step_output = self.get_config(
                'general:time_step_output_minutes')
            if time_step_output is None:
                self.time_step_output = self.time_step
            else:
                self.time_step_output = timedelta(minutes=time_step_output)
        else:
            if type(time_step_output) is timedelta:
                self.time_step_output = time_step_output
            else:
                self.time_step_output = timedelta(seconds=time_step_output)
            if self.time_step_output.days >= 0 and self.time_step.days < 0:
                self.time_step_output = -self.time_step_output

        time_step_ratio = self.time_step_output.total_seconds() / \
            self.time_step.total_seconds()
        if time_step_ratio < 1:
            raise ValueError('Output time step must be equal or larger '
                             'than calculation time step.')
        if not time_step_ratio.is_integer():
            raise ValueError('Ratio of calculation and output time steps '
                             'must be an integer - given ratio is %s' %
                             time_step_ratio)
        ########################
        # Simulation duration
        ########################
        if time_step.days < 0:
            logger.info(
                'Backwards simulation, starting from last seeded element')
            self.start_time = self.elements_scheduled_time.max()
        if (duration is not None and end_time is not None) or \
            (duration is not None and steps is not None) or \
                (steps is not None and end_time is not None):
            raise ValueError('Only one of "steps", "duration" and "end_time" '
                             'may be provided simultaneously')
        if duration is None and end_time is None:
            if steps is not None:
                duration = steps * self.time_step
            else:
                for reader in self.readers.values():
                    if reader.end_time is not None:
                        if end_time is None:
                            end_time = reader.end_time
                        else:
                            end_time = min(end_time, reader.end_time)
                    logger.info('Duration, steps or end time not specified, '
                                'running until end of first reader: %s' %
                                (end_time))
        if duration is None:
            duration = end_time - self.start_time

        if time_step.days < 0 and duration.days >= 0:
            # Duration shall also be negative for backwards run
            duration = -duration

        if np.sign(duration.total_seconds()) * np.sign(
                time_step.total_seconds()) < 0:
            raise ValueError(
                "Time step must be negative if duration is negative.")

        self.expected_steps_output = duration.total_seconds() / \
            self.time_step_output.total_seconds() + 1  # Includes start and end
        self.expected_steps_calculation = duration.total_seconds() / \
            self.time_step.total_seconds()
        self.expected_steps_output = int(self.expected_steps_output)
        self.expected_steps_calculation = int(self.expected_steps_calculation)
        self.expected_end_time = self.start_time + self.expected_steps_calculation * self.time_step

        ##############################################################
        # Prepare readers for the requested simulation domain/time
        ##############################################################
        max_distance = \
            self.get_config('drift:max_speed')*self.expected_steps_calculation * \
            np.abs(self.time_step.total_seconds())
        deltalat = max_distance / 111000.
        deltalon = deltalat / np.cos(
            np.radians(np.mean(self.elements_scheduled.lat)))
        # TODO: extent should ideally be a general polygon, not only lon/lat-min/max
        # TODO: Should also take into account eventual lifetime of elements
        simulation_extent = [
            np.maximum(-360,
                       self.elements_scheduled.lon.min() - deltalon),
            np.maximum(-89,
                       self.elements_scheduled.lat.min() - deltalat),
            np.minimum(360,
                       self.elements_scheduled.lon.max() + deltalon),
            np.minimum(89,
                       self.elements_scheduled.lat.max() + deltalat)
        ]
        if simulation_extent[2] == 360 and simulation_extent[0] < 0:
            simulation_extent[0] = 0
        logger.debug(
            'Finalizing environment and preparing readers for simulation coverage (%s) and time (%s to %s)'
            % (simulation_extent, self.start_time, self.expected_end_time))
        
        # Store expected simulation extent, to check if new readers have coverage
        self.simulation_extent = simulation_extent
        self.env.finalize(self.simulation_extent)
        
        ####################################################################
        # Preparing history array for storage in memory and eventually file
        ####################################################################
        if export_buffer_length is None:
            self.export_buffer_length = self.expected_steps_output
        else:
            self.export_buffer_length = export_buffer_length
        if steps > 98:
            self.export_buffer_length = steps + 2
        if self.time_step.days < 0:
            # For backwards simulation, we start at last seeded element
            logger.info('Backwards simulation, starting at '
                        'time of last seeded element')
            self.time = self.elements_scheduled_time.max()
            # Flipping ID array, so that lowest IDs are released first
            self.elements_scheduled.ID = \
                np.flipud(self.elements_scheduled.ID)
        else:
            # Forward simulation, start time has been set when seeding
            self.time = self.start_time

        # Add the output variables which are always required
        if export_variables is not None:
            export_variables = list(
                set(export_variables + ['lon', 'lat', 'ID', 'status']))
        self.export_variables = export_variables
        # Initialise array to hold history (element properties and environment)
        # for export to file.
        history_dtype_fields = [(name,
                                 self.ElementType.variables[name]['dtype'])
                                for name in self.ElementType.variables]
        # Add environment variables
        self.history_metadata = self.ElementType.variables.copy()
        for env_var in self.required_variables:
            history_dtype_fields.append((env_var, np.dtype('float32')))
            self.history_metadata[env_var] = {}

        # Remove variables from output array, if only subset is requested
        if self.export_variables is not None:
            history_dtype_fields = [
                f for f in history_dtype_fields
                if f[0] in self.export_variables
            ]
            for m in list(self.history_metadata):
                if m not in self.export_variables:
                    del self.history_metadata[m]

        history_dtype = np.dtype(history_dtype_fields)
        self.history = np.ma.array(np.zeros(
            (len(self.elements_scheduled), self.export_buffer_length)),
                                   dtype=history_dtype)
        self.history.mask = True
        self.steps_exported = 0

        if outfile is not None:
            self.outfile = True
        else:
            self.outfile = None

        # Move point seeded on land to ocean
        if self.get_config('seed:ocean_only') is True and \
            ('land_binary_mask' in self.required_variables):
            #('land_binary_mask' not in self.fallback_values) and \
            self.timer_start('preparing main loop:moving elements to ocean')
            self.elements_scheduled.lon[np.invert(np.isnan(self.elements_scheduled.lon))], self.elements_scheduled.lat[np.invert(np.isnan(self.elements_scheduled.lat))] = \
                self.closest_ocean_points(self.elements_scheduled.lon[np.invert(np.isnan(self.elements_scheduled.lon))],
                                          self.elements_scheduled.lat[np.invert(np.isnan(self.elements_scheduled.lat))])
            self.timer_end('preparing main loop:moving elements to ocean')

        #############################
        # Check validity domain
        #############################
        validity_domain = [
            self.get_config('drift:deactivate_west_of'),
            self.get_config('drift:deactivate_east_of'),
            self.get_config('drift:deactivate_south_of'),
            self.get_config('drift:deactivate_north_of')
        ]
        if validity_domain == [None, None, None, None]:
            self.validity_domain = None
        else:
            self.validity_domain = validity_domain

        #############################
        # Model specific preparation
        #############################
        self.prepare_run()

        ##########################
        # Main loop
        ##########################
        self.add_metadata('simulation_time', datetime.now())
        self.timer_end('preparing main loop')
        self.timer_start('main loop')
        for i in range(self.expected_steps_calculation):
            try:
                # Release elements
                self.release_elements()

                if self.num_elements_active(
                ) == 0 and self.num_elements_scheduled() > 0:
                    self.steps_calculation += 1
                    logger.info(
                        'No active but %s scheduled elements, skipping timestep %s (%s)'
                        % (self.num_elements_scheduled(),
                           self.steps_calculation, self.time))
                    self.state_to_buffer()  # Append status to history array
                    if self.time is not None:
                        self.time = self.time + self.time_step
                    continue

                self.increase_age_and_retire()

                self.interact_with_seafloor()

                if self.show_continuous_performance is True:
                    logger.info(self.performance())
                # Display time to terminal
                logger.debug('===================================' * 2)
                logger.info('%s - step %i of %i - %i active elements '
                            '(%i deactivated)' %
                            (self.time, self.steps_calculation + 1,
                             self.expected_steps_calculation,
                             self.num_elements_active(),
                             self.num_elements_deactivated()))
                logger.debug('%s elements scheduled.' %
                             self.num_elements_scheduled())
                logger.debug('===================================' * 2)

                if len(self.elements.lon) > 0:
                    lonmin = self.elements.lon.min()
                    lonmax = self.elements.lon.max()
                    latmin = self.elements.lat.min()
                    latmax = self.elements.lat.max()
                    zmin = self.elements.z.min()
                    zmax = self.elements.z.max()
                    if latmin == latmax:
                        logger.debug('\t\tlatitude =  %s' % (latmin))
                    else:
                        logger.debug('\t\t%s <- latitude  -> %s' %
                                     (latmin, latmax))
                    if lonmin == lonmax:
                        logger.debug('\t\tlongitude = %s' % (lonmin))
                    else:
                        logger.debug('\t\t%s <- longitude -> %s' %
                                     (lonmin, lonmax))
                    if zmin == zmax:
                        logger.debug('\t\tz = %s' % (zmin))
                    else:
                        logger.debug('\t\t%s   <- z ->   %s' % (zmin, zmax))
                    logger.debug('---------------------------------')

                self.environment, self.environment_profiles, missing = \
                    self.env.get_environment(list(self.required_variables),
                                         self.time,
                                         self.elements.lon,
                                         self.elements.lat,
                                         self.elements.z,
                                         self.required_profiles,
                                         self.profiles_depth)

                self.store_previous_variables()

                if self.pause:
                    break

                self.calculate_missing_environment_variables()

                if any(missing):
                    self.report_missing_variables()

                self.interact_with_coastline()

                self.interact_with_seafloor()

                self.deactivate_elements(missing, reason='missing_data')

                self.state_to_buffer()  # Append status to history array

                self.remove_deactivated_elements()

                # Propagate one timestep forwards
                self.steps_calculation += 1

                if self.num_elements_active(
                ) == 0 and self.num_elements_scheduled() == 0:
                    raise ValueError(
                        'No more active or scheduled elements, quitting.')

                # Store location, in case elements shall be moved back
                self.store_present_positions()

                #####################################################
                if self.num_elements_active() > 0:
                    logger.debug('Calling %s.update()' % type(self).__name__)
                    self.timer_start('main loop:updating elements')
                    self.update()
                    self.timer_end('main loop:updating elements')
                else:
                    logger.info('No active elements, skipping update() method')
                #####################################################

                self.horizontal_diffusion()

                if self.num_elements_active(
                ) == 0 and self.num_elements_scheduled() == 0:
                    raise ValueError(
                        'No active or scheduled elements, quitting simulation')

                logger.debug('%s active elements (%s deactivated)' %
                             (self.num_elements_active(),
                              self.num_elements_deactivated()))
                # Updating time
                if self.time is not None:
                    self.time = self.time + self.time_step
                
            except Exception as e:
                message = ('The simulation stopped before requested '
                           'end time was reached.')
                logger.warning(message)
                self.store_message(message)
                logger.info('========================')
                logger.info('End of simulation:')
                logger.info(e)
                logger.info(traceback.format_exc())
                logger.info(self.get_messages())
                if not hasattr(self, 'environment'):
                    sys.exit('Simulation aborted. ' + self.get_messages())
                logger.info('========================')
                if stop_on_error is True:
                    sys.exit('Stopping on error. ' + self.get_messages())
                if self.steps_calculation <= 1:
                    raise ValueError('Simulation stopped within '
                                     'first timestep. ' + self.get_messages())
                break

        self.timer_end('main loop')
        self.timer_start('cleaning up')
        logger.debug('Cleaning up')

        if outfile is not None:
            self.io_init(outfile)

        self.interact_with_coastline(final=True)

        self.state_to_buffer()  # Append final status to buffer

        #############################
        # Add some metadata
        #############################
        for var in self.required_variables:
            keyword = 'reader_' + var
            if var not in self.env.priority_list:
                fallback = self.get_config(f'environment:fallback:{var}')
                if fallback is not None:
                    self.add_metadata(keyword, fallback)
                else:
                    self.add_metadata(keyword, None)
            else:
                readers = self.env.priority_list[var]
                if readers[0].startswith(
                        'constant_reader') and var in self.env.readers[
                            readers[0]]._parameter_value_map:
                    self.add_metadata(
                        keyword, self.env.readers[
                            readers[0]]._parameter_value_map[var][0])
                else:
                    self.add_metadata(keyword, self.env.priority_list[var])

        self.timer_end('cleaning up')
        self.timer_end('total time')
        if outfile is not None:
            logger.debug('Writing and closing output file: %s' % outfile)
            # Write buffer to outfile, and close
            if self.steps_output >= self.steps_exported:
                # Write last lines, if needed
                self.io_write_buffer()
            self.io_close()

        # Remove any elements scheduled for deactivation during last step
        self.remove_deactivated_elements()

        if export_buffer_length is None:
            # Remove columns for unseeded elements in history array
            if self.num_elements_scheduled() > 0:
                logger.info(
                    'Removing %i unseeded elements from history array' %
                    self.num_elements_scheduled())
                mask = np.ones(self.history.shape[0], dtype=bool)
                mask[self.elements_scheduled.ID - 1] = False
                self.history = self.history[mask, :]

            # Remove rows for unreached timsteps in history array
            self.history = self.history[:, range(self.steps_output)]
        else:  # If output has been flushed to file during run, we
            # need to reimport from file to get all data in memory
            del self.environment
            if hasattr(self, 'environment_profiles'):
                del self.environment_profiles
            self.io_import_file(outfile)

        if not self.pause:
            self.resume = False
    
    def ram_check(self):
        return psutil.Process().memory_info().rss / (10 ** 9) >= self.max_ram


