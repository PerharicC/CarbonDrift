import os

from model.areadecay.carbondrift import *
import numpy as np
from datetime import datetime, timedelta
from copy import copy
from netCDF4 import Dataset
import numpy.ma as ma

# import model.areadecay.carbondriftopen


class DivideAndConquerCarbonDrift:

    """A class that runs the CarbonDrift module in "economic" mode, by splitting the simulation into multiple subtasks, until
    either the ram and number of particles in subtask do not exceed the specified values, or the number of particles is < than a specified value.
    """
    
    def __init__(self, **kwargs):
        
        """
        Parameters
        ----------
        starttime: datetime object
            Beginning of simulation (default datetime.now())
        reader: list of readers
            The readers to be used in the simulation (defualt None)
        lon: 1D/2D ndarray or float
            Initial longtitude of seeded particles (default 0)
        lat: 1D/2D ndarray or float
            Initial latitude of seeded particles (default 0)
        configure: dict
            Dictionary of all properties, which should be configured.
        min_num: int
            Number of particles in a subtask, which cancel further splitting.
        deactivate_fragmentation: bool
            Boolean statement if fragmentation should be deactivated (default False)
        deactivate_horizontal_advection: bool
            Boolean statement if horizontal advection should be deactivated (default False)
        + All key arguments for CarbonDrift.
        """
        
        self.time = kwargs.pop("starttime", datetime.now())
        self.reader = kwargs.pop("reader", None)
        lon = kwargs.pop("lon", 0)
        lat = kwargs.pop("lat", 0)
        self.config = kwargs.pop("configure", {})
        self.min_particles_infile = kwargs.pop("min_num", 100)
        self.deactivate_frag = kwargs.pop("deactivate_fragmentation", False)
        self.deactivate_ha = kwargs.pop("deactivate_horizontal_advection", False)

        self.o_kwargs = kwargs

        self.o = CarbonDrift(**kwargs)

        if self.deactivate_frag:
            self.o.deactivate_fragmentation()
        
        if self.deactivate_ha:
            self.o.deactivate_horizontal_advection()
        
        self.o_kwargs["plot_distribution"] = False

        if self.reader is not None:
            self.o.add_reader(self.reader)

        self.configure()
        self.o.seed_elements(lon = lon, lat = lat, z = 0, time = self.time, mass = kwargs["mass"], number = len(kwargs["mass"]))

    def economic_run(self, **kwargs):

        """Run the simulation and split it into multiple instances, if the ram exceeds the defined value.
        
        Parameters
        ----------
        Same key arguments as for the method run in OceanDrift."""

        self.outfile_main = kwargs["outfile"].split(".")[0]
        self.o.run(**kwargs)
        self.r_kwargs = kwargs
        if type(self.r_kwargs["time_step"]) == timedelta:
            self.r_kwargs["time_step"] = self.r_kwargs["time_step"].total_seconds()
        self.steps = kwargs["steps"]

        if not self.o.resume:
            print("The RAM and number of allowed elements per simulation never exceeded the set limit.")
        else:
            self.tracker = self.outfile_main + "_tracker.txt"
            f = open(self.tracker, "w")
            f.write(kwargs["outfile"])
            f.close()

            self.recursive_run_divide(kwargs["outfile"])
            print("Merging files.")
            self.recursive_run_conquer()
            self.delete_in_between_files()
    
    def configure(self):
        for key, value in self.config.items():
            self.o.set_config(key, value)
        
    def remaining_steps(self, current_time):
        steps_done = (current_time - self.time).total_seconds() / self.r_kwargs["time_step"]
        return int(self.steps - steps_done + 1)
    
    def new_object_instance(self, file, k):

        """Reads a .nc file and creates a new CarbonDrift object with half the particles, 
        which are seeded at the last positions of the original .nc file.
        
        Parameters
        ----------
        file: str
            File name of .nc file to be read.
        k: int
            Which half of the particles should be taken from the original dataset (0, 1).
        
        Returns
        -------
        int: Number of reamining simulation steps."""

        #Open file and read the last time_step.
        o = carbondriftopen.open(file, self.o_kwargs["distribution"])
        elements_final = o.history[:, -1]
        mask = elements_final.mask

        elements_final = np.ma.filled(elements_final, np.nan)
        time_final = o.time
        num = len(o.history)

        dict = {}
        stop_split = False

        #Read and write the necessary variables to dictionary dict.
        if k == 0:
            for attribute in ["lat", "lon", "z", "mass", "age_seconds"]:
                dict[attribute] = elements_final[attribute][:num//2]
        else:
            for attribute in ["lat", "lon", "z", "mass", "age_seconds"]:
                dict[attribute] = elements_final[attribute][num//2:]
        
        dict["age_seconds"] -= o.time_step.total_seconds()
        
        # if len(dict["lat"]) <= self.min_particles_infile:
        #     stop_split = True
        
        del o
        
        kwargs = copy(self.o_kwargs)

        #If splitting should finish, because the number of particles is small enough.
        if stop_split:
            kwargs["max_ram"] = None
        
        #Create new object instance.
        self.o = CarbonDrift(**kwargs)

        if self.deactivate_frag:
            self.o.deactivate_fragmentation()
        
        if self.deactivate_ha:
            self.o.deactivate_horizontal_advection()
        
        if self.reader is not None:
            self.o.add_reader(self.reader)
        
        self.o.task_count = 1
        
        self.configure()
        self.o.seed_elements(lon = dict["lon"], lat = dict["lat"], z = dict["z"], mass = dict["mass"], age_seconds = dict["age_seconds"], time = time_final)
        return self.remaining_steps(time_final)-1
    
    def recursive_run_divide(self, file, l = 1, k = 0, split_number = 1, split_done = 0):

        """Recursive dividing of particles into two halfes and runing the new contracted dataset until the inital simulation time is reached.
        
        The devision of particles occurs, if the ram exceeds the max_ram key argument in the object initialisation and the number of particles is greater than 4.
        
        Parameters
        ----------
        file: str
            File name of previously ran dataset.
        l: int
            Vertical position in the division tree.
        k: int
            Horizontal position with respect to the parent in the division tree (0,1).
        split_number: int
            Current number of splits performed.
        split_done: int
            Curent number of finished splittings.
        
        Raises
        ------
        RecursionError: If l>10.
            """
        
        if l > 15:
            raise RecursionError("Stopping simulation, because the D&C tree has over 15 levels, i.e. the number of output files could exceed {0}.\n".format(2 ** l) +
                                 "Raise the max RAM value allowed or the min number of elements allowed in one file.")
        
        if (split_done == split_number + 1):
            print("Splitting tasks finished. The particles have been split into {} different tasks.".format(2 * split_number + 1))
            return 0
        
        kwargs = self.r_kwargs
        kwargs["steps"] = self.new_object_instance(file, k)
        file = file.split(".")[0] + str(k) + ".nc"
        
        kwargs["outfile"] = file
        
        # if self.reader is not None:
        #     self.o.add_reader(self.reader)
        
        self.o.run(**kwargs)

        f = open(self.tracker, "r+")
        current_files = f.read().split()
        if file not in current_files:
            f.write("\n" + file)
        f.close()

        if not self.o.resume:
            file = file.split(".")[0]
            if k  == 0:
                file = file[:-1] + ".nc"
                return self.recursive_run_divide(file, l, k + 1, split_number, split_done + 1)
            else:
                tree_raise = 1

                for num in file[::-1]:
                    if num == "1":
                        tree_raise += 1
                    else:
                        break

                file = file[:-tree_raise] + ".nc"
                return self.recursive_run_divide(file, l - 1, k, split_number, split_done + 1)
        
        return self.recursive_run_divide(file, l + 1, 0, split_number + 1, split_done)
    
    @staticmethod
    def merge_files(file1, file2):
        """Merge two child datasets d1 and d2 along the trajectory dimension into dataset d, 
        open their parent dataset d3, merge it with d along the time dimension and save as new dataset.
        
        Parameters
        ----------
        file1: str
            File name of first child dataset.
        file2: str
            File name of second child dataset.
            """
        
        d1 = Dataset(file1)
        d2 = Dataset(file2)
        file3 = file1.split(".")[0][:-1] + ".nc"
        d3 = Dataset(file3)
        trajectories_1 = len(d1["trajectory"])
        trajectories_2 = len(d2["trajectory"])
        trajectories_3 = len(d3["trajectory"])
        new_elements = trajectories_1 + trajectories_2 - trajectories_3
        #Split d1 and d2 to original particles and newly formed particles.
        original_num1 = trajectories_3 // 2
        original_num2 = trajectories_3 - original_num1

        variables12 = {}
        variables3 = {}

        d4 = Dataset("inbetween.nc", "w")

        for variable in d1.variables.keys():
            value = d1.variables[variable]
            value2 = d2[variable]

            if variable == "time":
                variables12[variable] = value[1:]
            
            elif variable == "trajectory":
                original1 = value[:original_num1]
                original2 = value2[:original_num2] + original1[-1]
                new1 = value[original_num1:]
                new2 = value2[original_num2:]
                original12 = np.append(original1, original2)
                if len(new1) > 0 and len(new2) > 0:
                    new1 = np.arange(original12[-1] + 1, len(new1) + 1 + original12[-1], 1)
                    new2 = np.arange(new1[-1] + 1, len(new2) + 1 + new1[-1], 1)
                    new12 = np.append(new1, new2)
                elif len(new1) > 0 or len(new2) > 0:
                    new12 = np.arange(original12[-1] + 1, max(len(new1), len(new2)) + 1 + original12[-1], 1)
                else:
                    new12 = []
                variables12[variable] = np.append(original12, new12)
            
            else:
                original1 = value[:original_num1, 1:]
                original2 = value2[:original_num2, 1:]
                new1 = value[original_num1:, 1:]
                new2 = value2[original_num2:, 1:]
                original12 = ma.append(original1, original2, axis = 0)
                variables12[variable] = ma.append(original12, ma.append(new1, new2, axis = 0), axis = 0)
            
            value = d3.variables[variable]

            if variable == "time":
                variables3[variable] = value
            elif variable == "trajectory":
                variables3[variable] = np.append(value, np.arange(value[-1] + 1, value[-1] + new_elements + 1, 1))
            else:
                mask = True
                if variable== "mass":
                    mask = False
                variables3[variable] = ma.append(value, ma.MaskedArray(np.zeros((new_elements, len(d3["time"]))), mask = mask), axis = 0)
            
            if variable == "time":
                dim = d3.dimensions[variable]
                size = len(dim) + len(variables12[variable])
                d4.createDimension(variable, size if not dim.isunlimited() else None)
            elif variable == "trajectory":
                dim = d3.dimensions[variable]
                size = trajectories_1 + trajectories_2
                d4.createDimension(variable, size if not dim.isunlimited() else None)
            
            if variable in ["time", "trajectory"]:
                new_variable = d4.createVariable(variable, value.dtype, value.dimensions)
                if variable == "time":
                    new_values = ma.append(variables3[variable], variables12[variable])
                else:
                    new_values = variables3[variable]
                new_variable.setncatts(value.__dict__)
                new_variable[:] = new_values
            else:
                new_variable = d4.createVariable(variable, value.dtype, value.dimensions)
                new_values = ma.append(variables3[variable], variables12[variable], axis = 1)
                new_variable.setncatts(value.__dict__)
                new_variable[:] = new_values
            
            variables12 = {}
            variables3 = {}
        
        d1.close()
        d2.close()
        d3.close()
        d4.close()
        os.remove(file3)
        os.rename("inbetween.nc", file3)
    
    def recursive_run_conquer(self):

        """Loop through all files, while merging them like a puzzle, starting at the lowest point in the divided tree and moving up step by step."""

        #Import all the created file names.
        f = open(self.tracker, "r")
        file_names = f.read().split()
        f.close()

        #Find the lowest point in the divided tree.
        l = len(max(file_names, key = len))

        y = []
        try:
            while len(file_names) > 1:

                for file in file_names:
                    if len(file) == l:
                        name = file.split(".")[0]
                        file2 = name[:-1] + str(int(not(int(name[-1])))) + ".nc"

                        if file2 in y:
                            if name[-1] == "0":
                                self.merge_files(file, file2)
                            else:
                                self.merge_files(file2, file)
                        
                        y.append(file)
                    
                for file in y:
                    file_names.remove(file)
                y = []
                l -= 1
            print("Merging files finished.")
        except Exception as e:
            logging.error(traceback.format_exc())
        
    def delete_in_between_files(self):

        """Delete all the created in between files, except the inital parent file."""

        f = open(self.tracker, "r")
        file_names = f.read().split()
        f.close()

        for file in file_names[1:]:
            os.remove(file)
        
        os.remove(self.tracker)
        print("All in between files deleted.")

