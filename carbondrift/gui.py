import numpy as np
import matplotlib.pyplot as plt
from PIL import ImageTk, Image
import tkinter
from tkinter import ttk
import threading
import time
import os
# from importlib import resources

from datetime import datetime, timedelta
# from copy import copy
from opendrift.readers import reader_global_landmask
from opendrift.readers import reader_netCDF_CF_generic
import carbondrift
# import carbondrift.models.massdecay.gridrun as mcdgrid
# import carbondrift.models.areadecay.gridrun as acdgrid
import carbondrift.models.massdecay.carbondrift as mcd
import carbondrift.models.areadecay.carbondrift as acd
# from carbondrift.simulation.param_classifier import Parameters
# from opendrift.readers import reader_global_landmask
# from opendrift.readers import reader_netCDF_CF_generic
from carbondrift.simulation.seeding import *

from carbondrift.models.logger import Logger

log = Logger("carbondrift.carbondrift_gui")
logger = log.LOGGER

class CarbonDriftGUI(tkinter.Tk):
    def __init__(self):
        tkinter.Tk.__init__(self)
        self.title("CarbonDrift " + carbondrift.__version__ + " GUI")

        self.supplementary_data_dir = os.path.join(
            os.path.dirname(os.path.dirname(carbondrift.__file__)),
            "supplementary_data")

        self.n = ttk.Notebook(self.master)
        self.n.grid()
        self.prepare = ttk.Notebook(self.n)
        self.run = ttk.Frame(self.n)
        self.plot = ttk.Frame(self.n)
        self.simseed = ttk.Frame(self.prepare)
        self.simreader = ttk.Frame(self.prepare)
        self.simconfigure = ttk.Frame(self.prepare)
        self.simtype = ttk.Frame(self.prepare)

        self.n.add(self.prepare, text='Prepare Simulation')
        self.n.add(self.run, text='Run Simulation')
        self.n.add(self.plot, text='Plot Simulation')
        self.prepare.add(self.simseed, text='Seed')
        self.prepare.add(self.simreader, text="Readers")
        self.prepare.add(self.simconfigure, text="Configure")
        self.prepare.add(self.simtype, text ="Simulation Type")

        self.frame_subtitle_font = ("Courier", 14, "bold")
        self.frame_font = None#("Helvetica", 12)

        self.setup_seed_frame()
        self.setup_reader_frame()
        self.setup_configure_frame()
        self.setup_simulation_type_frame()

        self.thread = None
        self.sdata = None

    def setup_seed_frame(self):
        self.seed_from_file = tkinter.Frame(self.simseed, bg='lightgray', bd=2,
                            relief=tkinter.SUNKEN, pady=25, padx=25)
        self.seed_from_file.grid(row=0, column=0, rowspan=1, padx=20, pady=(20,0), sticky="we")
        
        seed_from_file_title = tkinter.Label(self.seed_from_file, text="Seed from pre-made file", 
                                             font=self.frame_subtitle_font).grid(row=0, pady=(0, 25), columnspan=2)
        
        seed_file_label = tkinter.Label(self.seed_from_file, text="File Path:", font=self.frame_font).grid(row=1, column=0, sticky="w")
        self.seed_file_entry = tkinter.Entry(self.seed_from_file)
        self.seed_file_entry.bind("<Key>", self.disable_new_seed)
        self.seed_file_entry.grid(row=1, column=1, sticky="w")


        self.seed_new = tkinter.Frame(self.simseed, bg='lightgray', bd=2,
                            relief=tkinter.SUNKEN, pady=25, padx=25)
        self.seed_new.grid(row=1, column=0, rowspan=1, padx=20, pady=(20,0), sticky="we")
        seed_from_file_title = tkinter.Label(self.seed_new, text="Make New Seed", 
                                             font=self.frame_subtitle_font).grid(row=0, pady=(0, 25), columnspan=2, sticky="w")
        
        self.RELATIVEMASS = tkinter.BooleanVar()
        relative_mass_button = tkinter.Checkbutton(self.seed_new, text="relative mass", command=self.on_check_relative_mass,
                                                            variable=self.RELATIVEMASS).grid(row=1, column=0, sticky="w")
        
        self.SAVESEED = tkinter.BooleanVar()
        save_seed_button = tkinter.Checkbutton(self.seed_new, text="save seed", command=self.on_check_save_seed,
                                                            variable=self.SAVESEED).grid(row=1, column=1, sticky="w")
        lonmin_label = tkinter.Label(self.seed_new, text="lonmin", font=self.frame_font).grid(row=2, column=0)
        self.lonmin_entry = tkinter.Entry(self.seed_new)
        self.lonmin_entry.insert(0, "-180")
        self.lonmin_entry.grid(row=3, column=0)

        lonmax_label = tkinter.Label(self.seed_new, text="lonmax", font=self.frame_font).grid(row=2, column=1)
        self.lonmax_entry = tkinter.Entry(self.seed_new)
        self.lonmax_entry.insert(0, "180")
        self.lonmax_entry.grid(row=3, column=1)
        
        latmin_label = tkinter.Label(self.seed_new, text="latmin", font=self.frame_font).grid(row=2, column=2)
        self.latmin_entry = tkinter.Entry(self.seed_new)
        self.latmin_entry.insert(0, "-90")
        self.latmin_entry.grid(row=3, column=2)

        latmax_label = tkinter.Label(self.seed_new, text="latmax", font=self.frame_font).grid(row=2, column=3)
        self.latmax_entry = tkinter.Entry(self.seed_new)
        self.latmax_entry.insert(0, "90")
        self.latmax_entry.grid(row=3, column=3)

        dx_label = tkinter.Label(self.seed_new, text="dx", font=self.frame_font).grid(row=2, column=4)
        self.dx_entry = tkinter.Entry(self.seed_new)
        self.dx_entry.insert(0, "1")
        self.dx_entry.grid(row=3, column=4)

        dy_label = tkinter.Label(self.seed_new, text="dy", font=self.frame_font).grid(row=2, column=5)
        self.dy_entry = tkinter.Entry(self.seed_new)
        self.dy_entry.insert(0, "1")
        self.dy_entry.grid(row=3, column=5)

        poc_label = tkinter.Label(self.seed_new, text="poc_type", font=self.frame_font).grid(row=4, column=0)
        self.pocvar = tkinter.StringVar()
        self.poctype = ["Eg", "M"]
        self.pocvar.set("M")
        self.poctype = tkinter.OptionMenu(self.seed_new, self.pocvar, *self.poctype)
        self.poctype.grid(row=5, column=0)

        phylum_label = tkinter.Label(self.seed_new, text="phylum", font=self.frame_font).grid(row=4, column=1)
        self.phylavar = tkinter.StringVar()
        self.phylum = ["cnidaria", "ctenophora", "chordata"]
        self.phylavar.set("cnidaria")
        self.phylum = tkinter.OptionMenu(self.seed_new, self.phylavar, *self.phylum)
        self.phylum.grid(row=5, column=1)

        bathymetry_label = tkinter.Label(self.seed_new, text="bathymetry", font=self.frame_font).grid(row=4, column=2)
        self.bathymetry_entry = tkinter.Entry(self.seed_new)
        self.bathymetry_entry.insert(0, os.path.join(self.supplementary_data_dir, "etopo2.nc"))
        self.bathymetry_entry.grid(row=5, column=2)

        m0_label = tkinter.Label(self.seed_new, text="initial mass", font=self.frame_font).grid(row=4, column=3)
        self.m0_entry = tkinter.Entry(self.seed_new)
        self.m0_entry.insert(0, os.path.join(self.supplementary_data_dir, "Luo_M_Eg_biome_data.json"))
        self.m0_entry.grid(row=5, column=3)

        area_label = tkinter.Label(self.seed_new, text="area", font=self.frame_font).grid(row=4, column=4)
        self.area_entry = tkinter.Entry(self.seed_new)
        self.area_entry.insert(0, os.path.join(self.supplementary_data_dir, "area_grid.npy"))
        self.area_entry.grid(row=5, column=4)

        biome_label = tkinter.Label(self.seed_new, text="biome", font=self.frame_font).grid(row=4, column=5)
        self.biome_entry = tkinter.Entry(self.seed_new)
        self.biome_entry.insert(0, os.path.join(self.supplementary_data_dir, "biomegrid2.npy"))
        self.biome_entry.grid(row=5, column=5)

        outfile_label = tkinter.Label(self.seed_new, text="outfile", font=self.frame_font).grid(row=6, column=0)
        self.outfile_entry = tkinter.Entry(self.seed_new, state="disabled")
        self.outfile_entry.grid(row=7, column=0)
        

        self.seed_start = tkinter.Frame(self.simseed, relief=tkinter.FLAT, pady=25, padx=25)
        self.seed_start.grid(row=2, column=0, rowspan=1)
        
        tkinter.Button(self.seed_start, text="SEED", bg='green', command=self.on_click_seed).grid(row=0)

    def setup_reader_frame(self):
        self.required_readers = tkinter.Frame(self.simreader, bg='lightgray', bd=2,
                            relief=tkinter.SUNKEN, pady=25, padx=25)
        self.required_readers.grid(row=0, column=0, rowspan=1, padx=20, pady=(20,0), sticky="we")
        
        req_readers_title = tkinter.Label(self.required_readers, text="Required Readers", 
                                             font=self.frame_subtitle_font).grid(row=0, pady=(0, 25))
        self.LANDMASK = tkinter.BooleanVar()
        self.LANDMASK.set(True)
        land_mask_button = tkinter.Checkbutton(self.required_readers, text = "OpenDrift LandMask",
                                                            variable=self.LANDMASK  ).grid(row=1, column=0, sticky="w")
        tmp_label = tkinter.Label(self.required_readers, text="Temperature:", font=self.frame_font).grid(row=2, column=0, sticky="w")
        self.tmp_entry = tkinter.Entry(self.required_readers, width=100)
        self.tmp_entry.insert(0, os.path.join(self.supplementary_data_dir, "tmp_luo.nc"))
        self.tmp_entry.grid(row=2, column=1, sticky="we")

        bmt_label = tkinter.Label(self.required_readers, text="Bathymetry:", font=self.frame_font).grid(row=3, column=0, sticky="w")
        self.bmt_entry = tkinter.Entry(self.required_readers, width=100)
        self.bmt_entry.insert(0, os.path.join(self.supplementary_data_dir, "etopo2.nc"))
        self.bmt_entry.grid(row=3, column=1, sticky="we")


        self.optional_readers = tkinter.Frame(self.simreader, bg='lightgray', bd=2,
                            relief=tkinter.SUNKEN, pady=25, padx=25)
        self.optional_readers.grid(row=1, column=0, rowspan=1, padx=20, pady=(20,0), sticky="we")
        
        opt_readers_title = tkinter.Label(self.optional_readers, text="Optional Readers", 
                                             font=self.frame_subtitle_font).grid(row=0, pady=(0, 25), columnspan=2)
        self.opt_reader_box = tkinter.Text(self.optional_readers)
        self.opt_reader_box.grid(row=1, column=0)
        self.opt_reader_box.insert(tkinter.END, "FORMAT:\nfilepath1\nfilepath2\n...")
        
    def setup_configure_frame(self):
        self.required_configs = tkinter.Frame(self.simconfigure, bg='lightgray', bd=2,
                            relief=tkinter.SUNKEN, pady=25, padx=25)
        self.required_configs.grid(row=0, column=0, rowspan=1, padx=20, pady=(20,0), sticky="we")
        
        req_configs_title = tkinter.Label(self.required_configs, text="Required Configures", 
                                             font=self.frame_subtitle_font).grid(row=0, pady=(0, 25))
        self.OCEANONLY = tkinter.BooleanVar()
        self.OCEANONLY.set(True)
        ocean_only_button = tkinter.Checkbutton(self.required_configs, text = "oceanonly (recommended on)",
                                                            variable=self.OCEANONLY).grid(row=1, column=0, sticky="w")
        
        self.AUTOLANDMASK = tkinter.BooleanVar()
        land_mask_button = tkinter.Checkbutton(self.required_configs, text = "auto landmask",
                                                            variable=self.AUTOLANDMASK).grid(row=1, column=1, sticky="w")

        solver_label = tkinter.Label(self.required_configs, text="advection_scheme", font=self.frame_font).grid(row=2, column=0, sticky="w")
        self.solvervar = tkinter.StringVar()
        self.solver = ["runge-kutta", "runge-kutta4", "euler"]
        self.solvervar.set("runge-kutta")
        self.solver = tkinter.OptionMenu(self.required_configs, self.solvervar, *self.solver)
        self.solver.grid(row=2, column=1)


        self.optional_configs = tkinter.Frame(self.simconfigure, bg='lightgray', bd=2,
                            relief=tkinter.SUNKEN, pady=25, padx=25)
        self.optional_configs.grid(row=1, column=0, rowspan=1, padx=20, pady=(20,0), sticky="we")
        
        opt_configs_title = tkinter.Label(self.optional_configs, text="Optional Configures", 
                                             font=self.frame_subtitle_font).grid(row=0, pady=(0, 25), sticky="w")
        self.opt_config_box = tkinter.Text(self.optional_configs)
        self.opt_config_box.grid(row=1, column=0)
        self.opt_config_box.insert(tkinter.END, "FORMAT:\nconfigure1Name configure1Value\nconfigure2Name configure2Value\n...")
        
    def setup_simulation_type_frame(self):
        self.simulationparams = tkinter.Frame(self.simtype, bg='lightgray', bd=2,
                            relief=tkinter.SUNKEN, pady=25, padx=25)
        self.simulationparams.grid(row=0, column=0, rowspan=1, padx=20, pady=(20,0), sticky="we")
        
        simtype_label = tkinter.Label(self.simulationparams, text="type", font=self.frame_font).grid(row=0, column=0, sticky="w", pady=(0, 10))
        self.typevar = tkinter.StringVar()
        self.type = ["normal", "grid"]
        self.typevar.set("normal")
        self.type = tkinter.OptionMenu(self.simulationparams, self.typevar, *self.type, command=self.on_check_sim_type)
        self.type.grid(row=0, column=1, padx=10, pady=(0, 10))

        sf_label = tkinter.Label(self.simulationparams, text="split factor:", font=self.frame_font).grid(row=0, column=2, sticky="w", pady=(0, 10))
        self.sf_entry = tkinter.Entry(self.simulationparams, state="disabled")
        self.sf_entry.grid(row=0, column=3, sticky="w", padx=10, pady=(0, 10))

        self.ADVECTION = tkinter.BooleanVar()
        advection_button = tkinter.Checkbutton(self.simulationparams, text = "horizontal advection",
                                                            variable=self.ADVECTION).grid(row=0, column=4, sticky="w", padx=10, pady=(0, 10))
        
        self.FRAGMENTATION = tkinter.BooleanVar()
        fragmentation_button = tkinter.Checkbutton(self.simulationparams, text = "fragmentation", command=self.on_check_fragmentation,
                                                            variable=self.FRAGMENTATION).grid(row=0, column=5, sticky="w", padx=10, pady=(0, 10))

        start_label = tkinter.Label(self.simulationparams, text="start time:", font=self.frame_font).grid(row=1, column=0, sticky="w", pady=(0, 10))
        self.start_entry = tkinter.Entry(self.simulationparams)
        self.start_entry.insert(0, "YYYY-MM-DD-h")
        self.start_entry.grid(row=1, column=1, sticky="w", padx=10, pady=(0, 10))

        dt_label = tkinter.Label(self.simulationparams, text="time step:", font=self.frame_font).grid(row=1, column=2, sticky="w", pady=(0, 10))
        self.dt_entry = tkinter.Entry(self.simulationparams)
        self.dt_entry.insert(0, "h:m:s")
        self.dt_entry.grid(row=1, column=3, sticky="w", padx=10, pady=(0, 10))

        dto_label = tkinter.Label(self.simulationparams, text="output time step:", font=self.frame_font).grid(row=1, column=4, sticky="w", pady=(0, 10))
        self.dto_entry = tkinter.Entry(self.simulationparams)
        self.dto_entry.insert(0, "h:m:s")
        self.dto_entry.grid(row=1, column=5, sticky="w", padx=10, pady=(0, 10))

        steps_label = tkinter.Label(self.simulationparams, text="steps:", font=self.frame_font).grid(row=1, column=6, sticky="w", pady=(0, 10))
        self.step_entry = tkinter.Entry(self.simulationparams)
        self.step_entry.grid(row=1, column=7, sticky="w", padx=10, pady=(0, 10))

        w0_label = tkinter.Label(self.simulationparams, text="init vertical\n velocity:", font=self.frame_font).grid(row=2, column=0, sticky="w", pady=(0, 10))
        self.w0_entry = tkinter.Entry(self.simulationparams)
        self.w0_entry.grid(row=2, column=1, sticky="w", padx=10, pady=(0, 10))

        w0type_label = tkinter.Label(self.simulationparams, text="velocity type", font=self.frame_font).grid(row=2, column=2, sticky="w", pady=(0, 10))
        self.w0var = tkinter.StringVar()
        self.w0type = ["variable", "constant"]
        self.w0var.set("variable")
        self.w0type = tkinter.OptionMenu(self.simulationparams, self.w0var, *self.w0type)
        self.w0type.grid(row=2, column=3, padx=10, pady=(0, 10))

        decaytype_label = tkinter.Label(self.simulationparams, text="decay rate", font=self.frame_font).grid(row=2, column=4, sticky="w", pady=(0, 10))
        self.decayvar = tkinter.StringVar()
        self.decaytype = ["exp", "linear"]
        self.decayvar.set("exp")
        self.decaytype = tkinter.OptionMenu(self.simulationparams, self.decayvar, *self.decaytype)
        self.decaytype.grid(row=2, column=5, padx=10, pady=(0, 10))

        mdecaytype_label = tkinter.Label(self.simulationparams, text="decay", font=self.frame_font).grid(row=2, column=6, sticky="w", pady=(0, 10))
        self.mdecayvar = tkinter.StringVar()
        self.mdecaytype = ["mass", "area"]
        self.mdecayvar.set("mass")
        self.mdecaytype = tkinter.OptionMenu(self.simulationparams, self.mdecayvar, *self.mdecaytype)
        self.mdecaytype.grid(row=2, column=7, padx=10, pady=(0, 10))

        out_label = tkinter.Label(self.simulationparams, text="outfile:", font=self.frame_font).grid(row=3, column=0, sticky="w", pady=(0, 10))
        self.simout_entry = tkinter.Entry(self.simulationparams)
        self.simout_entry.insert(0, "netCDF file")
        self.simout_entry.grid(row=3, column=1, sticky="w", padx=10, pady=(0, 10))

        ev_label = tkinter.Label(self.simulationparams, text="export\nvariables", font=self.frame_font).grid(row=3, column=2, sticky="w", pady=(0, 10))
        self.ev_list = tkinter.Listbox(self.simulationparams, selectmode=tkinter.MULTIPLE)
        valid_vars = ["trajectory", "time", "status", "moving", "age_seconds", "origin_marker",
                  "lon", "lat", "z", "wind_drift_factor", "current_drift_factor",
                  "terminal_velocity", "mass", "x_sea_water_velocity", "y_sea_water_velocity",
                  "land_binary_mask", "sea_water_temperature", "depth"]
        for i, var in enumerate(valid_vars):
            self.ev_list.insert(i + 1, var)
            if var in ["trajectory", "time", "status", "age_seconds", "lon", "lat", "z", "mass"]:
                self.ev_list.selection_set(i)
        self.ev_list.grid(row=3, column=3, sticky="w", padx=10, pady=(0, 10))

        fragfunc_label = tkinter.Label(self.simulationparams, text="fragmentation\nfunction:", font=self.frame_font).grid(row=3, column=4, sticky="w", pady=(0, 10))
        self.ffunc_entry = tkinter.Entry(self.simulationparams, state="disabled")
        self.ffunc_entry.grid(row=3, column=5, sticky="w", padx=10, pady=(0, 10))

        self.make_object = tkinter.Frame(self.simtype, relief=tkinter.FLAT, pady=25, padx=25)
        self.make_object.grid(row=1, column=0, rowspan=1)
        
        tkinter.Button(self.make_object, text="INITIALIZE", bg='green', command=self.on_click_initialize).grid(row=0)

    def on_check_save_seed(self):
        if self.SAVESEED.get():
            self.outfile_entry["state"] ="normal"
        else:
            self.outfile_entry["state"] ="disabled"
    
    def on_check_relative_mass(self):
        if self.RELATIVEMASS.get():
            self.m0_entry["state"]="disabled"
        else:
            self.m0_entry["state"] = "normal"
    
    def on_check_fragmentation(self):
        if self.FRAGMENTATION.get():
            self.ffunc_entry["state"] = "normal"
        else:
            self.ffunc_entry["state"] = "disabled"
    
    def on_check_sim_type(self, event):
        if event == "normal":
            self.sf_entry["state"]="disabled"
        else:
            self.sf_entry["state"]="normal"

    def on_click_initialize(self):
        self.update_idletasks()
        if self.typevar.get() == "normal":
            self.initialize_normal_run()
        else:
            self.initialize_grid_run()

    def disable_new_seed(self, event):
        if len(self.seed_file_entry.get()) > 1:
            self.set_frame_state(self.seed_from_file, "normal")
            self.set_frame_state(self.seed_new, "disable")
        else:
            self.set_frame_state(self.seed_from_file, "normal")
            self.set_frame_state(self.seed_new, "normal")

    def set_frame_state(self, frame, state):
        for widget in frame.winfo_children():
            widget.config(state=state)
        if frame==self.seed_new:
            if state=="disable":
                pass
            else:
                self.on_check_relative_mass()
                self.on_check_save_seed()
    
    def launch_error_window(self, error):
        self.error_window = tkinter.Toplevel(self)
        self.error_window.title(type(self.error).__name__)
        self.error_window.geometry("250x100")
        self.error_label = tkinter.Label(self.error_window, text=error)
        self.error_label.pack(pady=5)
    
    def launch_progress_window(self, title):
        self.progress_window = tkinter.Toplevel(self)
        self.progress_window.title("Progress")
        self.progress_window.geometry("250x100")
        self.progress_window.resizable(False, False)
        self.progress_label = tkinter.Label(self.progress_window, text=title)
        self.progress_label.pack(pady=5)
        self.progress = ttk.Progressbar(self.progress_window, mode="indeterminate")
        self.progress.pack(pady=5, padx=20, fill="x")
        self.progress.start(10)

    def on_click_seed(self):
        self.update_idletasks()
        self.launch_progress_window("Seeding")
        self.thread = threading.Thread(target=self.seed, daemon=True)
        self.thread.start()
        self.after(100, self.check_thread_finished)
    
    def check_thread_finished(self):
        if self.thread.is_alive():
            self.after(100, self.check_thread_finished)    
    
    def seed(self):
        if self.lonmin_entry["state"] == "disabled":
            entry = self.seed_file_entry.get()
            if len(entry)==0:
                self.error = ValueError("File path is empty.")
                self.progress_window.destroy()
                self.launch_error_window(self.error)
            else:
                try:
                    self.sdata = SeedFromFile(entry)
                except FileNotFoundError:
                    self.error = FileNotFoundError(f"File {entry} does not exist.")
                    self.launch_error_window(self.error)
        else:
            latmin = self.latmin_entry.get()
            lonmin = self.lonmin_entry.get()
            lonmax = self.lonmax_entry.get()
            latmax = self.latmax_entry.get()
            dx = self.dx_entry.get()
            dy = self.dy_entry.get()
            b = self.bathymetry_entry.get()
            area = self.area_entry.get()
            m0 = self.m0_entry.get()
            biome = self.biome_entry.get()
            poctype = self.pocvar.get()
            phylum = self.phylavar.get()
            outfile = None if self.outfile_entry.get()=="" else self.outfile_entry.get()
            data = {"latmin":latmin, "lonmin":lonmin, "lonmax":lonmax, "latmax":latmax, "dx":dx, "dy":dy,
                    "bathymetry":b, "area":area, "initial_mass":m0, "biome":biome}
            for item_name, item in data.items():
                if item_name=="initial mass" and self.RELATIVEMASS.get():
                    data[item_name] = None
                elif len(item) == 0:
                    self.error = ValueError(f"{item_name} field is empty.")
                    self.launch_error_window(self.error)
                else:
                    if item_name in ["initial_mass", "biome", "area", "bathymetry"]:
                        pass
                    else:
                        try:
                            data[item_name] = float(item)
                        except ValueError:
                            self.error = ValueError(f"{item} in {item_name} is not a valid float format.")
                            self.launch_error_window(self.error)
            try:
                self.sdata = Seed(latmin=data["latmin"], latmax=data["latmax"], lonmin=data["lonmin"], lonmax=data["lonmax"],
                                  dx=data["dx"], dy=data["dy"], bathymetrypath=data["bathymetry"], poctype=poctype, phylum=phylum,
                                  areapath=data["area"], biomegridpath=data["biome"], initialmassdata=data["initial_mass"],outfile=outfile)
            except Exception as e:
                self.launch_error_window(e)

        self.progress_window.destroy()

    def add_readers(self):
        readers = []
        bat = reader_netCDF_CF_generic.Reader(self.bmt_entry.get(), standard_name_mapping={"topo":"depth"})
        readers.append(bat)
        tmp = reader_netCDF_CF_generic.Reader(self.tmp_entry.get())
        tmp.verticalbuffer = 100
        tmp.always_valid = True
        readers.append(tmp)
        if self.LANDMASK.get() and not self.AUTOLANDMASK.get():
            readers.append(reader_global_landmask.Reader())
        opt_readers = self.opt_reader_box.get("1.0",'end-1c').split("\n")
        if "FORMAT" not in opt_readers[0]:
            for reader in opt_readers:
                rd = reader_netCDF_CF_generic.Reader(reader)
                readers.append(rd)
        return readers

    def add_configures(self):
        config = {}
        config['drift:advection_scheme'] = self.solvervar.get()
        config['general:use_auto_landmask'] = self.AUTOLANDMASK.get()
        config['seed:ocean_only']=self.OCEANONLY.get()
        opt_configs = self.opt_config_box.get("1.0",'end-1c').split("\n")
        if "FORMAT" not in opt_configs[0]:
            for opt in opt_configs:
                configname, configvalue = opt.split(" ")
                config[configname] = configvalue
        return config 

    def initialize_normal_run(self):
        params = {"m0":self.sdata.mass, "loglevel":0, "decay_type":self.decayvar.get(), "initial_velocity":float(self.w0_entry.get()),
                  "vertical_velocity_type": self.w0var.get(), "ffunction":self.ffunc_entry.get()}
        if self.mdecayvar.get() == "mass":
            self.obj = mcd.CarbonDrift(**params)
        else:
            self.obj = acd.CarbonDrift(**params)
        if not self.ADVECTION.get():
            self.obj.deactivate_horizontal_advection()
        if not self.FRAGMENTATION.get():
            self.obj.deactivate_fragmentation()
        self.obj.add_reader(self.add_readers())
        for key, val in self.add_configures().items():
            self.obj.set_config(key, val)
        if hasattr(self.sdata, "origin_marker"):
            origin_marker = self.sdata.origin_marker
        if hasattr(self.sdata, "biome"):
            origin_marker = self.sdata.biome
        else:
            origin_marker = 0
        
        self.obj.seed_elements(lon=self.sdata.lon, lat = self.sdata.lat, z=self.sdata.z,
                                   mass = self.sdata.mass, time=self.valid_date(), origin_marker = origin_marker)
        self.obj.run(time_step=timedelta(minutes=30), time_step_output=timedelta(hours=1), steps=500, outfile="home/peharicc/Documents/test.nc")
        
    def valid_date(self):
        try:
            return datetime.strptime(self.start_entry.get(), "%Y-%m-%d-%H")
        except ValueError:
            raise ValueError(f"not a valid date: {self.start_entry.get()!r}")

if __name__ == "__main__":
    GUI = CarbonDriftGUI()
    GUI.mainloop()