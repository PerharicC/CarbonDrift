import tkinter.messagebox
import tkinter.scrolledtext
import numpy as np
import matplotlib.pyplot as plt
from PIL import ImageTk, Image
import tkinter
from tkinter import ttk, filedialog
import threading
import time
import os
import logging
from PIL import Image, ImageTk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from datetime import datetime, timedelta
from opendrift.readers import reader_global_landmask
from opendrift.readers import reader_netCDF_CF_generic
import carbondrift
import carbondrift.models.massdecay.gridrun as mcdgrid
import carbondrift.models.areadecay.gridrun as acdgrid
import carbondrift.models.massdecay.carbondrift as mcd
import carbondrift.models.areadecay.carbondrift as acd
from carbondrift.simulation.seeding import *
from carbondrift.simulation.run import valid_date, valid_fragmentation_function, valid_time_step

from carbondrift.models.logger import Logger
from carbondrift.models import plots
from carbondrift.plotting.plot_run import configure_locations

log = Logger("carbondrift.carbondrift_gui")
logger = log.LOGGER


class MissingDataError(Exception):
    def __init__(self, msg="Please provide data."):
        self.msg = msg
        super().__init__(self.msg)

    def __str__(self):
        return f'{self.msg}'

class LoggerHandler(logging.Handler):
    def __init__(self, text_widget):
        logging.Handler.__init__(self)
        self.text_widget = text_widget

    def emit(self, record):
        log_entry = self.format(record) + "\n"
        self.text_widget.after(0, self.append_log, log_entry)

    def append_log(self, log_entry):
        self.text_widget.config(state="normal")
        self.text_widget.insert(tkinter.END, log_entry)
        self.text_widget.see(tkinter.END)  # Auto-scroll to the latest entry
        self.text_widget.config(state="disabled")

class CarbonDriftGUI(tkinter.Tk):
    def __init__(self):
        tkinter.Tk.__init__(self)
        self.title("CarbonDrift " + carbondrift.__version__ + " GUI")
        style = ttk.Style()
        style.theme_use("clam")
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
        self.setup_simulation_run_frame()
        self.setup_plotting_frame()

        self.seed_file_entry = ""
        self.SEEDCLICKED=False
        self.INITIALIZECLICKED = False
        
        self.saved_ev_selection = set()
        self.ev_list.bind("<FocusOut>", self.save_ev_selection)
        self.ev_list.bind("<FocusIn>", self.restore_ev_selection)

        self.thread = None
        self.sdata = None
        self.tmp_entry=None
        self.bmt_entry = None
        self.opt_readers = []
        self.obj = None
        self.simout_entry = None
        self.oplot_entry = None
        self.plotfilesplus = []
        self.plotfilesminus = []
        self.SIMDISPLAYMAXWIDTH = 800
        self.SIMDISPLAYMAXHEIGHT = 800

        style.configure(
            "Success.Horizontal.TProgressbar",
            troughcolor="white",
            background="green",
            bordercolor="gray",
            lightcolor="lightgreen",
            darkcolor="darkgreen",
            pbarrelief = "groove"
        )

    def setup_seed_frame(self):
        self.seed_from_file = tkinter.Frame(self.simseed, bg='lightgray', bd=2,
                            relief=tkinter.SUNKEN, pady=25, padx=25)
        self.seed_from_file.grid(row=0, column=0, rowspan=1, padx=20, pady=(20,0), sticky="w")
        
        seed_from_file_title = tkinter.Label(self.seed_from_file, text="Seed from existing file", 
                                             font=self.frame_subtitle_font).grid(row=0, pady=(0, 25), columnspan=2)
        
        self.seed_file_label = tkinter.Label(self.seed_from_file, text="Import seed:", font=self.frame_font)
        self.seed_file_label.grid(row=1, column=0, sticky="w")
        tkinter.Button(
            self.seed_from_file, text = "Browse Files", command = self.browse_seed
            ).grid(row = 1, column = 1, sticky="w", padx=10, pady=(0, 10))
        self.remove_seed_button = tkinter.Button(
            self.seed_from_file, text = "Remove Import", command = self.remove_seed, state="disabled",
            )
        self.remove_seed_button.grid(row = 1, column = 2, sticky="w", padx=10, pady=(0, 10))


        self.seed_new = tkinter.Frame(self.simseed, bg='lightgray', bd=2,
                            relief=tkinter.SUNKEN, pady=25, padx=25)
        self.seed_new.grid(row=1, column=0, rowspan=1, padx=20, pady=(20,0), sticky="we")
        seed_from_file_title = tkinter.Label(self.seed_new, text="Make New Seed", 
                                             font=self.frame_subtitle_font).grid(row=0, pady=(0, 25), columnspan=2, sticky="w")
        
        self.RELATIVEMASS = tkinter.BooleanVar()
        relative_mass_button = tkinter.Checkbutton(self.seed_new, text="relative mass", command=self.on_check_relative_mass,
                                                            variable=self.RELATIVEMASS).grid(row=1, column=0, sticky="w")
        
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
        

        self.seed_start = tkinter.Frame(self.simseed, relief=tkinter.FLAT, pady=25, padx=25)
        self.seed_start.grid(row=2, column=0, rowspan=1)
        
        tkinter.Button(self.seed_start, text="SEED", bg='green', command=self.on_click_seed).grid(row=0)
        self.reset_seed_button = tkinter.Button(self.seed_start, text="RESET", bg='yellow', command=self.on_click_reset_seed, state = "disabled")
        self.reset_seed_button.grid(row=0, column=1)

    def setup_reader_frame(self):
        self.required_readers = tkinter.Frame(self.simreader, bg='lightgray', bd=2,
                            relief=tkinter.SUNKEN, pady=25, padx=25)
        self.required_readers.grid(row=0, column=0, rowspan=1, padx=20, pady=(20,0), sticky="we")
        
        req_readers_title = tkinter.Label(self.required_readers, text="Required Readers", 
                                             font=self.frame_subtitle_font).grid(row=0, pady=(0, 25), sticky="w")
        self.LANDMASK = tkinter.BooleanVar()
        self.LANDMASK.set(True)
        land_mask_button = tkinter.Checkbutton(self.required_readers, text = "OpenDrift LandMask",
                                                            variable=self.LANDMASK).grid(row=1, column=0, pady=(0,10), sticky="w")
        
        self.tmp_label = tkinter.Label(self.required_readers, text="Import temperature:", font=self.frame_font)
        self.tmp_label.grid(row=2, column=0, padx=10, sticky="w")
        tkinter.Button(
            self.required_readers, text = "Browse Files", command = lambda:self.browse_required_readers("tmp")
            ).grid(row = 3, column = 0, sticky="w", padx=10, pady=(0, 10))
        
        self.bmt_label = tkinter.Label(self.required_readers, text="Import bathymetry:", font=self.frame_font)
        self.bmt_label.grid(row=2, column=1, padx=10, sticky="w")
        tkinter.Button(
            self.required_readers, text = "Browse Files", command = lambda:self.browse_required_readers("bmt")
            ).grid(row = 3, column = 1, sticky="w", padx=10, pady=(0, 10))

        self.optional_readers = tkinter.Frame(self.simreader, bg='lightgray', bd=2,
                            relief=tkinter.SUNKEN, pady=25, padx=25)
        self.optional_readers.grid(row=1, column=0, rowspan=1, padx=20, pady=(20,0), sticky="we")
        
        opt_readers_title = tkinter.Label(self.optional_readers, text="Optional Readers", 
                                             font=self.frame_subtitle_font).grid(row=0, pady=(0, 25), columnspan=2, sticky="w")
        self.opt_label = tkinter.Label(self.optional_readers, text="Import other raeders:", font=self.frame_font)
        self.opt_label.grid(row=1, column=0, sticky="w")
        tkinter.Button(
            self.optional_readers, text = "Browse Files", command = self.browse_optional_readers
            ).grid(row = 2, column = 0, sticky="w", padx=10, pady=(0, 10))
        
    def setup_configure_frame(self):
        self.required_configs = tkinter.Frame(self.simconfigure, bg='lightgray', bd=2,
                            relief=tkinter.SUNKEN, pady=25, padx=25)
        self.required_configs.grid(row=0, column=0, rowspan=1, padx=20, pady=(20,0), sticky="we")
        
        req_configs_title = tkinter.Label(self.required_configs, text="Required Configures", 
                                             font=self.frame_subtitle_font).grid(row=0, pady=(0, 25), sticky="w")
        
        self.OCEANONLY = tkinter.BooleanVar()
        self.OCEANONLY.set(True)
        ocean_only_button = tkinter.Checkbutton(self.required_configs, text = "oceanonly (recommended on)",
                                                            variable=self.OCEANONLY).grid(row=1, column=0, sticky="w")
        
        self.AUTOLANDMASK = tkinter.BooleanVar()
        land_mask_button = tkinter.Checkbutton(self.required_configs, text = "auto landmask",
                                                            variable=self.AUTOLANDMASK).grid(row=1, column=1, sticky="w")

        solver_label = tkinter.Label(self.required_configs, text="advection scheme:", font=self.frame_font).grid(row=1, column = 2, sticky="w", padx=(20, 0))
        self.solvervar = tkinter.StringVar()
        self.solver = ["runge-kutta", "runge-kutta4", "euler"]
        self.solvervar.set("runge-kutta")
        self.solver = tkinter.OptionMenu(self.required_configs, self.solvervar, *self.solver)
        self.solver.grid(row=2, column=2, pady=(10,0), padx=(20, 0), sticky="w")


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

        simtype_label = tkinter.Label(self.simulationparams, text="type", font=self.frame_font).grid(row=0, column=0)
        self.typevar = tkinter.StringVar()
        self.type = ["normal", "grid"]
        self.typevar.set("normal")
        self.type = tkinter.OptionMenu(self.simulationparams, self.typevar, *self.type, command=self.on_check_sim_type)
        self.type.grid(row=1, column=0, pady=(0, 20))

        sf_label = tkinter.Label(self.simulationparams, text="split factor:", font=self.frame_font).grid(row=0, column=1)
        self.sf_entry = tkinter.Entry(self.simulationparams, state="disabled")
        self.sf_entry.grid(row=1, column=1, pady=(0, 20))

        self.ADVECTION = tkinter.BooleanVar()
        advection_button = tkinter.Checkbutton(self.simulationparams, text = "horizontal advection",
                                                            variable=self.ADVECTION).grid(row=0, column=2, padx=10, sticky="w")
        
        self.FRAGMENTATION = tkinter.BooleanVar()
        fragmentation_button = tkinter.Checkbutton(self.simulationparams, text = "fragmentation", command=self.on_check_fragmentation,
                                                            variable=self.FRAGMENTATION).grid(row=1, column=2, pady=(0, 20), sticky="w", padx=10)

        fragfunc_label = tkinter.Label(self.simulationparams, text="fragmentation\nfunction:", font=self.frame_font).grid(row=0, column=3)
        self.ffunc_entry = tkinter.Entry(self.simulationparams, state="disabled")
        self.ffunc_entry.grid(row=1, column=3, pady=(0, 20))

        start_label = tkinter.Label(self.simulationparams, text="start time:", font=self.frame_font).grid(row=2, column=0)
        self.start_entry = tkinter.Entry(self.simulationparams)
        self.start_entry.insert(0, "YYYY-MM-DD-h")#"1993-01-01-0"
        self.start_entry.grid(row=3, column=0, pady=(0, 20))

        dt_label = tkinter.Label(self.simulationparams, text="time step:", font=self.frame_font).grid(row=2, column=1)
        self.dt_entry = tkinter.Entry(self.simulationparams)
        self.dt_entry.insert(0, "h:m:s")#"0:30:0"
        self.dt_entry.grid(row=3, column=1, pady=(0, 20))

        dto_label = tkinter.Label(self.simulationparams, text="output time step:", font=self.frame_font).grid(row=2, column=2)
        self.dto_entry = tkinter.Entry(self.simulationparams)
        self.dto_entry.insert(0, "h:m:s")#"1:0:0"
        self.dto_entry.grid(row=3, column=2, pady=(0, 20))

        steps_label = tkinter.Label(self.simulationparams, text="steps:", font=self.frame_font).grid(row=2, column=3)
        self.step_entry = tkinter.Entry(self.simulationparams)
        self.step_entry.insert(0, 500)
        self.step_entry.grid(row=3, column=3, pady=(0, 20))

        w0_label = tkinter.Label(self.simulationparams, text="init vertical\n velocity:", font=self.frame_font).grid(row=4, column=0)
        self.w0_entry = tkinter.Entry(self.simulationparams)
        self.w0_entry.insert(0, -0.01)
        self.w0_entry.grid(row=5, column=0, pady=(0, 20))

        w0type_label = tkinter.Label(self.simulationparams, text="velocity type", font=self.frame_font).grid(row=4, column=1)
        self.w0var = tkinter.StringVar()
        self.w0type = ["variable", "constant"]
        self.w0var.set("variable")
        self.w0type = tkinter.OptionMenu(self.simulationparams, self.w0var, *self.w0type)
        self.w0type.grid(row=5, column=1, pady=(0, 20))

        decaytype_label = tkinter.Label(self.simulationparams, text="decay rate", font=self.frame_font).grid(row=4, column=2)
        self.decayvar = tkinter.StringVar()
        self.decaytype = ["exp", "linear"]
        self.decayvar.set("exp")
        self.decaytype = tkinter.OptionMenu(self.simulationparams, self.decayvar, *self.decaytype)
        self.decaytype.grid(row=5, column=2, pady=(0, 20))

        mdecaytype_label = tkinter.Label(self.simulationparams, text="decay", font=self.frame_font).grid(row=4, column=3)
        self.mdecayvar = tkinter.StringVar()
        self.mdecaytype = ["mass", "area"]
        self.mdecayvar.set("mass")
        self.mdecaytype = tkinter.OptionMenu(self.simulationparams, self.mdecayvar, *self.mdecaytype)
        self.mdecaytype.grid(row=5, column=3, pady=(0, 20))
        

        output_frame_ev = tkinter.Frame(self.simulationparams, relief=tkinter.FLAT)
        output_frame_ev.grid(row=6, column=0, rowspan=1, columnspan=2)
        output_frame_of = tkinter.Frame(self.simulationparams, relief=tkinter.FLAT)
        output_frame_of.grid(row=6, column=2, rowspan=1, columnspan=2)

        self.out_label = tkinter.Label(output_frame_of, text="outfile:", font=self.frame_font)
        self.out_label.grid(row=0, column=0)
        tkinter.Button(
            output_frame_of, text = "Select Location", command = self.save_simulation
            ).grid(row=1, column=0, pady=(0, 20))

        ev_label = tkinter.Label(output_frame_ev, text="export variables:", font=self.frame_font).grid(row=0, column=0)
        self.ev_list = tkinter.Listbox(output_frame_ev, selectmode=tkinter.MULTIPLE, exportselection=False)
        valid_vars = ["trajectory", "time", "status", "moving", "age_seconds", "origin_marker",
                  "lon", "lat", "z", "wind_drift_factor", "current_drift_factor",
                  "terminal_velocity", "mass", "x_sea_water_velocity", "y_sea_water_velocity",
                  "land_binary_mask", "sea_water_temperature", "depth"]
        self.current_ev_selection = ["trajectory", "time", "status", "age_seconds", "lon", "lat", "z", "mass", "sea_water_temperature"]
        for i, var in enumerate(valid_vars):
            self.ev_list.insert(i + 1, var)
            if var in self.current_ev_selection:
                self.ev_list.selection_set(i)
        self.ev_list.grid(row=1, column=0, pady=(0, 20))

        self.make_object = tkinter.Frame(self.simtype, relief=tkinter.FLAT, pady=25, padx=25)
        self.make_object.grid(row=1, column=0, rowspan=1)
        
        tkinter.Button(self.make_object, text="INITIALIZE", bg='green', command=self.on_click_initialize).grid(row=0)
    
    def setup_simulation_run_frame(self):
        self.simulationlog = tkinter.Frame(self.run, bg='lightgray', bd=2,
                            relief=tkinter.SUNKEN, pady=25, padx=25, width=100)
        self.simulationlog.grid(row=0, column=0, rowspan=1, padx=20, pady=(20,0))
        
        simlog_label = tkinter.Label(self.simulationlog, text="Output Log", font=self.frame_subtitle_font).grid(row=0, column=0, pady=(0, 10))
        self.log_text = tkinter.scrolledtext.ScrolledText(self.simulationlog, state="disabled")
        self.log_text.grid(row=1, column=0, padx=10, pady=(0, 10))
        self.logger = logger
        text_handler = LoggerHandler(self.log_text)
        text_handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)-7s: %(message)s', "%H:%M:%S"))
        logging.getLogger().addHandler(text_handler)
        tkinter.Button(self.simulationlog, text="RUN", bg='green', command=self.on_click_run).grid(row=3)
    
    def setup_plotting_frame(self):
        self.plotparams = tkinter.Frame(self.plot, bg='lightgray', bd=2,
                            relief=tkinter.SUNKEN, pady=25, padx=25, width=100)
        self.plotparams.grid(row=0, column=0, rowspan=1, padx=20, pady=(20,0))
        
        plotmethod_label = tkinter.Label(self.plotparams, text="method", font=self.frame_font).grid(row=0, column=0)
        self.methodvar = tkinter.StringVar()
        self.method = np.sort(["mass_map", "current_strength", "drifter_properties", "mass_flux_map", "animate_3D",
                       "mean_lon_mass_flux", "mass_flux_distribution", "animate_current_3D", "vertical_particle_distribution"])
        self.methodvar.set("mass_flux_map")
        self.method = tkinter.OptionMenu(self.plotparams, self.methodvar, *self.method, command=self.on_check_plot_method)
        self.method.grid(row=1, column=0)

        self.simulation_label_plus = tkinter.Label(self.plotparams, text="upload simulations", font=self.frame_font)
        self.simulation_label_plus.grid(row=0, column=1)
        simulation_explore_plus = tkinter.Button(self.plotparams, text = "Browse Files", command = lambda:self.browse_simulations("plus")).grid(row = 1, column = 1)
        
        self.DIFF = tkinter.BooleanVar()
        difference_button = tkinter.Checkbutton(self.plotparams, text="diff", command=self.on_click_diff,
                                                            variable=self.DIFF).grid(row=0, column=2, sticky="S")
        
        self.ADD = tkinter.BooleanVar()
        add_button = tkinter.Checkbutton(self.plotparams, text="add",
                                                            variable=self.ADD).grid(row=1, column=2)
        
        self.ABS = tkinter.BooleanVar()
        add_button = tkinter.Checkbutton(self.plotparams, text="abs",
                                                            variable=self.ABS).grid(row=2, column=2, sticky="N", pady=(0,20))
        
        depth_label = tkinter.Label(self.plotparams, text="depth:", font=self.frame_font).grid(row=0, column=4)
        self.depth_entry= tkinter.Entry(self.plotparams)
        self.depth_entry.grid(row=1, column=4)

        self.simulation_label_minus = tkinter.Label(self.plotparams, text="upload simulations\nwith minus sign", font=self.frame_font, state="disabled")
        self.simulation_label_minus.grid(row=0, column=3)
        self.simulation_explore_minus = tkinter.Button(self.plotparams, text = "Browse Files", command = lambda:self.browse_simulations("minus"), state="disabled")
        self.simulation_explore_minus.grid(row = 1, column = 3)


        row2 = tkinter.Frame(self.plotparams, relief=tkinter.FLAT, bg="lightgray")
        row2.grid(row=3, column=0, columnspan=5, sticky="we")
        for i in range(6):  # Total columns
            row2.columnconfigure(i, weight=1)
        self.p1_label = tkinter.Label(row2, text="x property:", font=self.frame_font, state="disabled")
        self.p1_label.grid(row=0, column=0)
        self.p1var = tkinter.StringVar()
        valid_vars = np.sort(["mass", "lon", "lat", "z", "time", "x_sea_water_velocity",
                  "y_sea_water_velocity", "sea_water_temperature", "tau", "total_horizontal_velocity"])
        self.p1var.set("mass")
        self.p1 = tkinter.OptionMenu(row2, self.p1var, *valid_vars)
        self.p1.config(state="disabled")
        self.p1.grid(row=1, column=0, pady=(0, 20))

        self.p2_label = tkinter.Label(row2, text="y property:", font=self.frame_font, state="disabled")
        self.p2_label.grid(row=0, column=2)
        self.p2var = tkinter.StringVar()
        self.p2var.set("z")
        self.p2 = tkinter.OptionMenu(row2, self.p2var, *valid_vars)
        self.p2.config(state="disabled")
        self.p2.grid(row=1, column=2, pady=(0, 20))

        self.loc_label = tkinter.Label(row2, text="locations:", font=self.frame_font, state="disabled")
        self.loc_label.grid(row=0, column=4)
        self.loc_entry= tkinter.Entry(row2, state="disabled")
        self.loc_entry.insert(0, "lon1:lat1,lon2:lat2...")
        self.loc_entry.grid(row=1, column=4, pady=(0, 20))

        row3 = tkinter.Frame(self.plotparams, relief=tkinter.FLAT, bg="lightgray")
        row3.grid(row=4, column=0, columnspan=5, sticky="we")
        for i in range(6):  # Total columns
            row3.columnconfigure(i, weight=1)

        title_label = tkinter.Label(row3, text="title:", font=self.frame_font).grid(row = 0, column=0)
        self.title_entry= tkinter.Entry(row3)
        self.title_entry.grid(row=1, column=0, pady = (0, 20))

        xlabel_label = tkinter.Label(row3, text="xlabel:", font=self.frame_font).grid(row = 0, column=2)
        self.xlabel_entry= tkinter.Entry(row3)
        self.xlabel_entry.grid(row=1, column=2, pady=(0, 20))

        ylabel_label = tkinter.Label(row3, text="ylabel:", font=self.frame_font).grid(row = 0, column=4)
        self.ylabel_entry= tkinter.Entry(row3)
        self.ylabel_entry.grid(row = 1, column = 4, pady = (0, 20))

        cblabel_label = tkinter.Label(row3, text="colorbar label:", font=self.frame_font).grid(row=0, column=6)
        self.cblabel_entry= tkinter.Entry(row3)
        self.cblabel_entry.grid(row=1, column=6, pady=(0, 20))

        self.outfile_plot_label = tkinter.Label(row2, text="outfile:", font=self.frame_font)
        self.outfile_plot_label.grid(row=0, column= 6)
        tkinter.Button(
            row2, text = "Select Location", command = self.save_plot
            ).grid(row=1, column=6, pady=(0, 20))
        
        row4 = tkinter.Frame(self.plotparams, relief=tkinter.FLAT, bg="lightgray")
        row4.grid(row=5, column=0, columnspan=5, sticky="we")
        for i in range(6):  # Total columns
            row4.columnconfigure(i, weight=1)
        
        tkinter.Button(row4, text="PLOT", bg='green', command=self.on_click_plot).grid(row=0, column = 2)
        self.viewbutton = tkinter.Button(row4, text="VIEW", bg='green', command=self.on_click_view_plot, state="disabled")
        self.viewbutton.grid(row=0, column = 3)
    
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
        if self.SEEDCLICKED:
            self.update_idletasks()
            try:
                if self.typevar.get() == "normal":
                    self.initialize_normal_run()
                else:
                    self.initialize_grid_run()
                self.INITIALIZECLICKED = True
            except Exception as e:
                self.launch_error_window(e)
                self.reset_object_init()
        else:
            self.launch_error_window(MissingDataError("Please seed data first."))

    def on_click_run(self):
        if not self.SEEDCLICKED:
            self.error = MissingDataError("Please seed data first.")
            self.launch_error_window(self.error)
        elif not self.INITIALIZECLICKED:
            self.error = MissingDataError("Please initialize Simulation first.")
            self.launch_error_window(self.error)
        else:
            try:
                run_params = {"time_step":valid_time_step(self.dt_entry.get()), "time_step_output":valid_time_step(self.dto_entry.get()),
                            "steps":int(self.step_entry.get()), "outfile":self.simout_entry,
                            "export_variables":[self.ev_list.get(i) for i in self.ev_list.curselection()]}
            except Exception as e:
                self.INITIALIZECLICKED = False
                self.reset_object_init()
                self.launch_error_window(e)
                return
            
            self.run_progress = ttk.Progressbar(
                self.simulationlog, mode="determinate",
                orient="horizontal", length=200, style='Success.Horizontal.TProgressbar'
                                                )
            self.run_progress.grid(row = 2, column=0, padx=10, pady=(0, 10))
            if self.typevar.get() =="normal":
                self.run_progress['value'] = self.obj.steps_calculation
            else:
                self.run_progress.config(mode="indeterminate")
                self.run_progress.start()
            self.simulation_complete_event = threading.Event()
            self.thread = threading.Thread(target=self.run_simulation, daemon=True, args=[run_params])
            self.thread.start()
            self.after(1000, self.get_simulation_progress)

    def on_check_plot_method(self, event):
        if event == "drifter_properties":
            self.p1_label.config(state="normal")
            self.p1.config(state="normal")
            self.p2_label.config(state="normal")
            self.p2.config(state="normal")
            self.loc_label.config(state="normal")
            self.loc_entry.config(state="normal")
        else:
            self.p1_label.config(state="disabled")
            self.p1.config(state="disabled")
            self.p2_label.config(state="disabled")
            self.p2.config(state="disabled")
            self.loc_label.config(state="disabled")
            self.loc_entry.config(state="disabled")
    
    def on_click_plot(self):
        self.update_idletasks()
        self.launch_progress_window("Plotting")
        self.thread = threading.Thread(target=self.plot_simulation, daemon=True)
        self.thread.start()
        self.after(1000, self.get_plot_progress)
    
    def plot_simulation(self):
        try:
            title = None if self.title_entry.get() == "" else self.title_entry.get()
            xlabel = None if self.xlabel_entry.get() == "" else self.xlabel_entry.get()
            ylabel = None if self.ylabel_entry.get() == "" else self.ylabel_entry.get()
            cblabel = None if self.cblabel_entry.get() == "" else self.cblabel_entry.get()
            outfile = self.oplot_entry
            depth = -200 if self.depth_entry.get() == "" else -abs(float(self.depth_entry.get()))
            if outfile is not None:
                plots.mpl.use("Agg")
            else:
                plots.mpl.use('TkAgg')
            if self.loc_label["state"] == "disabled":
                locs = None
            else:
                locs = configure_locations(self.loc_entry.get())
            
            p = plots.Plot(*np.append(self.plotfilesplus, self.plotfilesminus), add = self.ADD.get(), diff=self.DIFF.get(),
                        diffidx=len(self.plotfilesplus), title = title, xlabel=xlabel,
                        ylabel=ylabel, colorbarlabel=cblabel, absolute=self.ABS.get(),
                        locations=locs, prop1=self.p1var.get(), prop2=self.p2var.get(),
                        depth=depth, outfile = outfile)
            getattr(p, self.methodvar.get())()
            if outfile is None:#TODO: Handle mainthread is not in main loop error.
                plt.close("all")
                display_window = tkinter.Toplevel(self)
                canvas = FigureCanvasTkAgg(p.fig, master=display_window)
                canvas.draw()
                canvas.get_tk_widget().pack()
                plt.close("all")

            self.progress_window.destroy()
            return 
        except Exception as e:
            self.progress_window.destroy()
            self.launch_error_window(e)
            return
    
    def get_plot_progress(self):
        if self.thread.is_alive():
            self.after(1000, self.get_plot_progress)
        else:
            if self.oplot_entry is not None:
                self.viewbutton.config(state="normal")
    
    def on_click_view_plot(self):
        filepath = self.oplot_entry
        format = filepath[-4:]
        if format[-4:] == ".gif":
            func = self.view_ANIM
        elif format in [".pdf", ".png", ".jpg"]:
            if format == ".pdf":
                try:
                    from pdf2image import convert_from_path
                except ImportError:
                    raise ImportError("Format pdf cannot be imported.\
                                       Please install pdf2image or convert the input to a different format.")
            func = self.view_IMG
        else:
            raise TypeError(f"Format {format} is not supported in this GUI.")
        self.launch_progress_window(f"Importing {filepath}")
        self.thread = threading.Thread(target=func, daemon=True, args=[filepath])
        self.thread.start()
    
    def view_ANIM(self, filepath):
        info = Image.open(filepath)
        frames = info.n_frames
        images = []
        size = self.get_display_geom(*info.size)
        
        for i in range(frames):
            info.seek(i)
            frame = info.copy().convert("RGBA")
            frame = frame.resize(size, Image.LANCZOS)

            tk_frame = ImageTk.PhotoImage(frame)
            images.append(tk_frame)
        display_window = tkinter.Toplevel(self)
        display_window.title(filepath)
        display_window.geometry(f"{size[0]}x{size[1]}")
        gif_label = tkinter.Label(display_window, image="")
        display_window.bind()
        self.progress_window.destroy()
        gif_label.pack()
        self.run_animation(images, gif_label, display_window)
    
    def run_animation(self, images, gif_label, display_window, current_frame=0):
        if not display_window.winfo_exists():
            return
        image = images[current_frame]

        gif_label.configure(image = image)
        current_frame = current_frame + 1

        if current_frame == len(images):
            current_frame = 0 # reset the current_frame to 0 when end is reached

        self.after(50, lambda: self.run_animation(images, gif_label, display_window, current_frame))
    
    def view_IMG(self, filepath):
        if filepath[-4:] == ".pdf":
            from pdf2image import convert_from_path

            image = convert_from_path(filepath)[0]
            image = image.convert("RGB")
        else:
            image = Image.open(filepath)
        size = self.get_display_geom(*image.size)
        image.thumbnail(size)
        image = ImageTk.PhotoImage(image)

        display_window = tkinter.Toplevel(self)
        display_window.title(filepath)
        display_window.geometry(f"{size[0]}x{size[1]}")
        image_label = tkinter.Label(display_window, image=image)
        image_label.image = image
        self.progress_window.destroy()
        image_label.pack(pady=5)

    def get_display_geom(self, width, height):
        if width <= self.SIMDISPLAYMAXWIDTH and height <= self.SIMDISPLAYMAXHEIGHT:
            return (width, height)
        else:
            if width > height:
                scale = width/self.SIMDISPLAYMAXWIDTH
            else:
                scale = height / self.SIMDISPLAYMAXHEIGHT
            return (round(width / scale), round(height / scale))
    
    def on_click_diff(self):
        if self.DIFF.get():
            self.simulation_label_minus.config(state = "normal")
            self.simulation_explore_minus.config(state = "normal")
        else:
            self.simulation_label_minus.config(state = "disabled")
            self.simulation_explore_minus.config(state = "disabled")

    def browse_seed(self):
        filetypes = (
            ('pickle files', '*.pkl'),
            ('All files', '*.*')
        )

        filename = filedialog.askopenfilename(initialdir = os.getcwd(),
                                          title = "Select a File",
                                          filetypes = filetypes)
        self.seed_file_entry = filename
        self.set_frame_state(self.seed_new, "disabled")
        self.seed_file_label.config(text=f"Imported from {filename}")
        self.remove_seed_button.config(state = "normal")
    
    def remove_seed(self):
        self.set_frame_state(self.seed_new, "normal")
        self.seed_file_entry = ""
        self.seed_file_label.config(text="Import seed:")
        self.remove_seed_button.config(state="disabled")

    def browse_required_readers(self, reader):
        filetypes = (
            ('netCDF files', '*.nc'),
            ('All files', '*.*')
        )

        filename = filedialog.askopenfilename(initialdir = self.supplementary_data_dir,
                                          title = "Select a File",
                                          filetypes = filetypes)
        if reader == "tmp":
            self.tmp_entry = filename
            self.tmp_label.config(text=f"Imported from {filename}")
        else:
            self.bmt_entry = filename
            self.bmt_label.config(text=f"Imported from {filename}")

    def browse_optional_readers(self):
        filetypes = (
            ('netCDF files', '*.nc'),
            ('All files', '*.*')
        )

        filenames = filedialog.askopenfilenames(initialdir = os.getcwd(),
                                          title = "Select files",
                                          filetypes = filetypes)
        self.opt_readers = filenames
        self.opt_label.config(text = f"Imported {len(filenames)} from {os.path.dirname(filenames[0])}.")

    def browse_simulations(self, sign):
        filetypes = (
            ('netCDF files', '*.nc'),
            ('All files', '*.*')
        )

        filenames = filedialog.askopenfilenames(initialdir = os.getcwd(),
                                          title = "Select a File",
                                          filetypes = filetypes)
        if sign=="plus":
            self.plotfilesplus = list(filenames)
            self.simulation_label_plus.config(text=f"Uploaded {len(filenames)}\nsimulations")
        else:
            self.plotfilesminus = list(filenames)
            self.simulation_label_minus.config(text=f"Uploaded {len(filenames)}\nsimulations")

    def get_simulation_progress(self):
        if self.typevar.get() =="normal":
            num_steps = int(self.step_entry.get())
            current_step = self.obj.steps_calculation
            progress_value = current_step / num_steps * 100
            self.run_progress['value'] = round(progress_value)
        
        if self.thread.is_alive():
            self.after(1000, self.get_simulation_progress)
        else:
            self.run_progress.destroy()

    def run_simulation(self, args):
        self.obj.run(**args)

    def save_ev_selection(self, event):
        """Store selected items when listbox loses focus."""
        self.saved_ev_selection = set(self.ev_list.curselection())
    
    def restore_ev_selection(self, event):
        """Restore selection when listbox regains focus."""
        for index in self.saved_ev_selection:
            self.ev_list.selection_set(index)

    def set_frame_state(self, frame, state):
        for widget in frame.winfo_children():
            widget.config(state=state)
        if frame==self.seed_new:
            if state=="disabled":
                pass
            else:
                self.on_check_relative_mass()
    
    def launch_error_window(self, error):
        self.error_window = tkinter.Toplevel(self)
        self.error_window.title(type(error).__name__)
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
        self.after(100, self.check_seed_progress)
    
    def check_seed_progress(self):
        if self.thread.is_alive():
            self.after(100, self.check_seed_progress)
        else:
            if type(self.sdata) == Seed:
                save = tkinter.messagebox.askyesno("Save seed?", "Would you like to save the seed for future use?")
                if save:
                    self.save_seed()
            self.freeze_seed_frame()
    
    def on_click_reset_seed(self):
        self.SEEDCLICKED = False
        self.sdata = None
        self.seed_file_entry=""
        self.seed_file_label.config(text="Import seed:")
        self.set_frame_state(self.seed_new, "normal")
        self.set_frame_state(self.seed_from_file, "normal")
        self.reset_seed_button.config(state="disabled")

    def freeze_seed_frame(self):
        self.SEEDCLICKED = True
        self.set_frame_state(self.seed_new, "disabled")
        self.set_frame_state(self.seed_from_file, "disabled")
        self.reset_seed_button.config(state = "normal")

    def save_seed(self):
        filetypes = (
            ('pickle files', '*.pkl'),
            ('All files', '*.*')
        )

        file = filedialog.asksaveasfile(initialdir = os.getcwd(),
                                          title = "Select a location",
                                          filetypes = filetypes, defaultextension = filetypes)
        
        if file is None: return
        try:
            self.sdata.save_seed(file.name)
        except Exception as e:
            self.launch_error_window(e)
    
    def seed(self):
        entry = self.seed_file_entry
        if len(entry) > 0:
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
                                  areapath=data["area"], biomegridpath=data["biome"], initialmassdata=data["initial_mass"])
            except Exception as e:
                self.launch_error_window(e)

        self.progress_window.destroy()

    def add_readers(self):
        readers = []
        if self.bmt_entry is None:
            raise MissingDataError("Bathymetry reader must be provided.")
        bat = reader_netCDF_CF_generic.Reader(self.bmt_entry, standard_name_mapping={"topo":"depth"})
        readers.append(bat)
        if self.tmp_entry is None:
            raise MissingDataError("Temperature reader must be provided.")
        tmp = reader_netCDF_CF_generic.Reader(self.tmp_entry)
        tmp.verticalbuffer = 100
        tmp.always_valid = True
        readers.append(tmp)
        if self.LANDMASK.get() and not self.AUTOLANDMASK.get():
            readers.append(reader_global_landmask.Reader())
        opt_readers = self.opt_readers
        for reader in opt_readers:
            readers.append(reader_netCDF_CF_generic.Reader(reader))
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
        ffunc = None if not self.FRAGMENTATION.get() else valid_fragmentation_function(self.ffunc_entry.get())
        params = {"m0":self.sdata.mass, "loglevel":0, "decay_type":self.decayvar.get(), "initial_velocity":-abs(float(self.w0_entry.get())),
                "vertical_velocity_type": self.w0var.get(), "ffunction":ffunc}
        
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
        elif hasattr(self.sdata, "biome"):
            origin_marker = self.sdata.biome
        else:
            origin_marker = 0
        
        self.obj.seed_elements(lon=self.sdata.lon, lat = self.sdata.lat, z=self.sdata.z,
                                   mass = self.sdata.mass, time=valid_date(self.start_entry.get()), origin_marker = origin_marker)
    
    def initialize_grid_run(self):
        if hasattr(self.sdata, "origin_marker"):
            origin_marker = self.sdata.origin_marker
        elif hasattr(self.sdata, "biome"):
            origin_marker = self.sdata.biome
        else:
            origin_marker = 0
        
        sf = 1 if self.sf_entry.get()=="" else int(self.sf_entry.get())
        ffunc = None if not self.FRAGMENTATION.get() else valid_fragmentation_function(self.ffunc_entry.get())
        params = {"m0":self.sdata.mass, "lon":self.sdata.lon, "lat": self.sdata.lat, "z": self.sdata.z,
                    "loglevel":0, "decay_type":self.decayvar.get(), "initial_velocity":-abs(float(self.w0_entry.get())),
                    "vertical_velocity_type": self.w0var.get(), "ffunction":ffunc, "split_factor":sf,
                    "reader": self.add_readers(), "configure":self.add_configures(), "origin_marker":origin_marker,
                    "deactivate_horizontal_advection": not self.ADVECTION.get(), "deactivate_fragmentation": not self.FRAGMENTATION.get(),
                    "starttime":valid_date(self.start_entry.get()), }
        
        if self.mdecayvar.get() == "mass":
            self.obj = mcdgrid.GridRun(**params)
        else:
            self.obj = acdgrid.GridRun(**params)
    
    def reset_object_init(self):
        self.obj = None
    
    def save_simulation(self):
        filetypes = (
            ('netCDF files', '*.nc'),
            ('All files', '*.*')
        )

        file = filedialog.asksaveasfile(initialdir = os.getcwd(),
                                          title = "Select a location",
                                          filetypes = filetypes, defaultextension = filetypes)
        
        if file is None: 
            self.simout_entry = None
        else:
            self.simout_entry = file.name
            os.remove(file.name)
            self.out_label.config(text=f"Saved to: {file.name}")
    
    def save_plot(self):
        filetypes = (
            ('png files', '*.png'),
            ('jpg files', '*.jpg'),
            ('pdf files', '*.pdf'),
            ('gif files', '*.gif'),
            ('All files', '*.*')
        )

        file = filedialog.asksaveasfile(initialdir = os.getcwd(),
                                          title = "Select a location",
                                          filetypes = filetypes, defaultextension = filetypes)
        
        if file is None: 
            self.oplot_entry = None
        else:
            self.oplot_entry = file.name
            os.remove(file.name)
            self.outfile_plot_label.config(text=f"Saved to: {file.name}")

def main():
    GUI = CarbonDriftGUI()
    GUI.mainloop()

if __name__ == "__main__":
    main()