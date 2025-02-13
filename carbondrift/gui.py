import numpy as np
import matplotlib.pyplot as plt
from PIL import ImageTk, Image
import tkinter
from tkinter import ttk
# from importlib import resources

# from datetime import datetime, timedelta
# from copy import copy
import carbondrift
# import carbondrift.models.massdecay.gridrun as mcdgrid
# import carbondrift.models.areadecay.gridrun as acdgrid
# import carbondrift.models.massdecay.carbondrift as mcd
# import carbondrift.models.areadecay.carbondrift as acd
# from carbondrift.simulation.param_classifier import Parameters
# from opendrift.readers import reader_global_landmask
# from opendrift.readers import reader_netCDF_CF_generic
# from carbondrift.simulation.seeding import *

from carbondrift.models.logger import Logger

log = Logger("carbondrift.carbondrift_gui")
logger = log.LOGGER

class CarbonDriftGUI(tkinter.Tk):
    def __init__(self):
        tkinter.Tk.__init__(self)
        self.title("CarbonDrift " + carbondrift.__version__ + " GUI")
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

        frame_subtitle_font = ("Courier", 18, "bold")
        frame_font = ("Helvetica", 12)

        self.seed_from_file = tkinter.Frame(self.simseed, bg='lightgray', bd=2,
                            relief=tkinter.SUNKEN, pady=25, padx=25)
        self.seed_from_file.grid(row=0, column=0, rowspan=1, padx=20, pady=(20,0), sticky="we")
        seed_from_file_title = tkinter.Label(self.seed_from_file, text="Seed from pre-made file", 
                                             font=frame_subtitle_font).grid(row=0, pady=(0, 25), columnspan=2)
        seed_file_label = tkinter.Label(self.seed_from_file, text="File Path:", font=frame_font).grid(row=1, column=0, sticky="w")
        seed_file_entry = tkinter.Entry(self.seed_from_file).grid(row=1, column=1, sticky="w")

        self.seed_new = tkinter.Frame(self.simseed, bg='lightgray', bd=2,
                            relief=tkinter.SUNKEN, pady=25, padx=25)
        self.seed_new.grid(row=1, column=0, rowspan=1, padx=20, pady=(20,0), sticky="we")
        seed_from_file_title = tkinter.Label(self.seed_new, text="Make New Seed", 
                                             font=frame_subtitle_font).grid(row=0, pady=(0, 25), columnspan=2)
        relative_mass = tkinter.IntVar()
        relative_mass_button = tkinter.Checkbutton(self.seed_new, text="relative mass",
                                                            variable=relative_mass).grid(row=1, column=0, sticky="w")
        save_seed = tkinter.IntVar()
        save_seed_button = tkinter.Checkbutton(self.seed_new, text="save seed",
                                                            variable=save_seed).grid(row=1, column=1, sticky="w")
        lonmin_label = tkinter.Label(self.seed_new, text="lonmin", font=frame_font).grid(row=2, column=0)
        lonmin_entry = tkinter.Entry(self.seed_new).grid(row=3, column=0)
        lonmax_label = tkinter.Label(self.seed_new, text="lonmax", font=frame_font).grid(row=2, column=1)
        lonmax_entry = tkinter.Entry(self.seed_new).grid(row=3, column=1)
        latmin_label = tkinter.Label(self.seed_new, text="latmin", font=frame_font).grid(row=2, column=2)
        latmin_entry = tkinter.Entry(self.seed_new).grid(row=3, column=2)
        latmax_label = tkinter.Label(self.seed_new, text="latmax", font=frame_font).grid(row=2, column=3)
        latmax_entry = tkinter.Entry(self.seed_new).grid(row=3, column=3)
        dx_label = tkinter.Label(self.seed_new, text="dx", font=frame_font).grid(row=2, column=4)
        dx_entry = tkinter.Entry(self.seed_new).grid(row=3, column=4)
        dy_label = tkinter.Label(self.seed_new, text="dy", font=frame_font).grid(row=2, column=5)
        dy_entry = tkinter.Entry(self.seed_new).grid(row=3, column=5)
        bathymetry_label = tkinter.Label(self.seed_new, text="bathymetry", font=frame_font).grid(row=4, column=0)
        bathymetry_entry = tkinter.Entry(self.seed_new).grid(row=5, column=0)
        m0_label = tkinter.Label(self.seed_new, text="initial mass", font=frame_font).grid(row=4, column=1)
        m0_entry = tkinter.Entry(self.seed_new).grid(row=5, column=1)
        area_label = tkinter.Label(self.seed_new, text="area", font=frame_font).grid(row=4, column=2)
        area_entry = tkinter.Entry(self.seed_new).grid(row=5, column=2)
        biome_label = tkinter.Label(self.seed_new, text="biome", font=frame_font).grid(row=4, column=3)
        biome_entry = tkinter.Entry(self.seed_new).grid(row=5, column=3)
        outfile_label = tkinter.Label(self.seed_new, text="outfile", font=frame_font).grid(row=4, column=4)
        outfile_entry = tkinter.Entry(self.seed_new).grid(row=5, column=4)
        
        self.seed_start = tkinter.Frame(self.simseed, relief=tkinter.FLAT, pady=25, padx=25)
        self.seed_start.grid(row=2, column=0, rowspan=1)
        tkinter.Button(self.seed_start, text="SEED", bg='green').grid(row=0)


if __name__ == "__main__":
    GUI = CarbonDriftGUI()
    GUI.mainloop()