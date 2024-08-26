from datetime import datetime, timedelta
import argparse

from simulation.prepare_run import PrepareSimulation

from model.logger import Logger

import json

def valid_date(s):
    try:
        return datetime.strptime(s, "%Y-%m-%d-%H")
    except ValueError:
        raise argparse.ArgumentTypeError(f"not a valid date: {s!r}")

def valid_time_step(s):
    try:
        t = datetime.strptime(s,"%H:%M:%S")
        return timedelta(hours = t.hour, minutes=t.minute, seconds = t.second)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not a valid time format: {s!r}")

def save_params(path, params):
    import time
    from copy import copy
    p = copy(params)
    p["TIME OF SIMULATION RUN"] = time.asctime()
    
    with open(path, "a") as file:
        json.dump(p, file, indent=4, default=str)

def main():
    p = argparse.ArgumentParser(fromfile_prefix_chars='@')
    p.add_argument("-sim", "--simulation_type", default="normal", choices=["normal", "grid", "DAC"])
    p.add_argument("-tmp", "--temperature", help="NetCDF temperature file.", type = str, required=True)
    p.add_argument("-c", "--current", help="NetCDF current file.", type = str, required=False)
    p.add_argument("-o", "--outfile", type = str, default = "outfile.nc")
    p.add_argument("-b", "--bathymetry", type = str)
    p.add_argument("-buff", "--temperaturebuffer", type = int, default=100)
    p.add_argument("-s", "--starttime", type = valid_date, required  = False,
                   help = "Start time of simulation - format YYYY-MM-DD-h")
    p.add_argument("-w0", "--initialvelocity", default = -0.01, type = float)
    p.add_argument("-dist", "--distribution", required = False, default = "mass_sqrt")
    p.add_argument("-f", "--fragmentation", action="store_false")
    p.add_argument("-a", "--advection", action="store_false")
    p.add_argument("-st", "--steps", default=200, type = int)
    p.add_argument("-dt", "--timestep",type = valid_time_step, default=timedelta(minutes = 30),
                   help = "Time step of simulation - format H:M:S" )
    p.add_argument("-d", "--decaytype", type = str, default="linear", choices=["linear", "exp"])
    p.add_argument("-log", "--loglevel", default = 0, type = int)
    p.add_argument("-skip", default = 1, type = int)
    p.add_argument("-seed", type = str, choices = ["line", "rectangle"], default = "rectangle")
    p.add_argument("-lon", required=False, type = str, help = "Longitude seeding - format min:max:seeednum. Use m for minus sign.")
    p.add_argument("-lat", required=False, type = str, help = "Latitude seeding - format min:max:seeednum. Use m for minus sign.")
    p.add_argument("-oceanonly", action="store_true")
    p.add_argument("-species", type = str, default = "Chordata")
    p.add_argument("-massdata", type = str, help = "Path to initaal mass data. If None, mass is set to one in all cells.")
    p.add_argument("-z", type = int, help = "Seeding depth.", default=0)
    p.add_argument("-mgtype", "--massgen_type", type = str, choices = ["from_file", "constant", "from_radius"])
    p.add_argument("-mdtype", "--microbialdecaytype", type = str, default="mass", choices = ["mass", "area"])
    args = p.parse_args()
    run_simulation(**vars(args))

def run_simulation(**kwargs):
    log = Logger("CarbonDrift.simulation.run")
    logger = log.LOGGER

    logger.info("Saving parameters.")
    save_params("saved_params.json", kwargs)

    logger.info("Preparing Simulation.")
    simulation = PrepareSimulation(**kwargs)
    logger.info("Simulation is ready.")
    o = simulation.obj
    logger.info("Starting simulation.")
    o.run(**simulation.p.simulation_run)
    logger.info("Simulation is finished.")



if __name__ == '__main__':
    main()