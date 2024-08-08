from datetime import datetime, timedelta
import argparse

from simulation.prepare_run import PrepareSimulation

from model.logger import Logger

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

def main():
    p = argparse.ArgumentParser()
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
    p.add_argument("-lon", required=False, type = str, help = "Longitude seeding - format min-max-seeednum.")
    p.add_argument("-lat", required=False, type = str, help = "Latitude seeding - format min-max-seeednum.")
    p.add_argument("-oceanonly", action="store_true")
    args = p.parse_args()
    run_simulation(**vars(args))

def run_simulation(**kwargs):
    log = Logger("CarbonDrift.scripts.run")
    logger = log.LOGGER
    logger.info("Preparing Simulation.")
    simulation = PrepareSimulation(**kwargs)
    logger.info("Simulation is ready.")
    o = simulation.obj
    logger.info("Starting simulation.")
    o.run(**simulation.p.simulation_run)
    logger.info("Simulation is finished.")



if __name__ == '__main__':
    main()