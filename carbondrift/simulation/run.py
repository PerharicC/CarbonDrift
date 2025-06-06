from datetime import datetime, timedelta
import argparse
import numpy as np #This must be here for proper evaluation of the fragmentation function!!!
from carbondrift.simulation.prepare_run import PrepareSimulation

from carbondrift.models.logger import Logger

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
    new_path = path.split(".")[0]+".json"
    with open(new_path, "w") as file:
        json.dump(p, file, indent=4, default=str)

def valid_export_variables(s):
    valid_vars = ["trajectory", "time", "status", "moving", "age_seconds", "origin_marker",
                  "lon", "lat", "z", "wind_drift_factor", "current_drift_factor",
                  "terminal_velocity", "mass", "x_sea_water_velocity", "y_sea_water_velocity",
                  "land_binary_mask", "sea_water_temperature", "depth"]
    
    try:
        level = int(s)
        if level == 0:
            return valid_vars
        elif level == 1:
            return ["trajectory", "time", "status", "origin_marker", "lon", "lat", "z", "mass",
                    "x_sea_water_velocity", "y_sea_water_velocity", "sea_water_temperature", "depth"]
        elif level == 2:
            return ["trajectory", "time", "status", "origin_marker", "lon", "lat", "z", "mass", "sea_water_temperature", "depth"]
        elif level >= 3:
            return ["trajectory", "time", "status", "origin_marker", "lon", "lat", "z", "mass", "sea_water_temperature"]
    except ValueError:
        pass

    try:
        variables = s.split(",")
    except Exception:
        raise argparse.ArgumentTypeError(f"not a valid export variables format: {s!r}")
    
    for  var in variables:
        if var not in valid_vars:
            raise argparse.ArgumentError(f"Export variable {var} is not a valid CarbonDrift variable.")
    return variables


def prepare_fragmentation(function_code):
    def get_variable_from_code(s):
        var = ""
        for k in s:
            if k == ";":
                break
            var+=k
        return var
    
    modified_code = ""
    dictionary_of_variables = {
        "T":"environment.sea_water_temperature",
        "m": "elements.mass",
        "u": "environment.x_sea_water_velocity",
        "v": "environment.y_sea_water_velocity",
        "uh": "np.sqrt(environment.x_sea_water_velocity**2 + environment.y_sea_water_velocity ** 2)",
        "z": "elements.z",
        "m0": "m0"
        }
    variable_keys = dictionary_of_variables.keys()
    found = False

    i = 0
    while i < len(function_code):
        s = function_code[i]
        if s != ";":
            modified_code += s
            i+=1
            continue
        var = get_variable_from_code(function_code[i+1:])
        
        for key in variable_keys:
            if key == var:
                found = True
                break
        if not found:
            raise NameError(f"Variable {var} not in variable list {variable_keys}.")
        
        modified_code += dictionary_of_variables[key]
        found = False
        i+= (2 + len(var))
    return lambda elements, environment, m0: eval(modified_code)

def valid_fragmentation_function(s):
    try:
        return prepare_fragmentation(s)
    except Exception as e:
        print(f"When trying to set-up the fragmentation function, the following error was raised: {type(e).__name__}")
        print(f"Error message: {str(e)}")
        raise


def main():
    p = argparse.ArgumentParser(fromfile_prefix_chars='@')
    p.add_argument("-sim", "--simulation_type", default="normal", choices=["normal", "grid"],
                   help = "Type of simulation, recommended normal unless you use fragmentation.")
    p.add_argument("-tmp", "--temperature", help="NetCDF temperature file.", type = str, required=True)
    p.add_argument("-c", "--current", help="NetCDF current file.", type = str, required=False)
    p.add_argument("-o", "--outfile", type = str, default = "outfile.nc", help = "NetCDF outfile path.")
    p.add_argument("-b", "--bathymetry", type = str, help = "NetCDF depth file.", required=True)
    p.add_argument("-buff", "--readerbuffer", type = int, default=100, help ="OceanDrift reader vertical buffer size.\
                   Only change if it is to small.")
    p.add_argument("-s", "--starttime", type = valid_date, required  = False,
                   help = "Start time of simulation - format YYYY-MM-DD-h")
    p.add_argument("-w0", "--initialvelocity", default = -0.01, type = float,
                   help = "Initial velocity in m/s, positive up. Deafult is -0.01 m/s.")
    p.add_argument("-ffunc", "--ffunction", required = False, type=valid_fragmentation_function,
                   help = "Fragmentation function, variables must be enclosed in ;; - format np.sqrt(;T;) / ;u; ** 2>np.random(len(;T;))")
    p.add_argument("-f", "--fragmentation", action="store_false", help = "Set if fragmentation should be allowed.")
    p.add_argument("-a", "--advection", action="store_false", help = "Set if horizontal advection should be allowed.")
    p.add_argument("-st", "--steps", default=200, type = int, help = "Simulation steps.")
    p.add_argument("-dt", "--timestep",type = valid_time_step, default=timedelta(minutes = 30),
                   help = "Time step of simulation - format H:M:S" )
    p.add_argument("-dto", "--timestepoutput",type = valid_time_step, default=timedelta(minutes = 30),
                   help = "Output time step of simulation - format H:M:S" )
    p.add_argument("-d", "--decaytype", type = str, default="linear", choices=["linear", "exp"])
    p.add_argument("-log", "--loglevel", default = 0, type = int)
    p.add_argument("-oceanonly", action="store_true")
    p.add_argument("-sdata", "--seeddata", type = str, help = "Path to initial seed data.")
    p.add_argument("-mdtype", "--microbialdecaytype", type = str, default="mass", choices = ["mass", "area"])
    p.add_argument("-wtype", "--vertical_velocity_type", type = str, choices = ["constant", "variable"], default = "variable")
    p.add_argument("-ev", "--export_variables", type = valid_export_variables, default=valid_export_variables("0"),
    help = "Export variables. Either separate variables with , or specify export level with number (0 for all variables (default)).")
    p.add_argument("-sf", "--splitfactor", type = int, default=1,
                   help="Split simulation to n subsimulations using grid simulation_type.")
    args = p.parse_args()
    run_simulation(**vars(args))

def run_simulation(**kwargs):
    log = Logger("CarbonDrift.simulation.run")
    logger = log.LOGGER

    logger.info("Saving parameters.")
    save_params(kwargs["outfile"], kwargs)

    logger.info("Preparing Simulation.")
    simulation = PrepareSimulation(**kwargs)
    logger.info("Simulation is ready.")
    o = simulation.obj
    logger.info("Starting simulation.")
    o.run(**simulation.p.simulation_run)
    logger.info("Simulation is finished.")



if __name__ == '__main__':
    main()