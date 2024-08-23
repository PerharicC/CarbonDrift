class Parameters:

    def __init__(self, params):
        
        self.readers = {}
        self.object_init = {}
        self.simulation_run = {}
        self.seeding = {"seedtype": params.pop("seed")}
        self.classify(params)
    
    def classify(self, params):

        for key, value in params.items():
            if key in ["temperature", "current", "bathymetry", "temperaturebuffer"]:
                self.readers[key] = value
                if key == "temperature":
                    self.seeding["ncfile"] = value
            elif key == "decaytype":
                self.object_init["decay_type"] = value
            elif key == "initialvelocity":
                self.object_init["initial_velocity"] = value
            elif key == "advection":
                self.object_init["deactivate_horizontal_advection"] = value
            elif key == "fragmentation":
                self.object_init["deactivate_fragmentation"] = value
            elif key in ["loglevel", "starttime", "distribution"]:
                self.object_init[key] = value
            elif key == "timestep":
                self.simulation_run["time_step"] = value
            elif key in ["steps", "outfile"]:
                self.simulation_run[key] = value
            elif key in ["lon", "lat", "skip", "massdata", "species", "z"]:
                self.seeding[key] = value
            else:
                raise Warning("Got an unexpected parameter: " + key)
            
