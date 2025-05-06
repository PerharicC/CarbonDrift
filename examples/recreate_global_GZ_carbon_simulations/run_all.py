import subprocess
import os
import numpy as np
import itertools
from datetime import timedelta
import carbondrift


curdir = os.getcwd()

os.chdir("/home/peharicc/CarbonDrift")
sup_dir = f"{os.path.dirname(os.path.dirname(carbondrift.__file__))}/supplementary_data"
exmpl_dir = f"{os.path.dirname(os.path.dirname(carbondrift.__file__))}/examples/recreate_global_GZ_carbon_simulations"
script = "run_simulation"

outfiles_path = exmpl_dir
parameters = f"@{exmpl_dir}/params_all.txt"

velocity_types = ["constant", "variable"]
phylum = ["cnidaria", "ctenophora", "chordata"]
decay_types = ["exp"]#["linear", "exp"]
seeddata = ["Eg", "M"]
variable_params = [velocity_types, phylum, decay_types, seeddata]
combinations = list(itertools.product(*variable_params))

def find_combination(wtype, w, decay, z):
    return combinations.index((wtype, w, decay, z))

def run():
    k = 0
    for wtype, p, decay, poctype  in combinations:
        folder = decay + "_" + wtype
        outfile = os.path.join(outfiles_path, f"{folder}/{p}_{poctype}.nc")
        sdata = os.path.join(sup_dir, f"{p}_{poctype}_seed.pkl")
        if poctype == "Eg":
            if p == "chordata":
                w0 = -700 / 24 / 3600
                stoutput = "1:0:0"
            else:
                w0 = -100 / 24 / 3600
                stoutput = "4:0:0"
        else:
            if p == "cnidaria":
                w0 = -1100 / 24 / 3600
                stoutput = "1:0:0"
            elif p == "ctenophora":
                w0 = -900 / 24 / 3600
                stoutput = "1:0:0"
            else:
                w0 = -1000 / 24 / 3600
                stoutput = "1:0:0"
        steps = int(round(10000 / abs(w0) / 3600 / 0.5))+100
        tmp = os.path.join(sup_dir, "tmp_luo_remaped1deg.nc")
        btm = os.path.join(sup_dir, "etopo2_remaped1deg.nc")
        params = f"{parameters} -st {steps} -o {outfile} -sdata {sdata} -wtype {wtype} --decaytype {decay} -w0 {w0} -dto {stoutput} -tmp {tmp} -b {btm}"
        command = f"{script} {params}"
        print(f"Running {k} out of {len(combinations)}")
        subprocess.run(command, shell = True)
        
        k += 1

run()
