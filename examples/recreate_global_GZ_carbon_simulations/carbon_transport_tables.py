import pandas as pd
import numpy as np
import subprocess
import pickle
import os

curdir = os.getcwd()

def save_masses(workdir, massfile):
    script = "plot_simulation"
    params = "get_biome_weighted_mass_at_depth"
    mass = {}
    for phylum in ["cnidaria", "ctenophora", "chordata"]:
        for poc_type in ["M", "Eg"]:
            mass[phylum + "_"+poc_type] = []
            for depth in [-100, -1000, -6000]:
                file = os.path.join(workdir, phylum + "_" + poc_type + ".nc")
                command = f"{script} {params} {file} -abs -d {depth}"
                result = subprocess.run(command, shell = True, capture_output=True, text=True)
                mass[phylum + "_"+poc_type].append([float(i) for i in result.stderr.split("\n")[-2].strip("[]").split(" ")])
    with open(massfile, "wb") as f:
        pickle.dump(mass, f)

def load_masses(massfile):
    with open(massfile, "rb") as f:
        data = pickle.load(f)
    return data

def run(decaytype):
    biomes = {"COAST":3, "HCPS":2, "HCSS":0, "LC":1}
    for decay in decaytype:
        workdir = os.path.join(curdir, decay)
        massfile = os.path.join(curdir, f"images/masses_{decay}.pkl")
        save_masses(workdir, massfile)
        mass = load_masses(massfile)
        Cnidaria_Eg = mass["cnidaria_Eg"]
        Cnidaria_M = mass["cnidaria_M"]
        Ctenophora_Eg = mass["ctenophora_Eg"]
        Ctenophora_M = mass["ctenophora_M"]
        Chordata_M = mass["chordata_M"]
        Chordata_Eg = mass["chordata_Eg"]

        Cnidaria = np.asarray(Cnidaria_Eg) + np.asarray(Cnidaria_M)
        Ctenophora = np.asarray(Ctenophora_Eg) + np.asarray(Ctenophora_M)
        Chordata = np.asarray(Chordata_Eg) + np.asarray(Chordata_M)
        Total = Cnidaria + Ctenophora + Chordata
        Total /= 10 ** 15
        rows = []
        values = np.zeros((20, 3))
        i = 0

        for phylum, fluxes in zip(["Cnidaria", "Ctenophora", "Chordata"], [Cnidaria, Ctenophora, Chordata]):
            fluxes/=10 ** 15
            for biome in biomes.keys():
                rows.append((phylum, biome))
                val = fluxes[:, biomes[biome]]
                values[i] = [round(val[i], 2) for i in range(len(val))]
                i+=1
            rows.append((phylum, "SUM"))
            values[i] = [round(x, 2) for x in np.sum(fluxes, axis=1)]
            i+=1

        for biome in biomes.keys():
            rows.append(("All GZ", biome))
            val = Total[:, biomes[biome]]
            values[i] = [round(val[i], 2) for i in range(len(val))]
            i+=1
        rows.append(("All GZ", "SUM"))
        values[i] = [round(x, 2) for x in np.sum(Total, axis=1)]

        columns = ["Export past 100 m [Pg/Y]", "Export past 1000 m [Pg/Y]", "Flux to Seafloor [Pg/Y]"]
        ind = pd.MultiIndex.from_tuples(rows, names = ["Phylum", "Biome"])
        data = pd.DataFrame(values, index = ind, columns= columns)
        outfile = f"images/{decay}.tex"
        data.to_latex(outfile, float_format="%.2f", index=False)

if __name__ == "__main__":
    run(["exp_constant", "exp_variable"])