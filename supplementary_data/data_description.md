# Supplementary Data Contents

Below we give a short description of the contents in this directory, together with appropriate references to the datasets.

## Luo_M_Eg_biome_data.json

A `JSON` file containing the **total GZ carbon content** that is ready for vertical export per each biome. The values are taken from Ref. [[1](#references)].

```json
{
    "cnidaria": {
        "M": {
            "COAST": 25000000000000,
            "HCPS": 243000000000000,
            "HCSS": 272000000000000,
            "LC": 136000000000000
        },
        "Eg": {
            "COAST": 19300000000000,
            "HCPS": 176000000000000, 
            "HCSS": 182000000000000,
            "LC": 108000000000000
        },
        "units": "gY-1"
    },
    "ctenophora": {
        "M": {
            "COAST": 1980000000000, 
            "HCPS": 6370000000000,
            "HCSS": 14400000000000, 
            "LC": 6680000000000
        },
        "Eg": {
            "COAST": 21900000000000, 
            "HCPS": 76700000000000, 
            "HCSS": 182000000000000, 
            "LC": 85100000000000
        },
        "units": "gY-1"
    },
    "chordata": {
        "M": {
            "COAST": 25600000000000,
            "HCPS": 181000000000000,
            "HCSS": 144000000000000,
            "LC": 158000000000000
        },
        "Eg": {
            "COAST": 163000000000000, 
            "HCPS": 1200000000000000, 
            "HCSS": 847000000000000, 
            "LC": 1170000000000000
        },
        "units": "gY-1"
    }
}
```

## area_grid.npy

`NUMPY` file containing an array of shape (180, 360) and data type `float64`, which contains the calculated **ocean area** values in each grid cell. The values were obtained with the `PYPROJ` [[2](#references)] and `SHAPELY` [[3](#references)] libraries.

## biome_grid.npy

`NUMPY` file containing an array of shape (180, 360) and data type `float64`, which contains **biome IDs** of each grid cell. The IDs correspond to the following biomes:

```json
{
    "HCSS": 0,
    "LC": 1,
    "HCPS": 2,
    "COAST": 3
}
```
Nan values represent grid cells, where no data was available. Biomes were computed from Refs. [[1](#references), [4](#references) - [6](#references)].

## cnidaria_M_seed.pkl, cnidaria_Eg_seed.pkl, chordata_M_seed.pkl...

`PICKLE` file retrieved from the `Seed` object in *carbondrift.simulation.seeding*, containing **initial conditions** and **biome IDs** for each grid cell for the carcasses (M) and egestion (Eg) of different GZ phyla.

```console
        lon     lat            mass    z   origin_marker
0      -163     -77    0.000000e+00  -20               3
1      -162     -77    1.478556e+09  -20               3
2      -156     -76    1.478670e+09  -20               3
...     ...     ...             ...  ...             ...
36596   -16      83    3.649607e+09  -20               0
```

## etopo2_remaped1deg.nc, tmp_luo_remaped1deg.nc

`NETCDF` files containing **bathymetry** and **ocean temperature** fields, respectively. The datasets were obtained from Refs. [[4](#references), [7](#references)]. The files were remapped to a $1^\circ \times 1^\circ$ grid resolution using climate data operators [[8](#references)].

## References

[1] Luo, J. Y., Condon, R. H., Stock, C. A., Duarte, C. M., Lucas, C. H., Pitt, K. A., et al. (2020). Gelatinous zooplankton-mediated carbon flows in the global oceans: A data-driven modeling study. Global Biogeochemical Cycles, 34, e2020GB006704. https://doi.org/10.1029/2020GB006704

[2] Whitaker, J. (2019). pyproj, cartographic projections and coordinate transformations library. Ver. 3.7.0. URL: https://pyproj4.github.io/pyproj/stable/index.html#

[3] Gillies, Sean et al. (2007). Shapely: manipulation and analysis of geometric objects.toblerity.org. URL: https://github.com/Toblerity/Shapely

[4] NOAA National Centers for Environmental Information (2022). ETOPO 2022 15 Arc-Second Global Relief Model. NOAA National Centers for Environmental Information. DOI: https://doi.org/10.25921/fd45-gt74.

[5] Boyer Montégut, Clément de, Gurvan Madec, Albert S. Fischer, Alban Lazar in Daniele Iudicone (2004). “Mixed layer depth over the global ocean: An examination of profile data and a profile-based climatology”. V: Journal of Geophysical Research: Oceans 109.C12. DOI: https://doi.org/10.1029/2004JC002378

[6] NASA Ocean Biology Processing Group. (2025). Sea-viewing wide field-of-view sensor (seawifs) level-3 mapped ocean color data. NASA Ocean Biology Distributed Active Archive Center. Retrieved from https://oceandata.sci.gsfc.nasa.gov/directdataaccess/Level-3%20Mapped/SeaWiFS/

[7] Locarnini, R. A. et al. (2013). World Ocean Atlas 2018, Volume 1: Temperature. Ur. S.Levitus Ed. and A. Mishonov Technical Ed. NOAA Atlas NESDIS 73, 40 pp.

[8] Schulzweida, Uwe (2023). CDO User Guide. Ver. 2.3.0. DOI: 10.5281/zenodo.10020800.