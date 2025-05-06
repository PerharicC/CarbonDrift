# Guide to reproduce global GZ carbon transport results

Open the terminal in this directory and run the command
```console
:~$ bash run.sh
```

**Note:** This will produce and plot simulations for the **mass-dependent decay with constant and variable sinking speeds with an exponential decay rate**.

If you would like to see **linear decay rate** results or **area_dependent decay** you should change the -mdtype in *params_all.txt* from mass to area and decay_types in *run_all.py** from ["exp"] to ["linear"].