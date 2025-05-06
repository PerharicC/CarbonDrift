#!/bin/bash

#Create directories to save simulation files
if [ ! -d "exp_constant" ]; then
	mkdir exp_constant
fi
if [ ! -d "exp_variable" ]; then
	mkdir exp_variable
fi

#Run simulations
echo "Start simulations."
python run_all.py

#Make plots
if [ ! -d "images" ]; then
	mkdir images
fi
echo "Running simulation analysis."

#Make carbon transport tables
echo "Calculating carbon transport (may take some time)"
python carbon_transport_tables.py

#Plotting
echo "Making carbon flux global maps"
plot_simulation @flux_map_params.txt -d -100 -o images/mass_exp_variable_100.pdf
plot_simulation @flux_map_params.txt -d -1000 -o images/mass_exp_variable_1000.pdf
plot_simulation @flux_map_params.txt -d -6000 -o images/mass_exp_variable_seafloor.pdf

echo "Making carbon flux global map differences"
plot_simulation @flux_map_params_diff.txt -d -100 -o images/mass_diff_100.pdf
plot_simulation @flux_map_params_diff.txt -d -1000 -o images/mass_diff_1000.pdf
plot_simulation @flux_map_params_diff.txt -d -6000 -o images/mass_diff_seafloor.pdf

echo "Plotting trajectories"
plot_simulation @z_m_params.txt
plot_simulation @z_t_params.txt

echo "Plotting longitude averaged fluxes"
plot_simulation @mean_lon_params.txt
