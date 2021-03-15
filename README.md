# Infeasible space reduction for QAOA through encoding change
Person responsible for data: *Adam Glos* (aglos [at] iitis.pl).

The scripts necessary for generating the results provided in the "Infeasible space reduction for QAOA through encoding change".

## Software installation

### Anaconda

Anaconda distribution can be downloaded from https://www.anaconda.com/products/individual

To create new environment and install required packages use

  `conda env create -f hobo_encoding_change.yml`

To activate this environment, use
  
  `conda activate hobo_encoding_change`

### Julia 

Julia can be downloaded from the official website https://julialang.org/. Version 1.5.2 was used.

In order to set up the environment, please go to the directory were *qaoa_efficiency_analysis*. Then run julia and use

  `] instantiate .`
  
This will install the required packages. The environment will be activated on its own when running scripts.

## Reproducing data

The repository already contains the data used in the publications. To generate new samples, please follow the instruction below. 

### QAOA performance


To generate data, run the following files in the *qaoa_efficiency_analysis* directory:
```
python permutation_cutter.py
julia tsp_generator.jl tsp4 4 40
julia hamilton_generator.jl hamilton 4
```
This will generate 40 tsp instances located at *tsp4*, and hamilton instaces for 3, 4  at *hamilton* number of cities used for plotting.

We use a improved version of QAOA which requires generating additional files. To generate those files, please run `julia` and type the following
```
using Pkg
Pkg.activate(".")
include("sparse_generator_loader.jl")
generator(9)
generator(6)
```


To produce data, run
```
julia tsp_qaoa_experiment.jl -p [x] -in tsp4 -outnew data
```
This will run emulator of QAOA and save the optimization results to *data/data\_[current\_time]*  together with the copy of relevant files. Instead of `[x]` please provide a number of cores. used for computation.


To plot the data please, run
```
julia tsp_qaoa_plot.jl generate plot
``` 
The line will generate and save the plot to `plots` directory.

### Error mitigation

Please enter *qaoa_efficiency_analysis* and run
```
julia tsp_generator_dict.jl 50
```
to generate samples. Then, go to the main directory and run
```
python energy_generator.py 0
python energy_generator.py 1
python energy_generator.py 2
```
To generate the energy data. Finally rn
```
python plotting_code.py
```
to generate plot.
