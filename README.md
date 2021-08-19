# [The title]

Person responsible for data: Adam Glos (aglos [at] iitis.pl).

The scripts necessary for generating the results provided in the "[The title]"

The code was used on Ubuntu OS 20.04. It is not guarantee it will work on other operating systems

# Software instalation

## Anaconda

Anaconda distribution can be downloaded from https://www.anaconda.com/products/individual

  `conda env create -f hobo_encoding_change.yml`

To activate this environment, use
  
  `conda activate hobo_encoding_change`

To deactivate an active environment, use
  
  `conda deactivate`

## Julia

Julia can be downloaded from the official website https://julialang.org/. Version 1.6.1 was used.

In order to set up the environment, please go to the directory ec-qaoa-code. Then run julia and use

```
julia> ]
(@v1.6) activate .
(@v1.6) instantiate
```

This will install the required packages based on _Manifest.toml_ file. The environment will be activated on its own when running scripts.

# Reproducing data

The repository already contains the data used in the publications. To generate new samples, please follow the instruction below. Occasionally the code omits generating new data if old already exist (the scripts will complain, or inform about such issue) - in such case please remove the old data.

## TSP generation
To generate data, run the following files in the _qaoa_efficiency_analysis_ directory:

```
julia tsp_generator.jl tsp_data_xy/tsp3 3 100
julia tsp_generator.jl tsp_data_xy/tsp4 4 100
```

This will generate 100 TSP for 3 and 4 cities located at the corresponding repositories.

## Experiment data generation
To generate data for random angles for XY-QAOA and GM-QAOA, please run in the main directory
```
./xy_experiment_generator.sh 
```
This will generate data of appropriate size used in the paper. To parallelize the compution, in the bash script we used `&` so that the next command runs after the previous would was _started_. This will likely produce to large number of processes for a regular PC. To surpass overloading the computer, please remove `&`. All data

The commands run through this script takes form
```
python energy_diff_generator.py $VAR $CITYNO 40 1 100 rand $GAMMA
```
where `VAR` is the type of error model used for computation, `CITYNO` is the number of cities, `40` is the number of layers considered, `1` is the number of times a single TSP instances should be considered `100` is the number of different TSP instances, `rand` is the type of angles used (the only option working) and `GAMMA` is the number of noise strength chosen.

Furthermore. the script runs
```
python optimization_qaoa.py 3 40
python optimization_qaoa.py 4 40
```
which runs the QAOA algorithm for 3/4 cities and 40 TSP instances (each considered once).  To reduce the number of parallel processes, please change the `pool_size` in the file _optimization_qaoa.py_.


For Energy difference computation, data will be saved to _data_noise_qaoa_xy_gm/tsp3_plotting_for_paper_  and _data_noise_qaoa_xy_gm/tsp4_plotting_for_paper_ depending on the number of cities considered. For QAOA optimization, the data will be saved to _data_noise_qaoa_xy_gm/optimization_.


## Plots generation

To generate plots, please run the following commands
```
python energy_diff_plot.py 40 1 100 rand
python optimization_plot.py
```
The last command may require slightly more time than the previous commands, as the energies has to be recomputed. The plots will be saved to _plots_.
