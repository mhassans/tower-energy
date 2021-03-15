# Find an optimized way to split the HGCal module energy over multiple "towers"

- python3 with the following packages is neeeded to run this code:
    - numpy
    - pandas
    - matplotlib
    - pyyaml
    - ROOT
    - itertools
    - math
    - root_numpy
    - scipy
    - sys
    - time
    - pathlib

(`pathlib`, `time` and probably `sys` are not very necessary and could be replaced/removed with minor changes in the code.)

## Configuration file

The config file is `config/conf_default.yaml`. In order to run or skip a function in `main.py` or `plotter.py`, only the booleans in the `mainFuncs` and/or `plotterFuncs` of the config file need to be changed, respectively. Other values in this file are just for fixing I/O names and addresses.

## `main.py`

`main.py` needs to run with the config file. Since the output of a function could be the input of another, I suggest to set all values in the `mainFuncs` of the config file to `True` when running for the first time. After that, each function can run independantly.

The command to run is:

`python3 main.py config/conf_default.yaml`

- The `param_mtx` function is the core of this analysis. It finds how the energy of HGCal modules (module sums) should split over multiple towers. The total energy is divided to (1/`N_div`)th's (e.g. 1/8th's for `N_div`=8) and split over the towers overlapping with the module. The larger the overlapping area, the more share a tower gets. The overlapping area is estimated using the number of module sensor cells enclosed by a tower. The outputs of this function are two dataframes (also called parameter matrices) showing the optimized split for all modules. A few plots are also made (explained in the code). `N_div` can be changed inside the function.

- `tower_per_module` function reads the parameter matrices produced in `param_mtx`. It then finds which towers correspond to each module and writes it into a txt file.

- Likewise, `module_per_tower` function also reads the parameter matrices from the output of `param_mtx`.  This function finds which modules correspond to each tower (hence, opposite to the `tower_per_module` function). The results are written into two txt files for CE-E and CE-H, seperately.


## `plotter.py`

This script is for producing some plots. It is better to run them after running all functions in `main.py` to make sure all of the ingredients for the plots are ready. The functions in this script can run independant of each other. One should set booleans of `plotterFuncs` in the config file to determine which functions to run.

The command to run is:

`python3 plotter.py config/conf_default.yaml`
