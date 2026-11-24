# Case 1: Run and visualize

Use MATLAB R2024a on Windows. Open this folder as the MATLAB current folder.

## Run

```matlab
result = run_case1('GFDHS', 0.05, 1, 'results/gfdhs_run1');
```

The arguments are the method, spatial setting in metres, random seed, and output folder. The output folder must not already exist.

Methods: `GFDHS`, `GBH`, `Astar`. Supported spatial settings: `0.01`, `0.05`, `0.10`, `0.15`, `0.20`, `0.25`.

```matlab
run_case1('GBH',   0.25, 1, 'results/gbh_run1');
run_case1('Astar', 0.25, 1, 'results/astar_run1');
```

The entry point loads the files in `data/` and selects the required functions automatically. Do not add all subfolders to the MATLAB path. Start with the examples above; fine-grid runs require more memory.

## Plot

```matlab
plot_case1('results/gfdhs_run1/result.mat');
```

Results are saved in the specified output folder as MAT/JSON files and route-coordinate CSV files. Exported coordinates are in metres.
