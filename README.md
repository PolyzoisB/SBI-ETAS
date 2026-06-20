# SBI-ETASI Inference Pipeline

Likelihood-free inference pipeline for estimating parameters of an incomplete ETAS-type seismicity model with short-term aftershock incompleteness (ETASI).
The workflow compares observed and simulated earthquake catalogs through normalized space-time-magnitude summary statistics rather than evaluating the ETAS/ETASI
 likelihood directly.

The code is written in Fortran and is organized into two main stages:

1. **Preprocessing** of the observed catalog.
2. **Simulation-based inference** using ETASI catalog simulations and summary-statistic matching.

---

## Repository structure

```text
.
├── config.f90                  # Global configuration, file paths, thresholds, model constants
├── si_declustering.f90          # Susceptibility-index / nearest-neighbor declustering tools
├── nn_cluster.f90               # Nearest-neighbor clustering and mainshock identification
├── space_time_mag_count.f90     # Space-time-magnitude foreshock/aftershock count statistics
├── model_sims.f90               # ETASI catalog simulation routines
├── preprocess_step.f90          # Preprocessing executable program
├── sbi_estimation.f90           # SBI / likelihood-free inference executable program
├── run_preprocess.sh            # Compile and run preprocessing
├── run_sbi_inference.sh         # Compile and run inference
├── datasets/                    # Input catalogs and background-coordinate files
└── results/                     # Output summary statistics and inferred parameter sets
```

---

## Method summary

The pipeline implements a likelihood-free inference strategy for incomplete seismicity catalogs.
The core idea is to avoid direct likelihood evaluation and instead compare simulated and observed catalogs through diagnostic summaries.

### Preprocessing stage

The preprocessing step:

1. Reads an earthquake catalog with columns

   ```text
   elapsed_time  latitude  longitude  depth  magnitude
   ```

2. Keeps events with magnitude `M >= mc`, where `mc` is defined in `config.f90`.
3. Estimates the Gutenberg-Richter `b`-value.
4. Applies a nearest-neighbor / susceptibility-index declustering procedure to estimate:
   - the nearest-neighbor threshold,
   - the number of background events,
   - the background rate,
   - background-event coordinates.
5. Identifies mainshock classes and computes normalized aftershock and foreshock summary statistics.
6. Exports the catalog-level quantities and summary statistics used by the inference step.

### Inference stage

The inference step:

1. Reads the observed aftershock and foreshock summary statistics.
2. Reads catalog-level quantities estimated during preprocessing, including:
   - `b`-value,
   - background rate,
   - nearest-neighbor threshold,
   - branching-ratio bounds,
   - accepted catalog-size bounds.
3. Initializes ETASI parameters under physical constraints, including a branching-ratio constraint.
4. Repeatedly simulates `K1` ETASI catalogs for each proposed parameter set.
5. Computes the same summary statistics for each simulated catalog.
6. Averages simulated summaries and compares them with the observed summaries using a relative discrepancy cost.
7. Updates parameter blocks and accepts proposals that reduce the cost.
8. Exports the final parameter set, initial parameter set, final cost, elapsed time, and number of iterations.

---

## Model parameters

The inference vector contains 11 parameters:

| Index | Name | Meaning | Status in current code |
|---:|---|---|---|
| 1 | `p` | Omori-Utsu temporal decay exponent | estimated |
| 2 | `c` | Omori-Utsu time offset, in days in the parameter vector | estimated |
| 3 | `alpha` | Magnitude-productivity exponent | estimated |
| 4 | `K` | Productivity scaling parameter | estimated |
| 5 | `d` | Spatial scale parameter, stored as `log10(d)` | estimated |
| 6 | `gamma` | Magnitude dependence of spatial scale | estimated |
| 7 | `q` | Spatial-kernel decay exponent | estimated |
| 8 | `tau-ETASI` | Short-term incompleteness time window, seconds | currently fixed at 0 s in initialization/update blocks |
| 9 | `dr-ETASI` | Short-term incompleteness spatial window, km | currently fixed at 50 km |
| 10 | `bg-rate` | Background rate, events/sec/deg² | fixed from preprocessing |
| 11 | `b-value` | Gutenberg-Richter b-value | fixed from preprocessing |

The temporal parameter `c` is stored in days in the parameter vector and converted to seconds inside the simulator.

---

## Configuration

Most user-facing settings are defined in `config.f90`.

Important configuration blocks include:

- **Input/output paths**
  - `input_catalog`
  - `input_aft_summary_stats`
  - `input_fore_summary_stats`
  - `input_catalog_stats`
  - `true_model_stats`
  - `bg_coords`

- **Magnitude and domain settings**
  - `mc`: magnitude cutoff
  - `msup`: upper magnitude bound
  - `lat_min`, `lat_max`, `lon_min`, `lon_max`: spatial domain
  - `tc`: auxiliary time window excluded from the target catalog
  - `tlast`: total simulated time window

- **Summary-statistic bins**
  - `thtime`: aftershock time windows
  - `thml`: aftershock magnitude bins
  - `thtimef`: foreshock time windows
  - `thspacef`: foreshock spatial windows
  - `thmlf`: foreshock magnitude bins

- **Monte Carlo / optimization settings**
  - `K0`: maximum optimization iterations
  - `K1`: number of simulations used to average summary statistics per parameter proposal
  - `epsilon`: convergence tolerance
  - `conv_thr`: number of stable iterations required for early stopping
  - `lr`: proposal step size in normalized parameter space

---

## Requirements

The pipeline requires:

- `gfortran`
- POSIX-compatible shell, e.g. `sh` or `bash`
- Input catalog placed under the path specified by `input_catalog`
- Writable output directories, especially `results/` and any subdirectories used in `config.f90`

Example installation on Ubuntu/Debian:

```bash
sudo apt update
sudo apt install gfortran
```

---

## Input data format

The earthquake catalog should be a whitespace-separated text file with five columns:

```text
elapsed_time  latitude  longitude  depth  magnitude
```

where:

- `elapsed_time` is in seconds,
- `latitude` and `longitude` are in decimal degrees,
- `depth` is read but not used in the current summary statistics,
- `magnitude` is used for catalog filtering, clustering, and ETASI simulation matching.

The default catalog path is set in `config.f90`:

```fortran
character(len=*), parameter :: input_catalog = 'datasets/sc_1981_31_8_2023_all_conv.txt'
```

---

## Running the pipeline

### 1. Create output directories

Before running, make sure that the directories referenced in `config.f90` exist. For example:

```bash
mkdir -p datasets
mkdir -p results/sc_mc2.5
```

Adjust the directory names if your `config.f90` uses different paths.

### 2. Run preprocessing

```bash
sh run_preprocess.sh
```

This compiles:

```text
config.f90
si_declustering.f90
nn_cluster.f90
space_time_mag_count.f90
preprocess_step.f90
```

then runs the preprocessing executable and removes temporary object and module files.

Expected preprocessing outputs include:

```text
true_aft_sum_stats.txt
results/sc_mc2.5/true_fore_sum_stats.txt
results/sc_mc2.5/true_catalog_stats.txt
results/sc_mc2.5/true_model_stats.txt
datasets/bg_coords_sc_mc3.0.txt
```

The exact paths depend on the values set in `config.f90`.

### 3. Run simulation-based inference

```bash
sh run_sbi_inference.sh
```

This compiles:

```text
config.f90
nn_cluster.f90
space_time_mag_count.f90
model_sims.f90
sbi_estimation.f90
```

then runs the SBI/ETASI inference executable and removes temporary object and module files.

The inference output is written to:

```text
results/sc_mc2.5/params_<seed>.txt
```

where `<seed>` is the seed read by `sbi_estimation.f90`. The provided shell script does not pass a seed explicitly, so the code uses its internal default seed unless you run the executable manually with a command-line argument.

To run manually with an explicit seed, compile the files and execute, for example:

```bash
./summary_stats.exe 4045
```

or rename the executable in the shell script to a clearer name such as `sbi_estimation.exe`.

---

## Output files

### `true_aft_sum_stats.txt`

Observed aftershock summary statistics. Each row contains:

```text
normalized_aftershock_count  number_of_mainshocks  mainshock_class  time_window_days  spatial_window_km  magnitude_bin_left_edge
```

### `true_fore_sum_stats.txt`

Observed foreshock summary statistics. Each row contains:

```text
normalized_foreshock_count  number_of_mainshocks  mainshock_class  time_window_days  spatial_window_km  magnitude_bin_left_edge
```

### `true_catalog_stats.txt`

Catalog-level quantities used by the inference step:

```text
b-value
background_rate
nearest_neighbor_threshold
branching_ratio_upper_bound
branching_ratio_lower_bound
catalog_size_upper_bound
catalog_size_lower_bound
```

### `true_model_stats.txt`

Human-readable summary of preprocessing choices and estimated catalog quantities.

### `params_<seed>.txt`

Final inference output. The row contains:

```text
final_parameter_set  initial_parameter_set  final_cost  elapsed_time_seconds  number_of_iterations
```

---

## Summary statistics

The current summaries are normalized counts of events around mainshocks, separated by:

- mainshock magnitude class:
  - `M4`
  - `M5`
  - `M6+`
- event magnitude class,
- time window,
- spatial window.

Aftershock summaries are computed forward in time from each mainshock. Foreshock summaries are computed backward in time. Counts are normalized by the number of mainshocks in each mainshock class.

A dynamic spatial window is also used for the largest aftershock spatial bin:

```fortran
R(M) = 0.01 * 10**(0.5*M)
```

where `M` is the mainshock magnitude.

---

## Inference objective

For each candidate parameter set, the code simulates `K1` catalogs and computes average summary statistics. The current cost compares simulated and observed aftershock summaries using a relative squared discrepancy:

```text
J(theta) = sum_i [ S_sim_i(theta) / S_obs_i - 1 ]^2
```

where the sum is over nonzero observed summary-statistic cells.

Foreshock summaries are read by the inference program but the foreshock contribution to the cost is currently commented out. It can be re-enabled in `sbi_estimation.f90` if foreshock statistics are intended to constrain the model.

---

## Branching-ratio and catalog-size constraints

Parameter initialization and block updates are constrained by the branching ratio. Candidate parameter sets are rejected if their branching ratio falls outside the interval defined by:

```text
br_inf <= n < br_sup
```

The simulation step also rejects simulated catalogs whose number of events falls outside the catalog-size bounds estimated during preprocessing:

```text
n_inf <= N_sim <= n_sup
```

These checks prevent unstable or physically unrealistic simulated catalogs from dominating the inference.

---

## Reproducibility

The random-number generator uses an integer seed. The inference program attempts to read the seed from the first command-line argument. If the argument is missing or cannot be read, a default seed is used.

For reproducible multi-start runs, launch the inference executable repeatedly with different seeds and keep the resulting `params_<seed>.txt` files.

Example:

```bash
for seed in 4001 4002 4003 4004 4005; do
    ./summary_stats.exe "$seed"
done
```

---

## Citation

If you use this code, please cite the associated SBI/ETASI manuscript and the related nearest-neighbor declustering and susceptibility-index papers:

1. Bountzis, P., Petrillo, G., & Lippiello, E. (2026). *A likelihood-free framework for seismic forecasting models: application to incomplete ETAS model via simulation-based inference*. Authorea Preprints / ESS Open Archive. Under major revision in *JGR: Solid Earth*. https://doi.org/10.22541/essoar.177315051.18116847/v1

2. Bountzis, P., Lippiello, E., Baccari, S., & Petrillo, G. (2026). Automatic earthquake declustering using the nearest-neighbor distance. *Earth and Space Science, 13*, e2025EA004539. https://doi.org/10.1029/2025EA004539

3. Lippiello, E., Baccari, S., & Bountzis, P. (2023). Determining the number of clusters, before finding clusters, from the susceptibility of the similarity matrix. *Physica A: Statistical Mechanics and its Applications, 616*, 128592. https://doi.org/10.1016/j.physa.2023.128592

### BibTeX

```bibtex
@misc{bountzis2026likelihoodfree,
  author       = {Bountzis, Polyzois and Petrillo, Giuseppe and Lippiello, Eugenio},
  title        = {A likelihood-free framework for seismic forecasting models: application to incomplete {ETAS} model via simulation-based inference},
  year         = {2026},
  publisher    = {Authorea Preprints / ESS Open Archive},
  doi          = {10.22541/essoar.177315051.18116847/v1},
  note         = {Under major revision in JGR: Solid Earth}
}

@article{bountzis2026automatic,
  author  = {Bountzis, Polyzois and Lippiello, Eugenio and Baccari, Silvio and Petrillo, Giuseppe},
  title   = {Automatic earthquake declustering using the nearest-neighbor distance},
  journal = {Earth and Space Science},
  year    = {2026},
  volume  = {13},
  pages   = {e2025EA004539},
  doi     = {10.1029/2025EA004539}
}

@article{lippiello2023determining,
  author  = {Lippiello, Eugenio and Baccari, Silvio and Bountzis, Polyzois},
  title   = {Determining the number of clusters, before finding clusters, from the susceptibility of the similarity matrix},
  journal = {Physica A: Statistical Mechanics and its Applications},
  year    = {2023},
  volume  = {616},
  pages   = {128592},
  doi     = {10.1016/j.physa.2023.128592}
}
```

---

## License

This project is released under the MIT License. See the [`LICENSE`](LICENSE) file for details.
