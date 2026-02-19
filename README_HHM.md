# Hybrid Harmonic Model (HHM)

This repository contains the MATLAB script `Parfor_HHM.m` used to generate long time series of Total Water Level (TWL) for a selected site.

## What the script does

For each SLR pathway and rate scenario, the script constructs TWL as:

Astronomical tides reconstructed from harmonic constituents

Plus mean sea level change from an SLR curve

Plus storm surge and non tidal residual (NTR) from 100 Monte Carlo simulations

Plus optional wind wave setup and runup (switch controlled)

## User editable settings

At the top of `Parfor_HHM.m`, edit only the block labelled `User settings`.

Key parameters

`t_tide_folder`

Path to the t tide toolbox containing `t_tide` and `t_predic`

`project_folder`

Working folder where the model files and MAT inputs exist

`computewindwaves`

Set to 1 to include wind wave setup and runup

Set to 0 to compute TWL as water level plus NTR only

`save_outputs`

Set to 1 to save a MAT output file in an `outputs` folder

## Expected folder structure

The script assumes it is running inside `project_folder` and uses relative paths.

Example

project_folder

GMSLR_1m

s2020_r20

output

FlowFM_0000_his.nc

s2040L_r20

output

FlowFM_0000_his.nc

...

MAT files in project_folder

all_montecarlo_SS.mat

montecarlo_time.mat

RSmat_1_2m_newbathy_v2.mat

RSmat_2_3point3m_newbathy_v2.mat

Inpaint_nans

Inpaint_nans

## Input data

### Delft3D FM history NetCDF

The script reads `waterlevel` from `FlowFM_0000_his.nc` and extracts a single mesh index `idx` for the selected location.

### Storm surge Monte Carlo

`all_montecarlo_SS.mat` must contain `ss_data_matrix`

`montecarlo_time.mat` must contain `time_ss`

### Wind wave setup and runup tables

Used only when `computewindwaves` equals 1.

`RSmat_1_2m_newbathy_v2.mat`

`RSmat_2_3point3m_newbathy_v2.mat`

## Outputs

The script creates the following MATLAB workspace variables

`WL`

Astronomical water level plus SLR and datum shift for each rate scenario

`Twl_r20`, `Twl_r37`, `Twl_r52`, `Twl_r100`

Cell arrays of TWL. Each cell corresponds to low, median, high SLR. Each cell contains an Nt by Nsims matrix.

Optional

`RS_timeseries`

Runup and setup time series used in TWL when wind waves are enabled

### Saved MAT file

If `save_outputs` equals 1, a MAT file is written to

`project_folder/outputs/HHM_outputs_<loc>_idx<idx>_<timestamp>.mat`

This MAT file contains `WL`, `t100`, TWL variables, and a `meta` struct documenting key settings.

## Dependencies

MATLAB

t tide toolbox

Parallel Computing Toolbox is optional. If not available, the script runs serially.

## Citation

If you use this code, cite the associated paper and include the GitHub commit hash used to generate results.
