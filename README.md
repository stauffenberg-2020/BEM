# Blade Element Momentum(BEM) Theory
## _Simple Steady State Solver for Wind Turbine Rotors_

[![Test BEM](https://github.com/stauffenberg-2020/BEM/actions/workflows/main_test.yml/badge.svg)](https://github.com/stauffenberg-2020/BEM/actions/workflows/main_test.yml) [![Test Build](https://github.com/stauffenberg-2020/BEM/actions/workflows/build_test.yml/badge.svg)](https://github.com/stauffenberg-2020/BEM/actions/workflows/build_test.yml)

This is a *MATLAB* based Blade Element Momentum (**BEM**) Theory solver for Wind Turbine **rotors**. The idea behind this is to have a simple steady state solver for initial estimation of the performance of the wind turbine and is mainly for the beginners to understand how the core physics is, behind the curtains. This is not to be compared against advanced and transcient solvers like NREL's FAST, DTU's HAWC2, etc

## Program Features
- MATLAB Coder Compatible (MEX and EXE can be built for distribution)
- Inputs and outputs are through text files
- Wrapper scripts are friendly with DTU's HAWC2 file format

## How to run
### using the EXE
If you do not have MATLAB available, this is perhaps the simplest way. Download the EXE from the [latest release](https://github.com/stauffenberg-2020/BEM/releases/latest/download/bem_calc.exe) and add it to your system path (Alternatively, you can always navigate to the folder or provide full path to where the EXE is saved in your machine). The synax to run the executable is
```sh
bem_calc input_file.txt
```
> Note : Details about the input file format will be available in wiki soon

To run the sample NREL 5MW turbine data given in this repo, download the NREL 5MW [input_file](https://raw.githubusercontent.com/stauffenberg-2020/BEM/main/Data/NREL_5MW.txt) and its [aerofoil_file](https://raw.githubusercontent.com/stauffenberg-2020/BEM/main/Data/NREL5MWRefTurb_v50/data/NREL_5MW_pc.txt). Make sure to update the path of the [aerofoil_file] in line number 55 of the [input_file] as per your local directory where you are downloading. Always specify full file paths to prevent any errors arising due to relative paths. To run, simply use the following syntax from the folder where you have saved the NREL_5MW.txt [input file]
```sh
bem_calc NREL_5MW.txt
```
The output will be written into a text file within the same folder of the input file with "_output.txt". The output will have details of what is being outputted and its units in the header.

### using MATLAB
If you have MATLAB and want to get into the details, clone (or download) this repo to your machine. Open the home folder of the repo and make sure to add to MATLAB's path.
```sh
addpath(genpath(pwd));
```
To run the NREL 5MW turbine example given in this repo, simply run the following command in MATLAB
```sh
[output, output_details, BEM] = bem_calc('NREL_5MW.txt');
```
While running using MATLAB, you can get the additional (optional) BEM output which has all the details about the core BEM calculations.

> Note : Details about core BEM calculations and assumptions will be updated in wiki soon

## Want to contribute?
You are always welcome, Create an issue and start developing
