# Inferring the natural history of HPV from global cancer registries: insights from a multi-country calibration

## Summary
This repository includes the code for generating the figures and results in the manuscript "Inferring the natural history of HPV from global cancer registries: insights from a multi-country calibration". The citation for the manuscript is:

> Inferring the natural history of HPV from global cancer registries: insights from a multi-country calibration (preprint). Stuart RM, Cohen JA, Abeysuriya RG, Sanz-Leon P, Kerr CC, Rao D, Klein DJ.


## Structure
- The scripts prefixed `run_` are used for running simulations and calibrations.
- The scripts prefixed `plot_` are used for plotting figures. The numbering of the figures corresponds to the numbering in the manuscript.
- The `data` folder contains all the input data used in the analyses
- The `results` folder contains the finalized parameter sets from the immunovarying calibration.
- The script `calibration.py` contains the methods and classes for running MultiCals, as used in the constrained calibrations within the manuscript.
- The scripts `locations.py` and `utils.py` contain the list of country locations and various utilities, respectively.

## Installating and running
- Ensure you have Python installed (if you haven't installed Python already, the easiest is to use Anaconda).
- Install HPVsim (hpvsim.org)
- Clone this repository and run the scripts as desired.
