# snowball-susceptibility

[![DOI](https://zenodo.org/badge/1427752136.svg)](https://doi.org/10.5281/zenodo.23522467)

This repository contains the code and data used to generate the figures for **"Relative Impact of Solar Luminosity vs. Paleogeography on Sensitivity to Snowball Transitions Over Geologic Time."** It combines PlaSim GCM simulations of modern, Carboniferous, and Neoproterozoic paleogeographies with a box model of the carbon cycle. As a matter of courtesy, we request that people using this code please cite Minsky & Wordsworth (in review).

## Installation

Ensure you have Python 3.13 and Jupyter installed. Install required dependencies using:

```bash
pip install -r requirements.txt
```

## Code structure

- `helpers.py` - Loading GCM output, background and Snowball threshold CO2, and plotting utilities
- `weathering.py` - WHAK and MAC silicate weathering
- `gcm_model.py` - Carbon cycle model coupled to the GCM output
- `box_model.py` - Idealized carbon cycle box model

## Notebooks and figure outputs

The numeric prefixes match figure numbers.

### Main text

- `1_maps.ipynb` - Produces Fig. 1: continental configurations
- `2_gcm_ensemble.ipynb` - Produces Fig. 2: GCM ensemble for the modern configuration
- `3_co2_results.ipynb` - Produces Fig. 3, Fig. 4, and Fig. 6, and the values quoted in the Results
- `5_carbon_cycle.ipynb` - Produces Fig. 5, Fig. A2, and Fig. A3: carbon cycle response to perturbations

### Appendix

- `A1_weathering_sensitivity.ipynb` - Produces Fig. A1: weathering at the background CO2
- `A4_land_area.ipynb` - Produces Fig. A4: land area

## Data

`data/` contains the 10-year mean surface temperature and runoff of every equilibrated PlaSim run, with one file per ocean heat transport case and paleogeography; the raw PlaSim output is not included due to its size. `topography/` contains the PlaSim surface files and the 720 Ma reconstructions; its README says where the data come from.

## Contact

Contact: cminsky@g.harvard.edu
