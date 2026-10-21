# Topography

- `surface_files/` - PlaSim surface geopotential and land-sea mask for each paleogeography
- `reconstructions/` - 720 Ma reconstructions of Merdith et al. (2021) and Li et al. (2023) on a 1 degree grid. The Merdith et al. (2021) continents were reconstructed to 720 Ma in GPlates from their plate model. The Li et al. (2023) continental outlines, from their supplementary material SM10 (https://ars.els-cdn.com/content/image/1-s2.0-S0012825223000259-mmc10.zip), were rotated to 720 Ma with its rotation file `EDRG_90W_2000-540Ma.rot` (Scenario Ia, 0-90W-0 reference frame). Both were then gridded with land at 500 m and ocean at -3000 m
- `paleodems/` - not included. To run `A4_land_area.ipynb`, download `Scotese_Wright_2018_Maps_1-88_1degX1deg_PaleoDEMS_nc.zip` from https://doi.org/10.5281/zenodo.5460860 and put its .nc files here
