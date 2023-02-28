Catchmentwide Erosion Rates with riversand 
------------------------------------------

`riversand` is a python package to calculate catchmentwide erosion rates from
cosmogenic nuclide concentrations in river sand samples. The program computes
the hypsometric statistics of the catchment area from a digital elevation model.
It uses the online erosion rate calculator by Greg Balco
(e.g. [http://stoneage.hzdr.de/](http://stoneage.hzdr.de/)) to determine predicted
nuclide concentrations $N$ for given erosion rates $\varepsilon$, and calculates
the erosion rate that corresponds to the measured nuclide concentration from a
polynomial fit to $N(\varepsilon)$.

![Rastergrafik](https://user-images.githubusercontent.com/73031498/221909077-601ccea1-880b-4738-89d8-2ff57a16c89b.png)


The method works for in situ Be-10 and Al-26 data. It is fast (few seconds for
one catchment) for all production scaling methods implemented in the online
calculator (St: Lal 1991/Stone 2000; Lm: Lal/Stone with a geomagnetic correction
after Nishiizumi et al. 1989; LSDn: Lifton et al. 2014) and independent of the
catchment size or the resolution of the digital elevation model. It is considered
robust for catchments up to approx. 600 km x 600 km; for larger catchments
the effect of latitude on cosmogenic production may become significant.

The approach is described in:

Stübner, K., Balco, G., and Schmeisser, N. (in review). Calculating catchmentwide
erosion rates using an existing online calculator. *Radiocarbon*. 

Documentation
-------------


Installation
------------
Install riversand by running:
```
$ pip install riversand
```
Requirements
------------
- numpy, scipy, pandas
- rasterio
- fiona
- pyproj
- xarray


To get started, copy the  folder [`example_scripts`](https://github.com/kstueb/riversand/tree/main/riversand/example_scripts)
to your computer and go through the example jupyter notebooks
[`quickstart.ipynb`](https://github.com/kstueb/riversand/blob/main/riversand/example_scripts/quickstart.ipynb) and
[`step_by_step.ipynb`](https://github.com/kstueb/riversand/blob/main/riversand/example_scripts/step_by_step.ipynb)
The subfolder `user_data` has some example datasets including geotiff's
of the catchment DEM and topographic shielding factors calculated with
[TopoToolbox](https://topotoolbox.wordpress.com),
shapefiles with catchment outlines and a spreadsheet with sample data).
