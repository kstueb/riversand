Catchmentwide Erosion Rate Calculator for in situ cosmogenic nuclides 
---------------------------------------------------------------------

This calculator processes geospatial data (digital elevation model, catchment
outline and other) to extract hypsometry statistics of the catchment and
determine a cosmogenic nuclide catchmentwide erosion rate. It uses the online
erosion rate calculator by Greg Balco (e.g. http://stoneage.hzdr.de/) to
calculate cosmogenic nuclide production.

The method works for in situ Be-10 and Al-26 data and is considered
robust for catchments up to approx. 600 km x 600 km; for larger catchments
the effect of latitude on cosmogenic production may become significant.

Citation
--------

The software is described in:

Stübner, K., Balco, G., and Schmeisser, N. (in review). Calculating catchmentwide erosion rates using an existing online calculator. *Radiocarbon*. 

Installation and Usage
----------------------

This `riversand` package has been developed in python 3.9. It needs the
geospatial data processing packages `rasterio`, `fiona` and `pyproj` as 
well as several other common python packages including `scipy`, `xarray`
and `pandas`.

The package is available at PyPI.org and can be installed with
```
> pip install riversand
```
A `conda` distribution will be released in the future.

To get started, copy the  folder [`example_scripts`](https://github.com/kstueb/riversand/tree/main/riversand/example_scripts)
to your computer and go through the example jupyter notebooks
[`quickstart.ipynb`](https://github.com/kstueb/riversand/blob/main/riversand/example_scripts/quickstart.ipynb) and
[`step_by_step.ipynb`](https://github.com/kstueb/riversand/blob/main/riversand/example_scripts/step_by_step.ipynb)
The subfolder `user_data` has some example datasets including geotiff's
of the catchment DEM and topographic shielding factors calculated with
[TopoToolbox](https://topotoolbox.wordpress.com),
shapefiles with catchment outlines and a spreadsheet with sample data).
