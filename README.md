Catchmentwide erosion rates with riversand 
------------------------------------------

`riversand` is a python package to calculate catchmentwide erosion rates from
cosmogenic nuclide concentrations in river sand samples. The program computes
the hypsometric statistics of the catchment area from a digital elevation model.
It uses the online erosion rate calculator by Greg Balco
(e.g. [http://stoneage.hzdr.de/](http://stoneage.hzdr.de/)) to determine predicted
nuclide concentrations $N$ for given erosion rates $E$, and calculates
the erosion rate that corresponds to the measured nuclide concentration from a
polynomial fit $N(E)$.

![Rastergrafik](https://user-images.githubusercontent.com/73031498/221909077-601ccea1-880b-4738-89d8-2ff57a16c89b.png)


The method works for in situ Be-10 and Al-26 data. It is fast (few seconds for
one catchment) for all production scaling methods implemented in the online
calculator (**St**: Lal 1991/Stone 2000; **Lm**: Lal/Stone with a geomagnetic correction
after Nishiizumi et al. 1989; **LSDn**: Lifton et al. 2014) and independent of the
catchment size or the resolution of the digital elevation model. It is considered
robust for catchments up to approx. 600 km x 600 km; for larger catchments
the effect of latitude on cosmogenic production may become significant.

The approach is described in:

Stübner, K., Balco, G., and Schmeisser, N. (in review). Calculating catchmentwide
erosion rates using an existing online calculator. *Radiocarbon*.

Definitely check out the documentation of the online calculator (e.g. [here](http://stoneage.ice-d.org/math/docs/v3/v3_input_explained.html)
or [here](https://sites.google.com/a/bgc.org/v3docs/)) and the publication
[Balco et al. (2008)](http://hess.ess.washington.edu/math/docs/al_be_v2/al_be_calc_2007.pdf)
before using this calculator.

Documentation
-------------
- [`quickstart.ipynb`](https://github.com/kstueb/riversand/blob/main/riversand/example_scripts/quickstart.ipynb)
- [`step_by_step.ipynb`](https://github.com/kstueb/riversand/blob/main/riversand/example_scripts/step_by_step.ipynb)
- [`test_data/`](https://github.com/kstueb/riversand/blob/main/riversand/example_scripts/test_data) : geotiffs of a 35m-resolution digital elevation model, a topographic shielding raster generated with [TopoToolbox](https://topotoolbox.wordpress.com/) and a binary raster indicating quartz-bearing and quartz-free lithologies; shapefiles with catchment outlines; a spreadsheet with sample data. 

Installation
------------
Install riversand by running:
```
$ pip install riversand
```
Requirements
------------
- numpy, scipy, pandas, xarray
- rasterio, fiona, pyproj
- matplotlib

License
-------
[GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html)
