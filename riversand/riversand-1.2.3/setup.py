#riversand setup.py
from distutils.core import setup   

from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()   
setup(                                
   name = 'riversand',
   version = '1.2.3',
   packages = ['riversand'],
   author = 'Konstanze Stuebner',
   author_email = 'kstueb@gmail.com',
   url = 'https://github.com/kstueb/riversand',
   description = 'Catchmentwide erosion rate calculator',
   classifiers = ["Programming Language :: Python",
                  "Programming Language :: Python :: 3",
                  "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
                  "Operating System :: OS Independent",
                  "Development Status :: 4 - Beta",
                  "Intended Audience :: Science/Research",
                  "Intended Audience :: Education",
                  "Topic :: Education",
                  "Topic :: Scientific/Engineering",
                  "Topic :: Scientific/Engineering :: GIS",
                  ],
   long_description = long_description,
   long_description_content_type = 'text/markdown',
   install_requires = ['numpy', 'scipy', 'pandas', 'xarray', 'rasterio', 'pyproj', 'fiona', 'matplotlib']
   ) 
