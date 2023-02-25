#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 11:00:51 2023

fixures for the unit tests
these can be overwritten by fixures defined in the test_*.py files

@author: Konstanze Stübner, kstueb@gmail.com

"""

import pytest
import riversand

@pytest.fixture
def path():
    """path to user input data"""
    path = "./tests/test_data"
    return path

@pytest.fixture
def rv_single(path):
    rv = riversand.Riversand()
    rv.add_raster('dem_utm_35m.tif', path, dtype='elevation')
    rv.add_raster('toposhielding_35m.tif', path, dtype='shielding')
    rv.add_raster('quartz_35m.tif', path, dtype='quartz')    
    rv.add_samples(fname='test_samples_single.ods', path=path)
    rv.add_catchments(fname='test_single_catchment.shp', path=path)
    return rv

@pytest.fixture
def rv_multi(path):
    rv = riversand.Riversand()
    rv.add_raster('dem_utm_35m.tif', path, dtype='elevation')
    rv.add_raster('toposhielding_35m.tif', path, dtype='shielding')
    rv.add_raster('quartz_35m.tif', path, dtype='quartz')    
    rv.add_samples(fname='test_samples1.ods', path=path)
    rv.add_catchments(fname='test_multi_catchment.shp', path=path)
    return rv
