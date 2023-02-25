#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 08:39:59 2023

@author: Konstanze Stübner, kstueb@gmail.com

"""

import pytest
from pytest import mark

import sys, os
import pandas as pd
import numpy as np

import rasterio
import fiona
from pyproj.crs import CRS
#import pyproj.exceptions
import xarray as xr


import riversand
from riversand.utils import validate_topo, validate_nuclide, validate_sample
from riversand.utils import feature_in_raster, clip_raster, eliminate_quartzfree
from riversand.geospatial import Raster, Catchment
from riversand.utils import get_xarray_centroid, get_polygon_centroid, projected_xy_to_longlat
from riversand.utils import get_bins

def test_params():
    assert riversand.params.url == "http://stoneage.hzdr.de/cgi/matweb"
    

# =============================================================================
# Validate geospatial objects
# =============================================================================

class Test_Raster():
    @mark.parametrize("fname, dtype, epsg, res",
                      [('test.tif', 'elevation', 32632, (500,500)),
                       ('test.tif', 'Shielding', 32632, (500,500)),
                       ('test.tif', 'QUARTZ',    32632, (500,500)),
                      ])
    def test_Raster(self, path, fname, dtype, epsg, res):
        """class Raster, valid input"""
        R = Raster(os.path.join(path, fname), dtype)
        assert R.dtype==dtype.lower() 
        assert R.fname==(os.path.join(path, fname))
        assert isinstance(R.src, rasterio.DatasetReader)
        assert R.epsg==epsg 
        assert R.res==res
    
    def test_Raster_Exceptions(self, path):
        """class Raster, invalid input"""
        
        with pytest.raises(TypeError) as e:
            Raster('test.tif', 'bogus')
        assert e.value.args[0] == "Invalid Raster dtype 'bogus' (must be 'elevation', 'shielding' or 'quartz')"
    
        with pytest.raises(TypeError) as e:
            Raster(0, 'quartz')
        assert e.value.args[0] == "Invalid Raster fname (must be string)"
    
        with pytest.raises(FileNotFoundError) as e:
            Raster(os.path.join('', 'test.tif'), 'quartz')
        assert e.value.args[0] == "test.tif: No such file or directory"
    
        with pytest.raises(OSError) as e:
            Raster(os.path.join(path, 'test.txt'), 'quartz')
        assert e.value.args[0] == "Cannot read geotiff ../test_data/test.txt" # raised by Raster.__init__()

    def test_add_raster_Exceptions(self, path):
        R = Raster(os.path.join(path, 'test.tif'), dtype='quartz')
        rv = riversand.Riversand()
        with pytest.raises(NotImplementedError) as e:
            rv.add_raster(R)
        assert e.value.args[0] == "add_raster() not yet implemented for Raster object"

        
class Test_Catchment():
    def test_Catchment1(self, path):
        """class Catchment, valid input"""
        C = Catchment(os.path.join(path, 'test_single_catchment.shp'))
        assert C.attrs==['name', 'id']
        assert len(C.catchments)==1
        assert isinstance(C.crs, CRS)
        assert C.dtype=='single'
        assert C.fname=='../test_data/test_single_catchment.shp'
        assert C.crs.to_epsg()==32632
        assert isinstance(C.src, fiona.Collection)
        
        assert C.cid is None
        C.set_cid()
        assert C.cid=='id'
        C.set_cid('name')
        assert C.cid=='name'
        C.set_cid('bogus')
        assert C.cid is None
        
    
    def test_Catchment2(self, path):
        """class Catchment, valid input"""
        C = Catchment(os.path.join(path, 'test_multi_catchment.shp'), 'multi')
        assert C.attrs==['name', 'id', 'area_km2', 'test']
        assert len(C.catchments)==8
        assert C.cid is None
        assert isinstance(C.crs, CRS)
        assert C.dtype=='multi'
        assert C.fname=='../test_data/test_multi_catchment.shp'
        assert C.crs.to_epsg()==32632
        assert isinstance(C.src, fiona.Collection)
            
    
    def test_Catchment_Exceptions(self, path):
        """class Catchment, invalid input"""
        with pytest.raises(TypeError) as e:
            Catchment('test.shp', 'bogus')
        assert e.value.args[0] == "Invalid Catchment dtype 'bogus' (must be 'single' or 'multi')"
    
        with pytest.raises(TypeError) as e:
            Catchment(0, 'single')
        assert e.value.args[0] == "Invalid Catchment fname (must be string)"
    
        with pytest.raises(FileNotFoundError) as e:
            Catchment(os.path.join('', 'test_multi_catchment.shp'), 'multi')
        assert e.value.args[0] == "test_multi_catchment.shp: No such file or directory"
    
        with pytest.raises(OSError) as e:
            Catchment(os.path.join(path, 'test.txt'), 'multi')
        assert e.value.args[0] == "Cannot read shapefile ../test_data/test.txt" # raised by Catchment.__init__()

    def test_add_catchments_Exceptions(self, path):
        C = Catchment(os.path.join(path, 'test_single_catchment.shp'))
        rv = riversand.Riversand()
        with pytest.raises(NotImplementedError) as e:
            rv.add_catchments(C)
        assert e.value.args[0] == "add_catchments() not yet implemented for Catchment object"
        
    def test_Catchment_set_cid(self, path):
        
        C = Catchment(os.path.join(path, 'test_single_catchment.shp'), 'single')
        assert C.cid is None
        
        C.set_cid('id')
        assert C.cid=='id'
        
        C.set_cid(None)
        assert C.cid is None
        
        C.set_cid('this')
        assert C.cid is None
    
        C.set_cid('id')
        C.set_cid('this')
        assert C.cid==None
        
        
    def test_Catchment_get_names1(self, path):
        """ .get_names() """
        
        # from different attribute fields
        C = Catchment(os.path.join(path, 'test_multi_catchment.shp'))
        c_names = C.get_names('name')
        assert c_names==['DB02', 'DB03', 'DB03', 'DB04', 'DB05', 'DB12', 'DB17', 'DB19']
        c_names = C.get_names('id')
        assert c_names==['12', '17', '19', '3', '33', '4', '5', 'None', ]
        c_names = C.get_names('test') # empty strings/white spaces are also 'None'
        assert c_names==['None', 'None', 'None', 'a', 'b', 'c', 'c', 'c']
        
    def test_Catchment_get_names2(self, path):
        """ raise Exception if cid not in attrs """
        
        C = Catchment(os.path.join(path, 'test_multi_catchment.shp'))
        with pytest.raises(ValueError) as e:
            C.get_names('bogus')
        assert e.value.args[0]==("Invalid catchment identifier 'bogus';\n"+
                                  "   attribute fields are: name, id, area_km2, test")
        
    def test_Catchment_get_names3(self, path):
        """ raise Exception if cid not in attrs """
        
        C = Catchment(os.path.join(path, 'test_multi_catchment.shp'))
        with pytest.raises(ValueError) as e: # accepts None as input
            C.get_names(None)
        assert e.value.args[0]==("Invalid catchment identifier 'None';\n" +
                                  "   attribute fields are: name, id, area_km2, test")
        
    def test_Catchment_get_valid_names1(self, path):
        """ .get_valid_names() """
        
        C = Catchment(os.path.join(path, 'test_multi_catchment.shp'))
        c_names = C.get_valid_names('test') # empty strings/white spaces are also 'None'
        assert c_names==['a', 'b'] # sorted, no nuplicates, no nan's
        
    def test_Catchment_get_valid_names2(self, path):
        """ Exception """
         
        C = Catchment(os.path.join(path, 'test_multi_catchment.shp'))
        with pytest.raises(ValueError) as e: # exception raised by .get_catchment_names()
            C.get_valid_names('bogus')
        assert e.value.args[0]==("Invalid catchment identifier 'bogus';\n" +
                                  "   attribute fields are: name, id, area_km2, test")
    
def test_Riversand_set_cid(path):

    rv = riversand.Riversand()
    rv.add_catchments('test_single_catchment.shp', path)
    assert rv.cid is None
    assert rv.catchments.cid is None
    
    # Catchment.set_cid() sets it for the catchment, not for the project rv
    rv.catchments.set_cid('id')
    assert rv.cid is None
    assert rv.catchments.cid=='id'
    
    # Riversand.set_cid() sets it for the project rv and for the catchment rv.catchments
    rv.set_cid('name')
    assert rv.cid=='name'
    assert rv.catchments.cid=='name'
    
    # Catchment.set_cid('bogus') RESETS it for the catchment, not leaves the project rv untouched
    rv.catchments.set_cid('bogus')
    assert rv.cid=='name'
    assert rv.catchments.cid is None
    
    rv.catchments.set_cid('name')
    rv.set_cid('id')
    assert rv.cid=='id'
    assert rv.catchments.cid=='id'

   
class Test_Riversand_add_samples():
    
    def test_Riversand_add_samples_Exception(self, path):
        """ wrapper for _from_dict() and _from_file() """
        
        rv = riversand.Riversand()
        
        with pytest.raises(TypeError) as e:
            rv.add_samples()
        assert e.value.args[0]=="add_samples() missing required keyword argument 'data' or 'fname'"
        assert rv.samples is None
        
    def test_Riversand_add_samples_Exception2(self, path):
        """ invalid kwarg 'file' instead of 'fname' """
        
        rv = riversand.Riversand()
        with pytest.raises(TypeError) as e:
            rv.add_samples(file='test_samples_invalid2.ods', path=path)
        assert e.value.args[0]=="add_samples() missing required keyword argument 'data' or 'fname'"


    def test_Riversand_add_samples_invalid1(self, path):
        """ invalid calls; prints error message, but no Exception is raised """

        rv = riversand.Riversand()
        
        rv.add_samples(data={'N': 10})
        assert rv.samples is None
        
        rv.add_samples(data={'N': 10, 'delN': 1, 'press_flag': 'pre'})
        assert rv.samples is None
        
    def test_Riversand_add_samples_valid1(self, path):
        """ valid calls; add data """
        
        rv = riversand.Riversand()
            
        rv.add_samples(data={'N': 10, 'delN': 1, 'press_flag': 'std'})
        assert len(rv.samples)==1
        
        rv.add_samples(data={'N': 10, 'delN': 1, 'comment': 'my comment'}, add=True)
        assert len(rv.samples)==2
        
        rv.add_samples(data={'N': 10, 'delN': 1, 'shielding': .99, 'new_col': 'something'})
        assert len(rv.samples)==1
        assert 'comment' not in rv.samples.keys()
        
    def test_Riversand_add_samples_invalid2(self, path):
        """ add data from dict; on error keep data untouched"""
        
        rv = riversand.Riversand()
        rv.add_samples(data={'N': 10, 'delN': 1, 'new_col': 'something'})
        assert len(rv.samples)==1
        
        rv.add_samples(data={'N': 10, 'delN': 1, 'density': -1}, add=True)
        assert len(rv.samples)==1
        assert 'new_col' in rv.samples.keys()
        
    def test_Riversand_add_samples_invalid3(self, path):
        """ add data from file; on error keep data untouched"""
        
        rv = riversand.Riversand()
        rv.add_samples(data={'N': 10, 'delN': 1, 'new_col': 'something'})
        assert len(rv.samples)==1
        
        rv.add_samples(fname='test.txt')
        assert len(rv.samples)==1
        assert 'new_col' in rv.samples.keys()

        rv.add_samples(fname='test_samples1.ods', path=path)
        assert len(rv.samples)==10
        
    def test_Riversand_add_samples_invalid4(self, path):
        """ save data to .samples even if mandatory col is missing """
        
        rv = riversand.Riversand()
        rv.add_samples(fname='test_samples_invalid1.ods', path=path)
        assert len(rv.samples)==9
        assert 'N' not in rv.samples.keys()
        
    def test_Riversand_add_samples_invalid5(self, path):
        """ empty columns in spreadsheet are removed """
        
        rv = riversand.Riversand()
        rv.add_samples(fname='test_samples_invalid2.ods', path=path)
        assert len(rv.samples)==7
        assert 'empty' not in rv.samples.keys()
        assert 'N' in rv.samples.keys()
        

    
    """ Add sample data from dict """
    
    def test_Riversand_add_samples_from_dict_exceptions(self):
        """invalid input data"""
        rv = riversand.Riversand()
        with pytest.raises(TypeError) as e:
            rv.add_samples_from_dict()
        assert e.value.args[0]=="add_samples_from_dict() missing keyword argument 'data'"
        
        df = pd.DataFrame({'N':[1e6, 1e6], 'delN': [143, 100]})
        with pytest.raises(NotImplementedError) as e:
            rv.add_samples_from_dict(data=df)
        
        with pytest.raises(TypeError) as e:
            rv.add_samples_from_dict(data=['N', 'delN'])
        
        with pytest.raises(ValueError) as e:
            rv.add_samples_from_dict(data={'N': 10})
        assert e.value.args[0]=="Invalid sample data: Uncertainty delN is not defined"
        
        with pytest.raises(ValueError) as e:
            rv.add_samples_from_dict(data={'N': 10, 'delN':1, 'density':0})
        assert e.value.args[0]=="Invalid sample data: Density is 0 or less"
        assert rv.samples is None
        
        rv.add_samples_from_dict(data={'N': 10 ,'delN': 1}) # create a valid entry
        assert len(rv.samples)==1
        rv.add_samples_from_dict(data={'N': 10, 'delN':1}, add=True) # add valid data
        assert len(rv.samples)==2
        with pytest.raises(ValueError) as e:
            rv.add_samples_from_dict(data={'N': 10}, add=True) # add invalid data
        assert len(rv.samples)==2, "with add=True existing data remains untouched"
    
        with pytest.raises(ValueError) as e:
            rv.add_samples_from_dict(data={'N': 10}) # invalid data, add=False
        assert rv.samples is None
        
    def test_Riversand_add_samples_from_dict_add(self):
        """option add=True"""
        # try add to non-existing rv.samples
        rv = riversand.Riversand()
        rv.add_samples_from_dict(data={'N': 10 ,'delN': 1}, add=True)
        assert len(rv.samples)==1, "override add=True"
        
        # try add to empty rv.samples
        rv = riversand.Riversand()
        rv.samples = pd.DataFrame()
        assert rv.samples is not None
        assert len(rv.samples)==0
        rv.add_samples_from_dict(data={'N': 10 ,'delN': 1}, add=True)
        assert len(rv.samples)==1, "can add to empty pd.DataFrame"
        
        # try add to invalid rv.samples
        rv.samples = 333
        rv.add_samples_from_dict(data={'N': 10 ,'delN': 1}, add=True)
        assert len(rv.samples)==1, "override add=True" 
        
    def test_Riversand_add_sample_from_dict_valid(self):
        """fill in default values if keys are not present"""
        rv = riversand.Riversand()
        rv.add_samples_from_dict(data={'N': 10 ,'delN': 1.2, # rounded and cast to int
                                 'comment': 'a comment',
                                 'long': 72., 'elevation': 0, 'lat': -33}) # exact
        assert rv.samples.equals(pd.DataFrame(
            {'name': ['Test'], 'press_flag': ['std'], 'thickness': [0],
             'density': [2.65], 'shielding': [1.], 'erate': [0], 'year': [2010],
             'nuclide': ['Be-10'], 'mineral': ['quartz'], 'N': [10], 'delN': [1],
             'standardization': ['07KNSTD'],
             'lat': [-33], 'long': [72.], 'elevation': [0],
             'comment': ['a comment']}))
        
    def test_Riversand_add_sample_from_dict_valid_add(self):
        """keep existing columns, but set default column ordering"""
        rv = riversand.Riversand()
        rv.samples = pd.DataFrame({'name': ['Test'], 'elevation': [0], 'mineral': ['quartz']})
        assert len(rv.samples)==1
        assert rv.samples.keys().equals(pd.Index(['name', 'elevation', 'mineral']))
        
        rv.add_samples_from_dict(data={'N':10, 'delN':1}, add=True)
        assert len(rv.samples)==2
        assert rv.samples.keys().equals(
            pd.Index(['name', 'press_flag', 'thickness', 'density', 'shielding', 
            'erate', 'year', 'nuclide', 'mineral', 'N', 'delN', 
            'standardization', 'elevation'])
            )
    
    def test_Riversand_add_sample_from_dict_bogus(self):
        """ignore bogus kwargs"""
        rv = riversand.Riversand()
        rv.add_samples_from_dict(data={'N':10, 'delN':1}, add=True, bogus=33)
        assert rv.samples.keys().equals(
            pd.Index(['name', 'press_flag', 'thickness', 'density', 'shielding', 
            'erate', 'year', 'nuclide', 'mineral', 'N', 'delN', 
            'standardization'])
            )


    """ Add sample data from file """
    
    def test_Riversand_add_sample_from_file_exceptions(self, path):
        """invalid input data"""
        rv = riversand.Riversand()
        with pytest.raises(TypeError) as e:
            rv.add_samples_from_file()
        assert e.value.args[0]=="add_samples_from_file() missing required keyword argument 'fname'"
        
        with pytest.raises(TypeError) as e:
            rv.add_samples_from_file(fname='test.txt', path=None)
        assert e.value.args[0]=="add_samples_from_file() 'path' must be string"

        with pytest.raises(FileNotFoundError) as e: # file doesnt exist in path=''
            rv.add_samples_from_file(fname='test.txt')
        assert e.value.args[0]=="test.txt: No such file or directory"
        
        # TODO: add test for IOError and ImportError
        
        with pytest.raises(ValueError) as e: # existing bogus file
            rv.add_samples_from_file(fname='test.txt', path=path)
        assert e.value.args[0]=="No sample data in file"
        
        with pytest.raises(ValueError) as e: # existing bogus file
            rv.add_samples_from_file(fname='test.tif', path=path)
        assert e.value.args[0]=="Excel file format cannot be determined, you must specify an engine manually."
                
    def test_Riversand_add_samples_from_file_invalid(self, path):
        """
        - set self.samples to sorted dataset; empty columns are not included
        - validate presence of columns name, N, delN
        - validate rows can be passed to validate_nuclide() and validate_sample()
          validations raise Exception but data is still stored in self.samples
        """
        rv = riversand.Riversand()
        # missing 'name', 'N'; empty column 'empty'
        with pytest.raises(ValueError) as e:
            rv.add_samples_from_file(fname='test_samples_invalid1.ods', path=path)
        assert e.value.args[0]=="Missing column(s): N", "Error detected by self.validate('sample')"
        assert rv.samples.keys().equals(
            pd.Index(['press_flag', 'density', 'shielding', 
            'year', 'nuclide', 'delN',
            'lat', 'long', 'sample', 'labID', 'comment'])
            )
        assert len(rv.samples)==9
        
        with pytest.raises(ValueError) as e:
            rv.add_samples_from_file(fname='test_samples_invalid2.ods', path=path)
        assert e.value.args[0]=="Invalid / missing data in lines: 3", "Error detected by self.validate('sample')"
        assert len(rv.samples)==7
        


class Test_Riversand_validate():
    
    def test_Riversand_validate_raster(self, path):
        """
        if 'raster' in args:
        - if self.elevation is not defined return ["No elevation raster defined"]
        - if all rasters have same epsg and res set self.epsg and self.res
          else return ["Conflicting projections in raster data"] and reset epsg, res
        """
        rv = riversand.Riversand()
        errors = rv.validate('raster', verbose=False)
        assert errors==["No elevation raster defined"]
        
        rv.add_raster('test.tif', path, 'shielding') 
        errors = rv.validate('raster', verbose=False)
        assert errors==["No elevation raster defined"]
        
        rv.add_raster('test.tif', path, 'elevation') # shielding = elevation
        errors = rv.validate('raster', verbose=False)
        assert errors==[]
        assert rv.epsg==32632
        assert rv.res==(500, 500)
        
        rv.add_raster('dem_utm_35m.tif', path, 'elevation') # different res
        assert rv.epsg is None, "Adding new raster should reset epsg and res"
        assert rv.res is None
        errors = rv.validate('raster', verbose=False)
        assert errors==["Conflicting projections in raster data"]
        assert rv.epsg is None
        assert rv.res is None
    
    
    def test_Riversand_validate_catchment_multi(self, path):
        """
        if 'catchment' in args:
        - validate projection:
        - C.dtype=='multi':
          - set self.valid_catchments if self.samples are defined and if there are
            valid catchment-sample pairs
          - self.valid_catchments==[] if samples defined but no valid data pairs    
        """
        rv = riversand.Riversand()
        errors = rv.validate('catchment', verbose=False)
        assert errors==["No catchment data defined"]
        
        rv.add_catchments('test_multi_catchment.shp', path) 
        errors = rv.validate('catchment', verbose=False)
        assert errors==["No elevation raster defined",
                        "Raster data invalid, cannot validate shapefile projection",
                        "Multi-catchment shapefile; use .set_cid() to set catchment identifier"]
        
        
        rv.add_raster('dem_WGS.tif', path, 'elevation')
        errors = rv.validate('catchment', verbose=False)
        assert errors==["Shapefile projection (epsg=32632) does not match raster projection (epsg=4326)",
                        "Multi-catchment shapefile; use .set_cid() to set catchment identifier"]
        
        rv.add_raster('dem_utm_35m.tif', path, 'elevation')
        errors = rv.validate('catchment', verbose=False)
        assert errors==["Multi-catchment shapefile; use .set_cid() to set catchment identifier"]
        
        rv.set_cid('name')
        assert rv.cid=='name'
        assert rv.catchments.cid=='name'
        errors = rv.validate('catchment', verbose=False)
        assert errors==[]
        assert rv.valid_catchments is None # no sample data defined
        
        rv.add_samples(fname='test_samples1.ods', path=path)
        assert isinstance(rv.samples, pd.DataFrame)
        errors = rv.validate(verbose=False)
        assert errors==[]
        assert rv.valid_catchments==['DB02', 'DB04', 'DB05', 'DB12', 'DB17', 'DB19']
        #DB03 is duplicate polygon
        #DB06, DB07 are sample data without polygon
        
        rv.set_cid('id')
        errors = rv.validate(verbose=False)
        assert errors==[]
        assert rv.valid_catchments==[], "cross-validated but no match found"
        
        rv.set_cid('name') #same for 'bogus', None
        assert rv.valid_catchments is None, "reset by self.set_cid()"
        
        
    def test_Riversand_validate_catchment_single(self, path):
        """
        if 'catchment' in args:
        - C.dtype=='single':
            - validate extent of all raster datasets
        """
        rv = riversand.Riversand()
        rv.add_raster('dem_utm_35m.tif', path, 'elevation')
        rv.add_catchments('test_single_catchment.shp', path) 
    
        rv.add_raster('toposhielding_35m_tooSmall.tif', path, 'shielding')    
        errors = rv.validate('catchment', verbose=False)
        assert errors==["Catchment polygon out of bounds of 1 raster dataset(s)"]
        
        rv.add_raster('quartz_35m_tooSmall.tif', path, 'quartz')
        errors = rv.validate('catchment', verbose=False)
        assert errors==["Catchment polygon out of bounds of 2 raster dataset(s)"]
        
        rv.add_raster('toposhielding_35m.tif', path, 'shielding')    
        errors = rv.validate('catchment', verbose=False)
        assert errors==["Catchment polygon out of bounds of 1 raster dataset(s)"]
        
        rv.quartz = None
        errors = rv.validate('catchment', verbose=False)
        assert errors==[]
    
    def test_Riversand_validate_kwargs(self, path):
        """ verbose=True / False; dtype='single' / 'multi'"""
        
        rv = riversand.Riversand()
        rv.add_raster('dem_utm_35m.tif', path, 'elevation')
        rv.add_catchments('test_single_catchment.shp', path, dtype='single')
        
        # valid, verbose=False
        errors = rv.validate('catchment', verbose=False)
        assert errors==[]
        # valid, verbose=True
        errors = rv.validate('catchment', verbose=True)
        assert errors is None

        # invalid, verbose=False
        rv.add_catchments('test_multi_catchment.shp', path, dtype='single')
        errors = rv.validate('catchment', verbose=False)
        assert errors==["Multi-catchment shapefile; use .set_cid() to set catchment identifier"]
        # invalid, verbose=True
        errors = rv.validate('catchment', verbose=True)
        assert errors is None
        
        # dtype='single'
        rv.add_catchments('test_multi_catchment.shp', path)
        with pytest.raises(ValueError) as e: 
            rv.validate('catchment', verbose=False, dtype='single')
        
        errors = rv.validate('catchment', verbose=False, dtype='multi')
        assert errors==["Multi-catchment shapefile; use .set_cid() to set catchment identifier"]
        
        rv.set_cid('name')
        errors = rv.validate('catchment', verbose=False, dtype='multi')
        assert errors==[]

        rv.add_catchments('test_single_catchment.shp', path, dtype='single')
        assert rv.catchments.dtype=='single'
        errors = rv.validate('catchment', verbose=False, dtype='multi')
        assert rv.catchments.dtype=='multi'
        
class Test_geospatial_functions():
    
    @pytest.fixture(autouse=True)
    def _RC(self, path):
        self._R = Raster(os.path.join(path, 'toposhielding_35m.tif'), 'shielding')
        self._C = Catchment(os.path.join(path, 'test_multi_catchment.shp'))
        
        self._cid = 'name'
        
        self._rv = riversand.Riversand()
        self._rv.add_raster(os.path.join(path, 'dem_utm_35m.tif'), dtype='elevation')
        self._rv.add_raster(os.path.join(path, 'toposhielding_35m.tif'), dtype='shielding')
        self._rv.add_raster(os.path.join(path, 'quartz_35m.tif'), dtype='quartz')
        
        self._rvF = riversand.Riversand() # faulty raster
        self._rvF.add_raster('dem_utm_35m.tif', path, 'elevation')
        self._rvF.add_raster('toposhielding_50m.tif', path, 'shielding')
        self._rvF.add_raster('quartz_35m.tif', path, 'quartz')
        
    
    def test_feature_in_raster_True(self):
        assert len(self._C.catchments)==8
        """n=3: DB05 in bounds"""
        n = 3
        assert self._C.catchments[n]['properties'][self._cid]=='DB05'
        polygon = self._C.catchments[n]['geometry']
        assert feature_in_raster(polygon, self._R.src)
    
    def test_feature_in_raster_False(self):
        """n=4: DB19 out of bounds"""
        n = 4
        assert self._C.catchments[n]['properties'][self._cid]=='DB19'
        polygon = self._C.catchments[n]['geometry']
        assert feature_in_raster(polygon, self._R.src)==False
        
    def test_feature_in_raster_Exception(self):
        """invalid polygon"""
        n = 4
        polygon = self._C.catchments[n]
        with pytest.raises(TypeError) as e:
            feature_in_raster(polygon, self._R.src)
        assert e.value.args[0]=="feature_in_raster() argument 'polygon' not a valid polygon dictionary"
        
    def test_clip_raster_1(self):
        """n=3: DB05 in bounds"""
        n = 3
        polygon = self._C.catchments[n]['geometry']
        Z = clip_raster(polygon, self._R.src, label='bla')
        assert isinstance(Z, xr.DataArray)
        assert Z.label=='bla'
        assert Z.shape==(598, 655)
        assert 1-(np.nanmin(Z.values)/0.58945)<1e-5 # min value ~0.58945
        
    def test_clip_raster_2(self, path):
        """different resolution raster"""
        n = 3
        R2 = Raster(os.path.join(path, 'toposhielding_50m.tif'), 'shielding')
        polygon = self._C.catchments[n]['geometry']
        Z = clip_raster(polygon, R2.src, label='bla')
        assert isinstance(Z, xr.DataArray)
        assert Z.label=='bla'
        assert Z.shape==(420, 458)
        assert 1-(np.nanmin(Z.values)/0.600668)<1e-5 # min value ~0.600668
        
    def test_clip_raster_Exception_out_of_bounds(self):
        """n=4: DB19 out of bounds"""
        n = 4
        polygon = self._C.catchments[n]['geometry']
        with pytest.raises(ValueError) as e:
            clip_raster(polygon, self._R.src, label='bla')
        assert e.value.args[0]=="clip_raster() : catchment polygon out of bounds"

    def test_clip_all_rasters_valid(self, path):
        """clip_all_rasters(n) validates raster and catchment data"""
        self._rv.add_catchments('test_multi_catchment.shp', path)
        self._rv.set_cid('name')
        n = 3
        clips = self._rv.clip_all_rasters(n)
        assert isinstance(clips, dict)
        assert set(clips.keys())=={'elevation', 'shielding', 'quartz', 'epsg'}
        assert clips['elevation'].shape==(598, 655)
    
    @mark.skip(reason="changed n=None, needs testing")
    def test_clip_all_rasters_Single(self, path):
        """Test the new feature n=None"""
        
    def test_clip_all_rasters_Exceptions(self, path):
        """exceptions raised by clip_all_rasters()"""
        n = 3 # in bound
        self._rv.add_catchments('test_multi_catchment.shp', path)
        with pytest.raises(ValueError) as e:
            self._rv.clip_all_rasters(n)
        assert e.value.args[0]=="Cannnot validate 'rv'; use 'rv.validate()' for details"
           
        assert self._rv.cid is None
        self._rv.add_catchments('test_single_catchment.shp', path, 'multi')
        assert self._rv.catchments.dtype=='multi'
        with pytest.raises(ValueError) as e:
            self._rv.clip_all_rasters(n)
        assert e.value.args[0]=="Cannnot validate 'rv'; use 'rv.validate()' for details"
        
        self._rv.add_catchments('test_single_catchment.shp', path)
        with pytest.raises(ValueError) as e:
            self._rv.clip_all_rasters(n)
        assert e.value.args[0]==("clip_all_rasters() : 'catchments' does not "+
                         "have n=3 polygons")
        
        n = 4 # out of bounds
        self._rv.add_catchments('test_multi_catchment.shp', path, 'multi')
        self._rv.set_cid('name')
        with pytest.raises(ValueError) as e:
            self._rv.clip_all_rasters(n)
        assert e.value.args[0]=="clip_all_rasters() : catchment polygon out of bounds"
        
    @mark.skip(reason="function modified; test needs to be adjusted")
    def test_eliminate_quartzfree(self, path):  
        """function eliminate_quartzfree(clips)"""
        self._rv.add_catchments('test_multi_catchment.shp', path)
        self._rv.set_cid('name')
        n = 3
        clips = self._rv.clip_all_rasters(n)
        outstr = eliminate_quartzfree(clips)
        assert outstr=="Removed 29.0 % of the catchment as quartz-free"
        
    def test_eliminate_quartzfree_Exception(self, path):  
        """no quartz raster defined"""
        self._rv.quartz=None
        self._rv.add_catchments('test_multi_catchment.shp', path)
        self._rv.set_cid('name')
        n = 3
        clips = self._rv.clip_all_rasters(n)
        with pytest.raises(ValueError) as e:
            eliminate_quartzfree(clips)
        assert e.value.args[0]=="eliminate_quartzfree() argument 'clips' missing required key 'quartz'"
        
    def test_get_xarray_centroid(self):
        """all data in utm"""
        n = 3
        polygon = self._C.catchments[n]['geometry']
        X = clip_raster(polygon, self._R.src, label='test')
        (x,y) = get_xarray_centroid(X)
        (x,y) = (int(np.round(x, 0)), int(np.round(y, 0)))
        assert (x,y)==(371169, 5050904) # confirmed in  QGIS; x, y

    def test_get_xarray_centroid_WGS(self, path):  
        """data in WGS84"""
        C = Catchment(os.path.join(path, 'DB05_WGS.shp'))
        R = Raster(os.path.join(path, 'dem_WGS.tif'), dtype='elevation')
        polygon = C.catchments[0]['geometry']
        X = clip_raster(polygon, R.src, label='test')
        (x,y) = get_xarray_centroid(X)
        (x,y) = (np.round(x, 3), np.round(y, 3))
        assert (x,y)==(7.348, 45.6) # confirmed in  QGIS; long, lat

    def test_get_polygon_centroid(self, path):  
        """data in WGS84"""
        C = Catchment(os.path.join(path, 'DB05_WGS.shp'))
        polygon = C.catchments[0]['geometry']
        (x,y) = get_polygon_centroid(polygon)
        (x,y) = (np.round(x, 3), np.round(y, 3))
        assert (x,y)==(7.348, 45.6) # confirmed in  QGIS; long, lat
        
    def test_get_polygon_centroid_Exception(self, path):
        """invalid polygon"""
        C = Catchment(os.path.join(path, 'DB05_WGS.shp'))
        polygon = C.catchments[0]#['geometry']        
        with pytest.raises(TypeError) as e:
            get_polygon_centroid(polygon)
        assert e.value.args[0]=="get_polygon_centroid() argument 'polygon' not a valid polygon dictionary"

    def test_projected_xy_to_longlat(self, path):  
        """projected_xy_to_longlat(xy, epsg)"""
        (x,y) = (371169, 5050904)
        epsg = self._rv.elevation.epsg
        (long, lat) = projected_xy_to_longlat((x,y), epsg)
        assert (np.round(long, 3), np.round(lat, 3))==(7.348, 45.6)


    def test_get_bins(self):  
        """function get_bins(Z, binsize=100)"""
        SF = clip_raster(self._C.catchments[0]['geometry'],
                         self._rv.shielding.src,
                         'bla')
        bins = get_bins(SF, .1)
        exp = np.array([0.5, 0.6, 0.7, 0.8, 0.9, 1. , 1.1])
        assert abs(np.sum(bins-exp))<1e-5
        assert np.nanmax(SF)<=1

    @mark.skip
    def test_get_topostats(self, path):  
        """
        # - epsg read from clips
        # - epsg specified
        # - centroid tuple
        # - centroid ='from_clipped"
        # - valid / invalid dict 'clips' # redundant
        # - clips has no elevation raster # redundant
        # - bins number / iterable
        # - centroid='from_polygon' is not yet implemented
        """
    

    
# =============================================================================
# Validate data for online calculator
# =============================================================================

class Test_validate_for_online_calculator():
    class Test_validate_topo():
        
        @mark.parametrize("item",
                          [({'lat': 32., 'long': 72.}), # missing mandatory key
                            ({'lat': 32., 'long': None, 'elevation': 1000}), # None
                            ({'lat': 32., 'long': 72., 'elevation': np.nan}), # np.nan
                            ({'lat': '', 'long': 72., 'elevation': 1000}), # empty string
                            ({'lat': 32., 'long': 'a', 'elevation': 1000}), # string
                            ({'lat': 32., 'long': [72.], 'elevation': 1000}), # strange data type
                          ])
        def test_validate_topo_exception(self, item):
            """validate_topo() with invalid data raises Exception"""
            with pytest.raises(ValueError) as e:
                validate_topo(item)
                
        def test_validate_topo_valid(self):
            """validate_topo() with valid data; does not alter input"""
            item = {'lat': 32., 'long': '-72', 'elevation': '0'}
            out = validate_topo(item)
            assert isinstance(out, dict)
            assert item['elevation']=='0'
            assert out['elevation']==0
            
        def test_validate_topo_Series(self):
            """validate_topo() accepts pd.Series"""
            item = pd.Series({'lat': 32., 'long': '-72', 'elevation': '0'})
            out = validate_topo(item)
            assert isinstance(out, dict)
            assert item['elevation']=='0'
            assert out['elevation']==0
            
        @mark.parametrize("shielding, expected",
                          [(None, "Shielding is not defined"),
                           (np.nan, "Shielding is not a number"),
                           ('', "Shielding is not a number"),
                           ('a', "Shielding is not a number"),
                           (1.1, "Shielding is out of bounds (0..1)")])
        def test_validate_topo_shielding_exception(self, shielding, expected):
            """validate_topo(shielding=True) with invalid data raises Exception"""
            
            item = {'lat': 32., 'long': 72., 'elevation': 1000}
            with pytest.raises(ValueError) as e:
                validate_topo(item, shielding=True)
            assert e.value.args[0] == "Shielding is not defined"
            
            item['shielding'] = shielding
            with pytest.raises(ValueError) as e:
                validate_topo(item, shielding=True)
            assert e.value.args[0]==expected
            
        def test_validate_topo_shielding_valid(self):
            """validate_topo(shielding=True/False)"""
            item = {'lat': 32., 'long': 72., 'elevation': 1000, 'shielding': '0.9'}
            assert validate_topo(item, shielding=True)['shielding']==0.9
            assert 'shielding' not in validate_topo(item).keys()
            
    class Test_validate_nuclide():
        
        @mark.parametrize("item",
                          [({'name': 'A', 'N': 10000}), # missing mandatory key
                           ({'name': 'A', 'N': 10000, 'delN': None}), # None
                           ({'name': 'A', 'N': 10000, 'delN': np.nan}), # np.nan
                           ({'name': 'A', 'N': 10000, 'delN': ''}), # empty string
                           ({'name': 'A', 'N': 10000, 'delN': 'a'}), # string
                           ({'name': 'A', 'N': 10000, 'delN': pd.DataFrame()}), # strange data type
                           ({'N':10, 'delN': 1, 'name': 'Ph.1'}),
                           ({'N':10, 'delN': 1, 'name': 'Ph/1'}),
                           ({'N':10, 'delN': 1, 'name': 'Ph\1'}),
                           ({'N':10, 'delN': 1, 'standardization': 'KNSTD'}),
                           ({'N':10, 'delN': 1, 'nuclide': 'Be-10', 'standardization': 'KNSTD'}),
                           ({'N':10, 'delN': 1, 'nuclide': 'Al-26', 'standardization': '07KNSTD'}),
                          ])
        def test_validate_nuclide_exception(self, item):
            """validate_nuclide() with invalid data raises Exception"""
            with pytest.raises(ValueError) as e:
                validate_nuclide(item)
                
        def test_validate_nuclide_valid(self):
            """validate_nuclide() with valid data; does not alter input"""
            item = {'name': 1, 'N': '1.2E+6', 'delN': 1.2E+4}
            out = validate_nuclide(item)
            assert isinstance(out, dict)
            assert item['N']=='1.2E+6'
            assert out['N']==1200000
            assert out['delN']==12000
            assert out['name']=='1'
            # default values
            assert out['nuclide']=='Be-10'
            assert out['mineral']=='quartz'
            assert out['standardization']=='07KNSTD'
            
        def test_validate_nuclide_nuclides(self):
            """validate_nuclide(); nuclides, mineral, standardization"""
            item = {'N': '1.2E+6', 'delN': 1.2E+4}
            out = validate_nuclide(item)
            assert out['nuclide']=='Be-10'
            assert out['mineral']=='quartz'
            assert out['standardization']=='07KNSTD'
            
            item = {'N': '1.2E+6', 'delN': 1.2E+4, 'nuclide': 'Al-26'}
            out = validate_nuclide(item)
            assert out['nuclide']=='Al-26'
            assert out['standardization']=='KNSTD'
            
            item['standardization'] = 'KNSTD'
            out = validate_nuclide(item)
            assert out['nuclide']=='Al-26'
            assert out['standardization']=='KNSTD'
            
            item['nuclide'] = 'Al-26'
            item['standardization'] = '07KNSTD'
            with pytest.raises(ValueError) as e:
                validate_nuclide(item)
            assert e.value.args[0]=="Illegal standardization"
            
            item['mineral'] = 'fsp'
            with pytest.raises(ValueError) as e:
                validate_nuclide(item)
            assert e.value.args[0]=="Illegal mineral"
            
        def test_validate_nuclide_Series(self):
            """validate_nuclide() accepts pd.Series"""
            item = pd.Series({'name': 'A', 'N': 10, 'delN': 1, 'nuclide': 'Al-26',
                              'comment': 'a comment'})
            out = validate_nuclide(item)
            assert 'comment' in item.keys()
            assert 'comment' not in out.keys()
            assert out['standardization']=='KNSTD'
           
        
    class Test_validate_sample():
        
        @mark.parametrize("item, expected",
                          [({'press_flag': '', 'density': 0, 'thickness': 0, 'erate': 0},
                            "Illegal pressure flag"),
                           ({'press_flag': 'std', 'density': None, 'thickness': 0, 'erate': 0},
                            "Density is not a number"),
                           ({'press_flag': 'std', 'density': np.nan, 'thickness': 0, 'erate': 0},
                            "Density is not a number"),
                           ({'press_flag': 'std', 'density': '', 'thickness': 0, 'erate': 0},
                            "Density is not a number"),
                           ({'press_flag': 'std', 'density': 2, 'thickness': 'a', 'erate': 0},
                            "Thickness is not a number"),
                           ({'press_flag': 'std', 'density': 2, 'thickness': pd.DataFrame(), 'erate': 0},
                            "Thickness is not a number"),
                           ({'press_flag': 'std', 'density': 0, 'thickness': 0, 'erate': 0},
                            "Density is 0 or less"), 
                           ({'press_flag': 'std', 'density': 2, 'thickness': 0, 'erate': -1},
                            "Erosion rate is less than 0"),
                          ])
        def test_validate_sample_exception(self, item, expected):
            """validate_sample() with invalid data raises Exception"""
            with pytest.raises(ValueError) as e:
                validate_sample(item)
            assert e.value.args[0]==expected
                
        def test_validate_sample_valid(self):
            """validate_sample() with valid data; does not alter input"""
            item = {'press_flag': 'ant', 'density': '1.8', 'thickness': 1e-1, 'erate': '8e-3'}
            out = validate_sample(item)
            assert item['density']=='1.8'
            assert item['thickness']==0.1
            assert item['erate']=='8e-3'
            assert out['density']==1.8
            assert out['thickness']==0.1
            assert out['erate']==0.008
            # default values
            assert out['year']==2010
        
        @mark.parametrize("year, expected",
                          [('1999.9', 2000),
                           (2001.1, 2001)
                          ])
        def test_validate_sample_year(self, year, expected):
            """validate_sample(); year"""
            out = validate_sample({'year': year})
            assert out['year']==expected
    
        def test_validate_sample_shielding_valid(self):
            """validate_sample(shielding=True/False)"""
            item = {'press_flag': 'ant', 'density': '1.8', 'thickness': 1e-1,
                    'shielding': 0.9, 'comment': 'some comment'}
            assert validate_sample(item, shielding=True)['shielding']==0.9
            assert 'shielding' not in validate_sample(item).keys()
            assert 'comment' not in validate_sample(item).keys()
 
    
def test_process_single_catchment(rv_single):
    """ handle RuntimeErrors raised by getE() """
    
    rv_single.add_samples(data={'N': 12000, 'delN': 400})
    rv_single.add_samples(data={'N': 2e+12, 'delN': 400, 'nuclide':'Al-26'}, add=True)
    rv_single.add_samples(data={'N': 12, 'delN': 400}, add=True)
    rv_single.add_samples(data={'N': 124000, 'delN': 400}, add=True)
    
    results = rv_single.process_single_catchment()
    expected = pd.Series({0: '',
                          1: 'sample appears to be saturated',
                          2: 'nuclide concentration too low for calculation',
                          3: ''})
    assert (results['error']==expected).all()
    
    expected = pd.Series({0: 'Be-10 quartz',
                          1: 'Al-26 quartz',
                          2: 'Be-10 quartz',
                          3: 'Be-10 quartz'})
    assert (results['nuclide']==expected).all()
    
def test_process_multi_catchment(rv_multi, path):
    """ handle RuntimeErrors raised by getE() """
    
    rv_multi.add_samples(fname='test_samples_multi.ods', path=path)
    rv_multi.set_cid('name')
    
    results = rv_multi.process_multi_catchment(unit='mm/kyr')
    expected = pd.Series({0: 'catchment out of bounds',
                          1: 'no catchment polygon', # in fact duplicate catchment
                          2: 'sample appears to be saturated',
                          3: 'nuclide concentration too low for calculation',
                          4: '',
                          5: 'no catchment polygon',
                          6: 'no catchment polygon',
                          7: '',
                          8: 'catchment out of bounds'
                         })
    assert (results['error']==expected).all()
    
    expected = pd.Series({0: 'Be-10 quartz',
                          1: 'Be-10 quartz',
                          2: 'Be-10 quartz',
                          3: 'Al-26 quartz',
                          4: 'Be-10 quartz',
                          5: 'Al-26 quartz',
                          6: 'Be-10 quartz',
                          7: 'Be-10 quartz',
                          8: 'Be-10 quartz'
                         })
    assert (results['nuclide']==expected).all()
