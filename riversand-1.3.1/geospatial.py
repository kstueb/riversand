#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

*******************************************************************************
geospatial.py  :  geospatial processing for the riversand package 

    Copyright (C) 2023  Konstanze Stübner, kstueb@gmail.com

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

*******************************************************************************

This module implements the 'Riversand' object as well as the 'Raster' and
'Catchment' objects.

A 'Riversand' object has attributes:
    - Raster objects .elevation, .shielding, .quartz (.add_raster())
    - a Catchment object .catchments (.add_catchment())
    - a pandas DataFrame .samples (.add_samples())
For multi-catchment datasets, the catchment identifier is set with .set_cid()
A list of catchment names is obtained with .catchments.get_name()
The project is validated with .validate()

Calculations are performed with:
results = rv.process_single_catchment()
results = rv.process_multi_catchment()
results = rv.catchment_stats()
    
The current version of the online calculator is obtained from calc.py
version = calc.get_version()

"""

import os

import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype

import collections
from numbers import Number

import rasterio
import rasterio.crs
import rasterio.mask

from pyproj.crs import CRS
import pyproj.exceptions

import fiona

from riversand.params import Params

import warnings
warnings.filterwarnings("error") # convert warnings to raise exceptions


class OutOfBoundsError(Exception):
    # used by clip_all_rasters() and caught by process_single/multi_catchment()
    pass

# =============================================================================
# Raster datasets
# =============================================================================

def get_geotiff(fname:str) -> rasterio.DatasetReader:
    """ Get file handle to geotiff. """

    src = None
    mn = None
    mx = None
    if not os.path.isfile(fname):
        raise FileNotFoundError("{}: No such file or directory".format(fname))
    
    try:
        src = rasterio.open(fname, 'r')
        mn = np.nanmin(src.read())
        mx = np.nanmax(src.read())
        src.close()
    except rasterio.errors.NotGeoreferencedWarning as e:
        raise pyproj.exceptions.CRSError(e.args[0]) # 'Dataset has no geotransform, gcps, or rpcs. The identity matrix will be returned.'
    except:# rasterio.RasterioIOError as error:
        pass
    
    if not src is None:
        src.min = mn
        src.max = mx
    return src


class Raster():
    """
    Elevation, shielding or lithology raster dataset.
    
    Attributes are:
    - self.fname: file name
    - self.src: source (closed DatasetReader)
    - self.crs: coordinate reference system
    - self.epsg: epsg code
    - self.res: resolution (tuple)
    - self.dtype: type ('elevation', 'shielding' or 'quartz')
    
    """
    
    dtypes = ('elevation', 'shielding', 'quartz') # valid Raster types
    
    def __init__(self, fname:str, dtype:str):
        
        dtype = dtype.lower()
        if dtype not in self.dtypes:
            raise TypeError("Invalid Raster dtype '{}' (must be 'elevation', 'shielding' or 'quartz')".format(dtype))
        
        if not isinstance(fname, str):
            raise TypeError("Invalid Raster fname (must be string)")
            
        src = get_geotiff(fname) # closed rasterio.DatasetReader
        # raises except pyproj.exceptions.CRSError  
                
        if src:
            self.dtype = dtype
            self.fname = fname
            self.src = src
            self.crs = CRS(src.crs)
            self.epsg = self.crs.to_epsg() # may be None
            self.res = src.res # resolution, tuple
            self.min = src.min
            self.max = src.max
            
            
            if self.crs.is_projected==False:
                raise pyproj.exceptions.CRSError(
                    ["Geotiff: no projected coordinate reference system",
                     self.crs])
            if self.epsg is None:
                raise pyproj.exceptions.CRSError(
                    ["Geotiff: cannot convert crs to epsg",
                     self.crs])
                
        else:
            raise IOError("Cannot read geotiff {}".format(fname))
        
    def __repr__(self):
       
        s = []
        s += ["dtype  : {}".format(self.dtype)]
        s += ["fname  : {}".format(self.fname)]
        s += ["src    : {}".format(self.src)]
        s += ["epsg   : {}".format(self.epsg)]
        s += ["res    : {}".format(self.res)]
        s += ["values : {} to {}".format(self.min, self.max)]
        return "\n".join(s)

# =============================================================================
# Catchment shapefiles
# =============================================================================

def get_shapefile(fname:str) -> fiona.Collection:
    """ Get file handle to shapefile. """

    src = None
    
    if not os.path.isfile(fname):
        raise FileNotFoundError("{}: No such file or directory".format(fname))
    
    try:
        with fiona.Env(OSR_WKT_FORMAT="WKT2_2018"):
            src = fiona.open(fname)
        src.close()
    except fiona.errors.DriverError as e:
        print(e)
            
    return src


class Catchment():
    """
    Attributes are:
    - self.fname: file name
    - self.src: source (closed Collection)
    - self.crs: coordinate reference system
    - self.epsg: epsg code
    - self.attrs: list of shapefile attribute fields
    - self.cid: attribute field used to match samples to catchments
    - self.catchments: list of catchment polygons (fiona.model.Feature)
    - self.len: number of catchment polygons in shapefile
    
    """
    
    def __init__(self, fname:str):
        
        if not isinstance(fname, str):
            raise TypeError("Invalid Catchment fname (must be string)")
        
        src = get_shapefile(fname) # closed fiona.collection
        if src is not None:
                    
            with fiona.open(fname) as src:
                schema = src.schema
                attrs = src.schema['properties'].keys()
                catchments = [f for f in src]
            
                geometry = schema['geometry']
                if geometry!='Polygon':
                    raise ValueError("Invalid shapefile geometry '{}' (must be 'Polygon')"
                                     .format(geometry))
    
                if len(catchments)==0:
                    raise ValueError("Shapefile empty, must have at least 1 polygon")
                    
                crs = CRS(src.crs)
                
                
            self.fname = fname
            self.src = src
            
            self.attrs = list(attrs)
            self.crs = crs
            self.epsg = crs.to_epsg() # may be None
            self.catchments = catchments
            
            self.cid = None
            
            if self.crs.is_projected==False:
                raise pyproj.exceptions.CRSError(
                    ["Shapefile: no projected coordinate reference system",
                     self.crs])
            if self.epsg is None:
                raise pyproj.exceptions.CRSError(
                    ["Shapefile: cannot convert crs to epsg",
                     self.crs])

        else:
            raise IOError("Cannot read shapefile {}".format(fname))
        
        
    def set_cid(self, cid:str='id'):
        """
        Set Catchment.cid (fyi only, Riversand.cid is relevant).
        
        Parameters
        ----------
        cid : str
            Shapefile attribute field used for catchment names.
            
        """
        
        if cid is None:
            self.cid = None
        else:
            if self.attrs is None:
                print("ValueError : Invalid cid='{}';\n".format(cid) +
                      "   shapefile has no attribute fields")
                self.cid = None
            if cid in self.attrs:
                self.cid = cid
            else:
                print("ValueError : Invalid cid='{}';\n".format(cid) +
                      "   shapefile attribute fields are: {}".format(", ".join(self.attrs)))
                self.cid = None


    def get_names(self, cid:str=None) -> list:
        """
        Get sorted list of all catchment names from attribute field 'cid'.
        Duplicates are included, missing values are shown as 'None'.
        
        Parameters
        ----------
        cid : str
            Shapefile attribute field used for catchment names.
            
        """
        
        if cid is None:
            cid = self.cid
        if cid is None:
            print("ValueError : use .set_cid() to set catchment identifier")
            return None
        
        if cid not in self.attrs:
            print("ValueError : Invalid cid='{}';\n".format(cid) +
                  "   shapefile attribute fields are: {}".format(", ".join(self.attrs)))
            return None
         
        try:
            c_names = [str(c['properties'][cid]) for c in self.catchments]
        except Exception as exc:
            raise exc
           
        c_names.sort()
        return c_names
        
    def get_valid_names(self, cid:str) -> list:
        """
        Get a sorted list of valid catchment names from attribute field 'cid'.
        Duplicates and missing values are excluded.
        
        Parameters
        ----------
        cid : str
            Shapefile attribute field used for catchment names.
        
        """
        
        c_names = self.get_names(cid) # returns None for invalid cid
        if c_names is None:
            print("ValueError : use .set_cid() to set catchment identifier")
        else:
            c_invalid = [itm for itm, cnt in collections.Counter(c_names).items() if cnt>1] # non-unique
            c_invalid += ['None', '']
            c_names = [c for c in c_names if c not in c_invalid]        
            c_names.sort()
        return c_names
        

    def __len__(self):
        """ Number of catchment polygons. """
        return len(self.catchments)
        
    def __repr__(self):
        s = []
        #s += ["dtype : {}".format(self.dtype)]
        s += ["fname  : {}".format(self.fname)]
        s += ["src    : {}".format(self.src)]
        
        s += ["attrs  : {}".format(self.attrs)]
        s += ["len    : {}".format(len(self.catchments))]
        s += ["epsg   : {}".format(self.epsg)]
        
        if self.cid:
            s += ["cid    : {}".format(self.cid)]
        
        return "\n".join(s)


# =============================================================================
# Riversand object
# =============================================================================
 

class Riversand():
    """
    A 'Riversand' object contains all the data needed to calculate
    catchmentwide erosion rates. Methods are:
    
    self.set_path_to_data() # define path to input data
    self.add_raster()  # add Raster objects (elevation, shielding, quartz)
    self.add_catchments()  # add Catchment object (catchments)
    self.add_samples()  # add sample data (samples)
 
    self.validate()   # validate projection and resolution of geospatial data
        (res, crs, epsg.)
        
    self.process_single_catchment()  # calculate erosion rate for a single catchment
    self.process_multi_catchment()  # calculate erosion rates for a shapefile with
        multiple catchments (use self.set_cid() to set the shapefile attribute
        to be used as catchment name)
        
    self.catchment_stats()  # calculate some catchment statistics
    self.restandardize()  # convert nuclide concentrations to defaults
        Be-10: 07KNSTD; Al-26: KNSTD
        
    """
    

    
    def __init__(self, path:str=None):
        self.path_to_data = None
        self.elevation = None # Raster object
        self.shielding = None # Raster object
        #self.snow = None # Raster object
        self.quartz = None # Raster object
        
        self.res = None # project resolution
        self.crs = None # project projection
        self.epsg = None # project projection
        self.cid = None # catchment identifier for multi-catchment datsets
        
        self.catchments = None # Catchment object

        self.samples = None # pd.DataFrame of sample data
        
        self.valid_catchments = None # set by self.validate()
        # not validated if within raster bounds or if sample data is valid
        # set to None by .add_catchments(), .set_cid(),
        #     .add_samples_from_file(), .add_samples_from_dict()
            
        if path is not None:
            try:
                self.set_path_to_data(path)
            except FileNotFoundError as e:
                print(type(e).__name__, ":", e.args[0])
        
    def set_path_to_data(self, path:str=None):
        """ Set path to input data. """    
        
        if os.path.isdir(path):
            self.path_to_data = path
        else:
            raise FileNotFoundError("{}: No such directory".format(path))
        
        
    def add_raster(self, fname:str, dtype:str='elevation'):
        """ Add raster dataset to the project. """
        
        if isinstance(fname, Raster):
            raise NotImplementedError("add_raster() not yet implemented for Raster object")
            
        if self.path_to_data is not None:
            fname = os.path.join(self.path_to_data, fname)
            
        dtype = dtype.lower()
        try:
            if dtype=='elevation':
                self.elevation = Raster(fname, dtype)
            elif dtype=='shielding':
                self.shielding = Raster(fname, dtype)
                if not self.elevation is None:
                    if self.shielding.epsg != self.elevation.epsg:
                        print("WARNING: Shielding raster projection (epsg) not matching DEM")
                    if self.shielding.res != self.elevation.res:
                        print("WARNING: Shielding raster resolution not matching DEM")
                    if (self.shielding.min is None or self.shielding.min <0 or
                        self.shielding.max is None or self.shielding.max >1):
                        print("WARNING: Shielding raster data outside of valid range (0..1)")
                    
            elif dtype=='quartz':
                self.quartz = Raster(fname, dtype)
                if not self.elevation is None:
                    if self.quartz.epsg != self.elevation.epsg:
                        print("WARNING: Quartz raster projection (epsg) not matching DEM")
                    if self.quartz.res != self.elevation.res:
                        print("WARNING: Quartz raster resolution not matching DEM")
                    if (self.quartz.min is None or self.quartz.min <0 or
                        self.quartz.max is None or self.quartz.max >1):
                        print("WARNING: Quartz raster data outside of valid range (0..1)")
            else:
                Raster(fname, dtype) # raises exceptions
        except FileNotFoundError as e:
                print("{} : {}: No such file"
                      .format(type(e).__name__, e.args[0].split(':')[0]))
        except pyproj.exceptions.CRSError as e:
            if isinstance(e.args[0], list):
                print(type(e).__name__, ":",
                      "Geotiff coordinate reference system ", end='')
                if e.args[0][1].is_projected==False:
                    print("is not a projected reference system.")
                if e.args[0][1].to_epsg() is None:
                    print("cannot be converted no EPSG code.")
                print("")
                print(e.args[0][1].__repr__())
            else:
                # "Dataset has no geotransform..." Warning raised as Exception by get_geotiff()
                print(type(e).__name__, ":", e.args[0])
        
        self.crs = None
        self.epsg = None
        self.res = None
        
    def add_catchments(self, fname:str):
        """ Add catchment shapefile to the project. """
        
        if isinstance(fname, Catchment):
            raise NotImplementedError("add_catchments() not yet implemented for Catchment object")
            
        if self.path_to_data is not None:
            fname = os.path.join(self.path_to_data, fname)
            
        try:
            self.catchments = Catchment(fname)
        except FileNotFoundError as e:
                print("{} : {}: No such file"
                      .format(type(e).__name__, e.args[0].split(':')[0]))
        except pyproj.exceptions.CRSError as e:
            print("Shapefile coordinate reference system ", end='')
            if e.args[0][1].is_projected==False:
                print("is not a projected reference system.")
            if e.args[0][1].to_epsg() is None:
                print("cannot be converted no EPSG code.")
            print("")
            print(e.args[0][1].__repr__())
        
        else: #only if catchment was added (no Error)
            if not self.elevation is None:
                if self.catchments.crs.to_epsg() != self.elevation.epsg:
                    print("WARNING: Catchment projection (epsg) not matching DEM")
                
        self.crs = None
        self.epsg = None
        self.res = None
        
        self.cid = None
        self.valid_catchments = None # set by self.validate()
        
        
    def set_cid(self, cid:str='id'):
        """ Set attribute 'cid'. """
        
        if self.catchments is None:
            print("ValueError : No catchments defined")
        else:
            if self.catchments.attrs is None:
                print("ValueError : Invalid cid='{}';\n".format(cid) +
                      "   shapefile has no attribute fields")
                self.cid = None
                self.catchments.set_cid(None)
            if cid in self.catchments.attrs:
                self.cid = cid
                self.catchments.set_cid(cid)
            else:
                print("ValueError : Invalid cid='{}';\n".format(cid) +
                      "   shapefile attribute fields are: {}"
                      .format(", ".join(self.catchments.attrs)))
                self.cid = None
                self.catchments.set_cid(None)
            
        self.valid_catchments = None # set by self.validate()
        
        
    def add_samples(self, samples, add:bool=False):
        """
        Add to or replace self.samples with sample data from dict, pd.Series,
        pd.DataFrame or spreadsheet file.
        
        Data from file or DataFrame cannot be added to existing data. Data is
        valildated and errors are printed but data are still stored in
        self.samples if possible.
        
        Prints error messages returned from self.add_samples_from_xxx().
        
        """
        
        if isinstance(samples, dict):
            try:
                self.add_samples_from_dict(samples, add)
            except Exception as e: # non-fatal exceptions printed as error message
                print("Error adding sample data from dictionary:\n   "+str(e))
            
        elif isinstance(samples, pd.Series):
            try:
                self.add_samples_from_dict(samples, add) # .add_samples_from_dict() accepts dict or pd.Series
            except Exception as e: # non-fatal exceptions printed as error message
                print("Error adding sample data from pandas Series:\n   "+str(e))

        elif isinstance(samples, pd.DataFrame):
            try:
                self.add_samples_from_dataframe(samples) # no adding to existing data
            except Exception as e: # non-fatal exceptions printed as error message
                print("Error adding sample data from pandas DataFrame:\n   "+str(e))
            
        elif isinstance(samples, str):
            try:
                self.add_samples_from_file(samples) # no adding to existing data
            except FileNotFoundError as e:
                    print("{} : {}: No such file"
                          .format(type(e).__name__, e.args[0].split(':')[0]))
            except Exception as e: # non-fatal exceptions printed as error message
                print("Error adding sample data from file:\n   "+str(e))
        else:
            raise TypeError("add_samples() missing mandatory argument "+
                            "(file name or dictionary)")

                
    def add_samples_from_file(self, fname:str):
        """
        Write sample data from spreadsheet to self.samples.
        Colunms sorted as in Params.all_keys with topo and other info at the end of the table.
        
        Raises ValueErrors "\n".join(errors) returned by:
        errR, errS, errC = self.validate(verbose=False)
        
        """
                            
        if self.path_to_data is not None:
            fname = os.path.join(self.path_to_data, fname)
            
        df = pd.DataFrame()
        
        try:
            df = pd.read_csv(fname, sep=None, engine='python') # python engine automatically detects separator
            
        except:
            try:
                df = pd.read_excel(fname) # imports the first sheet no matter the name
            except:
                #print("Failing to read input as csv or excel file.")
                pass
        
        if df.empty:
            raise ValueError("Failing to read file / no sample data in file")
        
        #strip leading/trailing spaces from all strings
        df_obj = df.select_dtypes('object')
        df[df_obj.columns] = df_obj.apply(lambda x: x.str.strip())
        
        # sorting as in Params.all_keys w/o topographic and other info at end 
        topo_keys = ['lat', 'latitude', 'lon', 'long', 'longitude', 'elev', 'elevation']
        mandatory_keys = [k for k in Params.all_keys if k not in topo_keys]
        other_keys = [k for k in topo_keys if k in df.keys()] # additional topographic info
        other_keys += [k for k in df.keys() if k not in Params.all_keys] # additional other
        
        set_default = [k for k in mandatory_keys if k not in df.keys()] # missing from data; set to default values 
        
        self.samples = pd.DataFrame(df, columns=mandatory_keys+other_keys) # sorted columns

        # set default values for missing columns 
        for d in set_default:
            if d in Params.default_values.keys(): # excluding standardization
                self.samples[d] = Params.default_values[d]
        if 'standardization' in set_default:
            self.samples['standardization'] = self.samples['standardization'].astype(str)
            self.samples.loc[self.samples['nuclide']=='Be-10', 'standardization']=Params._default_standards['Be-10']
            self.samples.loc[self.samples['nuclide']=='Al-26', 'standardization']=Params._default_standards['Al-26']

        self.samples.dropna(axis=1, how='all', inplace=True) # remove columns that are all nan
        
        if len(set_default)>0:
            print("Using default values for missing columns:\n   {}".format(', '.join(set_default)))
            
        self.valid_catchments = None # set by self.validate()

        _, errS, _ = self.validate(verbose=False)
        if errS:
            raise ValueError(errS[0]) # e.g. 'Invalid / missing data in line(s): 6, 7'


    def add_samples_from_dataframe(self, data:pd.DataFrame):
        """
        Write sample data from pd.DataFrame to self.samples.
        Colunms sorted as in Params.all_keys with topo and other info at the end of the table.
        
        Raises ValueErrors "\n".join(errors) returned by:
        errR, errS, errC = self.validate(verbose=False)
    
        """
        
        if not isinstance(data, (dict, pd.DataFrame)):
            raise TypeError("add_samples_from_dataframe() 'data' must be pandas DataFrame")
    
        if data.empty:
            raise ValueError("No sample data in DataFrame")
        
        #strip leading/trailing spaces from all strings
        df_obj = data.select_dtypes('object')
        data[df_obj.columns] = df_obj.apply(lambda x: x.str.strip())
        
        # sorting as in Params.all_keys w/o topographic and other info at end 
        topo_keys = ['lat', 'latitude', 'lon', 'long', 'longitude', 'elev', 'elevation']
        mandatory_keys = [k for k in Params.all_keys if k not in topo_keys]
        other_keys = [k for k in topo_keys if k in data.keys()] # additional topographic info
        other_keys += [k for k in data.keys() if k not in Params.all_keys] # additional other
    
        set_default = [k for k in mandatory_keys if k not in data.keys()] # missing from data; set to default values 
    
        self.samples = pd.DataFrame(data, columns=mandatory_keys+other_keys) # sorted columns
       
        # set default values for missing columns 
        for d in set_default:
            if d in Params.default_values.keys(): # excluding standardization
                self.samples[d] = Params.default_values[d]
        if 'standardization' in set_default:
            self.samples.loc[self.samples['nuclide']=='Be-10', 'standardization']=Params._default_standards['Be-10']
            self.samples.loc[self.samples['nuclide']=='Al-26', 'standardization']=Params._default_standards['Al-26']
    
        self.samples.dropna(axis=1, how='all', inplace=True) # remove columns that are all nan
    
        if len(set_default)>0:
            print("Using default values for missing columns:\n   {}".format(', '.join(set_default)))
        
        self.valid_catchments = None # set by self.validate()
        
        _, errS, _ = self.validate(verbose=False)
        if errS:
            raise ValueError(errS[0]) # e.g. 'Invalid / missing data in line(s): 6, 7'
    
    
    def add_samples_from_dict(self, data:dict, add:bool=False):
        """
        Write sample data from dict or pd.Series to self.samples.
        Colunms sorted as in Params.all_keys with topo and other info at the
        end of the table.
        
        """
        
        from riversand.utils import validate_nuclide, validate_sample
                        
        if not isinstance(data, (dict, pd.Series)):
            raise TypeError("add_samples_from_dict() 'data' must be dictionary")
        
        if add==False:
            self.samples = None
                
        if add==True and (
                (self.samples is None) or
                (not isinstance(self.samples, pd.DataFrame))):
            add = False # can only add to existing data
        
        topo_keys = ['lat', 'latitude', 'lon', 'long', 'longitude', 'elev', 'elevation']
        mandatory_keys = [k for k in Params.all_keys if k not in topo_keys]  
        
        if 'shielding' not in data.keys(): #default shielding=1 if not defined 
            data['shielding'] = 1.
        
        exc = None
        try:
            nuclide_dict = validate_nuclide(data)
            sample_dict = validate_sample(data, shielding=True)
        except ValueError as exc:
            # reports error string exc returned by validate_nuclide(), validate_sample() (first error encountered if multiple)
            raise ValueError("Invalid sample data: {}".format(exc)) from exc
        else:
            # wrap all non-mandatory data into a dict starting with the topograpic data
            other_keys = [k for k in topo_keys if k in data.keys()] # additional topographic info
            other_keys += [k for k in data.keys() if k not in Params.all_keys] # additional other
            other_dict = dict((k, data[k]) for k in other_keys)
        
            
            if add==False: # new dataset or adding to empty Dataframe
                self.samples = pd.DataFrame(nuclide_dict|sample_dict|other_dict,
                                            columns=mandatory_keys+other_keys,
                                            index=[0])       
            else:
                self.samples.reset_index(inplace=True, drop=True)
                # add new data to the end but use column ordering off added (validated) data
                df = pd.DataFrame(nuclide_dict|sample_dict|other_dict,
                                  columns=mandatory_keys+other_keys,
                                  index=[len(self.samples)])
                self.samples = (pd.concat([df, self.samples])
                                .sort_index().reset_index(drop=True))        
                
            self.valid_catchments = None # set by self.validate()
        
    def restandardize(self, verbose=True):
        """
        Convert self.samples to Be-10: 07KNSTD; Al-26: KNSTD.
        Requires columns 'nuclide', 'standardization', 'N', 'delN'.
        
        Note that utils.restandardize() is not caught up by samples with
        invalid standardization and/or N=0; it sets 'N'='delN'='standardization'
        to NaN). self.validate() will call out these problems.

        """
        from riversand import utils
        
        if self.samples is None:
            raise ValueError("No samples defined, use .add_samples() before re-standardizing")
        
        self.samples = utils.restandardize(self.samples)
        #raises KeyError if N, delN, nuclide, standardization are missing; this should never happen for self.samples
        if verbose:
            print("Nuclide data converted to default standardizations Be-10: {}; Al-26: {}."
                  .format(Params._default_standards['Be-10'],
                          Params._default_standards['Al-26']))
        
        
    
    def get_valid_catchments(self):
        self.validate(verbose=False)
        if self.valid_catchments is None:
            return []
        else:
            return self.valid_catchments
    
        
    def validate(self, multi:bool=None, verbose=True):
        """ Validate Riversand object and print error report. """
        
        from riversand.utils import validate_nuclide, validate_sample
        from riversand.utils import feature_in_raster

        class ValidationError(Exception):
            pass
        
        def val_R(self):
            # adding a raster resets epsg and res to None anyways
            self.epsg = None
            self.res = None
            
            if self.elevation is None:
                raise ValidationError("No elevation raster defined")
                
            if (self.elevation.epsg is None) or (self.elevation.res is None):
                raise ValidationError("Cannot determine projection or resolution of elevation raster")
            
            proj_errs = 0
            if self.shielding is not None:
                proj_errs += int(self.shielding.epsg!=self.elevation.epsg)
                proj_errs += int(self.shielding.res!=self.elevation.res)
                #print("checking shielding "+str(proj_errs))
                #if (self.shielding.min is None or self.shielding.min <0 or
                #    self.shielding.max is None or self.shielding.max >1):
                #    raise ValidationError("Shielding raster data outside of valid range (0..1)")
                    
            #if self.snow is not None:
            #    proj_errs += int(self.snow.epsg!=self.elevation.epsg)
            #    proj_errs += int(self.snow.res!=self.elevation.res)
            #    #print("checking quartz "+str(proj_errs))
            if self.quartz is not None:
                proj_errs += int(self.quartz.epsg!=self.elevation.epsg)
                proj_errs += int(self.quartz.res!=self.elevation.res)
                #print("checking quartz "+str(proj_errs))
                #if (self.quartz.min is None or self.quartz.min <0 or
                #    self.quartz.max is None or self.quartz.max >1):
                #    raise ValidationError("Quartz raster data outside of valid range (0..1)")
            if proj_errs>0:
                raise ValidationError("Conflicting projections in raster data")
            
            self.epsg = self.elevation.epsg
            self.crs = self.elevation.crs
            self.res = self.elevation.res
            return
        
        def val_S(self):
            if (not isinstance(self.samples, pd.DataFrame)
                or self.samples.empty):
                raise ValidationError("No sample data defined")
            
            tmp = self.samples#.reset_index(drop=True) # do not modify self.samples 
        
            # (1) Validate that mandatory columns 'N','delN' exist
            missing_columns = []
            for mc in ['N', 'delN']:
                if mc not in tmp.keys():
                    missing_columns += [mc]
            if len(missing_columns)>0:
                raise ValidationError("Missing column(s): {}"
                                      .format(", ".join(missing_columns)))
                
            # (2) Check row-wise that validation for online calculator will pass (only if no missing columns)
            errors_in_line = []
            for idx, item in tmp.iterrows():
                this_line_error = False
                try:
                    _ = validate_nuclide(item)
                except:
                    this_line_error = True
                try:
                    _ = validate_sample(item)
                except:
                    this_line_error = True
                if 'shielding' in tmp.keys():
                    try:
                        sf = float(item['shielding'])
                        if not np.isnan(sf): # nan are allowed in shielding
                            if not (0 <= sf <= 1): 
                                this_line_error = True 
                    except:
                        this_line_error = True
                if this_line_error:
                    errors_in_line += [str(idx)]
            if errors_in_line:
                raise ValidationError("Invalid / missing data in line(s): {}"
                                      .format(", ".join(errors_in_line)))
            return
        
        def val_C(self, outR, multi=False): # outR is the outcome of the raster validation; should be []
            self.valid_catchments = None
            
            if (not isinstance(self.catchments, Catchment)
                or self.catchments is None):
                raise ValidationError("No catchment data defined")
            
            if not outR==[]: # outR=[err mssg] if invalid or undefined ["No elevation raster defined"]
                raise ValidationError("Raster data undefined or invalid, cannot validate shapefile projection")
        
            C_epsg = self.catchments.crs.to_epsg()
            if C_epsg!=self.epsg:
                temp = self.epsg # for display
                self.epsg = None
                self.crs = None
                raise ValidationError("Shapefile projection (epsg={}) "
                                      .format(C_epsg) +
                                      "does not match raster projection (epsg={})"
                                      .format(temp))
            
            nC = len(self.catchments.catchments)
            if nC==0:
                raise ValidationError("No polygons in shapefile")
            
            # single-catchment dataset
            if multi==False:
                #print('single')
                if nC>1:
                    raise ValidationError("Not a single-catchment dataset ({} polygons)"
                                          .format(nC))
            
                out_of_bounds = 0
                rr = [r.src for r in [self.elevation, self.shielding, self.quartz]
                     if r is not None]
                for r in rr:
                    try:
                        fir = feature_in_raster(
                            self.catchments.catchments[0]['geometry'], r)
                        if not fir:
                            out_of_bounds += 1
                    except:
                        pass
                if out_of_bounds>0:
                    raise ValidationError("Catchment polygon out of bounds of raster data")
                
            # multi-catchment dataset
            if multi==True:
                #print('multi')
                if self.cid is None:
                    raise ValidationError("No catchment identifier defined; use .set_cid()")
                if self.cid not in self.catchments.attrs: # self.set_cid() does not allow to set an invalid cid
                    raise ValidationError("Invalid catchment identifier cid='{}'; use .set_cid()"
                                          .format(self.cid))
                
                c_names = (self.catchments.get_valid_names(self.cid))
                try:
                    s_names = list(self.samples['name'].values)
                except:
                    s_names = [] # empty list if no column 'name'
                self.valid_catchments = [c for c in c_names if c in s_names]
            return
                
                
        # determine parameter 'multi' from the number of catchment polygons
        if multi is None:
            try:
                nC = len(self.catchments.catchments)
            except:
                multi = False # no catchment defined; doesn't matter
            else:
                if nC==1:
                    multi = False
                else:
                    multi = True
                    
        outR = outS = outC = None
        try:
            val_R(self)
        except ValidationError as e:
            outR = [e.args[0]]
        else:
            outR = []
        
        try:
            val_S(self)
        except ValidationError as e:
            outS = [e.args[0]]
        else:
            outS = []
        
        try:
            val_C(self, outR, multi=multi)
        except ValidationError as e:
            outC = [e.args[0]]
        else:
            outC = []
                    
        if verbose: # outX should never be None; [] or ["error message"]
            print("")
            if outR==[]:
                print("Raster data valid")
            else:
                #print("\nErrors in raster data:")
                print("{}".format(outR[0]))
            if outS==[]:
                print("Sample data valid")
            else:
                #print("\nErrors in sample data:")
                print("{}".format(outS[0]))
            if outC==[]:
                print("Catchment data valid")
            else:
                #print("\nErrors in catchment data:")
                print("{}".format(outC[0]))
            if self.valid_catchments is not None:
                print("\nValid catchments / samples:")
                if len(self.valid_catchments)==0:
                    print("   No matches found")
                else:
                    print("   Found {} match(es)".format(len(self.valid_catchments)))

        else:
            return outR, outS, outC
        
    
    def process_single_catchment(self,
            bins=500, scaling='LSDn', shielding=1, unit='mm/yr',
            plot=None,
            url=None, verbose=True) -> pd.DataFrame:
        """
        Calculate catchmentwide erosion rates in cm/yr for a single-catchment
        shapefile. Nuclide concentrations in self.samples are restandardized.

        Parameters
        ----------
        bins : float, optional
            Bin size in metres for elevation binning. The default is 500.
        scaling : str, optional
            Scaling method 'St', 'Lm' or 'LSDn'. The default is 'LSDn'.
        shielding : str or float, optional
            Shielding method, 'sample', 'topo' or numeric (0 to 1).
            The default is 1.
        unit : str, optional
            Unit for plotting. The default is 'mm/yr'.
        plot : str, optional
            Extension for saving plots, e.g. 'jpg', 'png', 'svg', 'eps'.
            The default is None.

        Returns
        -------
        results : pd.DataFrame
            Table of results for each sample.
            
        Raises
        ------
        ValueError
            self.elevation, self.catchments or self.samples are not defined.
            More than one polygon in self.catchments.
            Invalid argument 'shielding' or 'unit'.

        """
        
        if url is None: url = Params.url
            
        #from riversand.utils import validate_topo, validate_nuclide, validate_sample
        from riversand.utils import get_topostats
        from riversand.calc import get_textline, get_E, guess_erates
        from riversand.calc import poly_E_results
        import riversand.plot
        
        if scaling not in Params._scalingmethods:
            raise ValueError("Invalid scaling '{}': must be 'St', 'Lm' or 'LSDn'")
        
        if self.elevation is None:
            raise ValueError("Missing elevation raster")
        if self.samples is None:
            raise ValueError("Missing sample data")           
        if self.catchments is None:
            raise ValueError("Missing catchment polygons") 
            
        if len(self.catchments.catchments)>1:
            raise ValueError("Use .process_multi_catchment() if shapefile has more than one polygon")
            
        if isinstance(shielding, str):
            if shielding not in {'topo', 'sample'}:
                raise ValueError("Invalid shielding: must be 'topo', 'sample' or numeric")
        if isinstance(shielding, Number):
            if not (0<=shielding<=1):
                raise ValueError("Invalid shielding: must be 0..1")
        if (shielding=='topo') and (self.shielding is None):
            raise ValueError("Invalid shielding='topo' but no shielding raster available")
        if (shielding=='sample') and ('shielding' not in self.samples.keys()):
            raise ValueError("Invalid shielding='sample' but no shielding in sample data")
            
        if isinstance(plot, str):
            if plot[0]=='.': plot = plot[1:]
            if str.lower(plot) not in ('eps', 'jpeg', 'jpg', 'pdf', 'pgf', 'png',
                            'ps', 'raw', 'rgba', 'svg', 'svgz', 'tif', 'tiff', 'webp'):
                plot = None
        else:
            plot = None
        
        if not plot is None:
            try:
                os.makedirs(Params.out_path)
            except FileExistsError:
                pass
            except PermissionError:
                print("No permission to make directory for saving plots")
                plot = None
            
            
        if unit not in Params.units.keys():
            raise ValueError("Invalid unit='{}': see riversand.print_units()"
                             .format(unit))
            
            
        result_cols = ['name', 'scaling', 'nuclide',
                       'E', 'delE-', 'delE+', 'NRMSE', 'Tavg', 'error']
        results = pd.DataFrame(columns=result_cols)
        
        errR, errS, errC = self.validate(verbose=False, multi=False)
        if len(errR+errS+errC)>0:
            print("Cannot validate dataset for single-catchment processing; "+
                             "use '.validate()' for details")
            return results
        
        if verbose:
            print("Processing single catchment")
            if isinstance(bins, Number):
                print("Bin size : {} m".format(bins))
            else:
                print("Custom bins")
            print("Scaling method : {}".format(scaling))
            if shielding in {'sample', 'topo'}:
                print("Topographic shielding from {} data".format(shielding))
            elif isinstance(shielding, Number):
                print("Topographic shielding : {}".format(shielding))
            if not self.quartz is None:
                print("Correcting for quartz-free lithologies")
            if plot is None:
                print("Not saving plots")
            else:
                print("Saving plots as .{} in '{}'".format(plot, Params.out_path))
            print("")

        try:
            clips = self.clip_all_rasters() # function does not call self.validate()
        except OutOfBoundsError:
            #    "clip_all_rasters() : catchment polygon out of bounds"
            print("Catchment is out of bounds")
            return results
        except ValueError as e:
            #"non-matching projection (epsg) of raster datasets"
            #"shapefile projection (epsg) not matching raster datasets"
            #"{} raster outside of valid range (0..1)".format(r)
            err_code = e.args[0]
            print("ERROR: {}\n".format(err_code[:1].upper() + err_code[1:]))
            # the first two of these should never happen
        except RuntimeError:
            #RuntimeError #can only happen if argument n in clip_all_rasters()
            #    "clip_all_rasters() : .catchments does not have n={} polygons".format(n)
            print("Indexing error; this should never happen")
            return results
            
        else:
                           
            if not plot is None:
                for label in clips.keys():
                    if label in Raster.dtypes:
                        fig, ax = riversand.plot.plot_clipped_raster(
                            clips, dtype=label)
                        if not fig is None:
                            try:
                                fig.savefig(os.path.join(Params.out_path, label +'.'+ plot))
                            except PermissionError:
                                print("No permission to save plot")
                                plot = None
                        
            topostats, summary = get_topostats(clips, bins=bins, centroid='from_clipped') # accepts iterable as bins
            
            if len(topostats)==0:
                err = ['no quartz in catchment']
                print("No quartz in catchment")
                return results
                
            else:
            
                try: # mock textline to make sure that 'shielding' is set correctly
                    summary['elevation'] = 100
                    _ = get_textline(self.samples.iloc[0],
                                     summary,
                                     shielding=shielding)
                except ValueError as e: # raised by get_textline() for invalid shielding
                    print(type(e).__name__, ":", e.args[0])
                    return results
                    
                self.restandardize(verbose=False)
                # modifies self.samples to restandardized values
                # in either case poly_E_results() uses restandardized values for calculation but does not touch self.samples
                
                for idx, sample_data in self.samples.iterrows(): # sample_data as pd.Series
                    
                    E = np.nan
                    delE = (np.nan, np.nan)
                    NRMSE = np.nan
                
                    # sample name from input or default value
                    try:
                        name = str(sample_data['name'])
                    except:
                        name = Params.default_values['name']
                    
                    # nuclide and mineral from input or default value
                    try:
                        nuclide = sample_data['nuclide']+' '
                    except:
                        nuclide = Params.default_values['nuclide']+' '
                    try:
                        nuclide += sample_data['mineral']
                    except:
                        nuclide += Params.default_values['mineral']
                                        
                    try:
                        # estimate minimum and maximum erosion rates
                        summary['elevation'] = summary['elevLo']
                        textline = get_textline(sample_data, summary, shielding=shielding)
                        E_Lo = get_E(textline) # erosion rates in cm/yr
                            
                        summary['elevation'] = summary['elevHi']
                        textline = get_textline(sample_data, summary, shielding=shielding)
                        E_Hi = get_E(textline) # erosion rates in cm/yr
                    except ValueError as e: # raised by get_textline() for invalid shielding
                        err = [e.args[0]]
                    except RuntimeError as e: # raised by get_E()
                        err = [e.args[0][10:]]
                        
                    else:
                    
                        # generate suitable erates
                        erates = guess_erates(E_Lo, E_Hi, scaling=scaling) # cm/yr
                        
                        E, delE, NofE, RMSE, err = (
                            poly_E_results(sample_data, topostats, 
                            shielding=shielding, erates=erates, scaling=scaling, url=url))
                        
                        # catch out-of-bound errors
                        while 'minE too high' in err:
                            erates = guess_erates(3*erates[0]-2*erates[1], erates[-1]) # cm/yr
                            E, delE, NofE, RMSE, err = (
                                poly_E_results(sample_data, topostats, 
                                shielding=shielding, erates=erates, scaling=scaling, url=url))
                            
                        while 'maxE too low' in err:
                            erates = guess_erates(erates[0], 3*erates[-1]-2*erates[-2]) # cm/yr
                            E, delE, NofE, RMSE, err = (
                                poly_E_results(sample_data, topostats, 
                                shielding=shielding, erates=erates, scaling=scaling, url=url))
                        
                        NRMSE = RMSE/sample_data['N']
                        
                        if not plot is None:
                            if 'linear fit' in err:
                                fig, ax = riversand.plot.plot_polyfit(
                                    E, delE, NofE, sample_data, unit=unit,
                                    linfit=True)
                            else:
                                fig, ax = riversand.plot.plot_polyfit(
                                    E, delE, NofE, sample_data, unit=unit)
                            if not fig is None:
                                try:
                                    fig.savefig(os.path.join(
                                        Params.out_path,
                                        str(idx)+'_'+name+'_'+scaling+'.'+plot))
                                except PermissionError:
                                    print("No permission to save plot")
                                    #plot = None
                            
                            
                    
                    finally:
                        results.loc[idx, 'name'] = name    
                        results.loc[idx, 'scaling'] = scaling
                        results.loc[idx, 'nuclide'] = nuclide
                        
                        results.loc[idx, 'E'] = E
                        results.loc[idx, 'delE-'] = delE[0]
                        results.loc[idx, 'delE+'] = delE[1]
                        results.loc[idx, 'NRMSE'] = NRMSE
                        
                        if np.isnan(E):
                            results.loc[idx, 'Tavg'] = E
                        else:
                            Tavg = 160/sample_data['density']/E
                            #f = int(np.log10(Tavg))-2
                            #Tavg = np.round(Tavg/10**f,)*10**f
                            results.loc[idx, 'Tavg'] = int(Tavg)
                        
                        if err==[]:
                            results.loc[idx, 'error'] = ''
                        else:
                            results.loc[idx, 'error'] = "; ".join(err)
                        ## possible error messages:
                        # 'NRMSE = ...'
                        # 'Root finding did not converge'
                        # 'Cannot compute uncertainty'
                        # 'Server cannot resolve erosion rates, dropping duplicates'
        
        if verbose:
            print("Processing finished.")
        
        return results
    
        
    def process_multi_catchment(self,
            bins=500, scaling='LSDn', shielding=1, unit='mm/yr',
            plot=None,
            url=None, verbose=True) -> pd.DataFrame:
        """
        Calculate catchmentwide erosion rates in cm/yr for a multi-catchment
        shapefile. Nuclide concentrations in self.samples are restandardized.
        
        Parameters
        ----------
        bins : float, optional
            Bin size in metres for elevation binning. The default is 500.
        scaling : str, optional
            Scaling method 'St', 'Lm' or 'LSDn'. The default is 'LSDn'.
        shielding : str or numeric, optional
            Shielding method, 'sample', 'topo' or numeric (0 to 1).
            The default is 1.
        unit : str, optional
            Unit for plotting. The default is 'mm/yr'.
        plot : str, optional
            Extension for saving plots, e.g. 'jpg', 'png', 'svg', 'eps'.
            The default is None.
            
        Returns
        -------
        results : pd.DataFrame
            Table of results for each sample.
            
        Raises
        ------
        ValueError
            self.elevation, self.catchments or self.samples are not defined.
            Invalid argument 'shielding' or 'unit'.
            
        """
        
        if url is None: url = Params.url
            
        #from riversand.utils import validate_topo, validate_nuclide, validate_sample
        from riversand.utils import get_topostats
        from riversand.calc import get_textline, get_E, guess_erates
        from riversand.calc import poly_E_results
        import riversand.plot
        
        if scaling not in Params._scalingmethods:
            raise ValueError("Invalid scaling='{}': must be 'St', 'Lm' or 'LSDn'")
            
        if self.elevation is None:
            raise ValueError("Missing elevation raster")
        if self.samples is None:
            raise ValueError("Missing sample data")           
        if self.catchments is None:
            raise ValueError("Missing catchment polygons")            
            
        if isinstance(shielding, str):
            if shielding not in {'topo', 'sample'}:
                raise ValueError("Invalid shielding: must be 'topo', 'sample' or numeric")
        if isinstance(shielding, Number):
            if not (0<=shielding<=1):
                raise ValueError("Invalid shielding: must be 0..1")
        if (shielding=='topo') and (self.shielding is None):
            raise ValueError("Invalid shielding='topo' but no shielding raster available")
        if (shielding=='sample') and ('shielding' not in self.samples.keys()):
            raise ValueError("Invalid shielding='sample' but no shielding in sample data")

        if isinstance(plot, str):
            if plot[0]=='.': plot = plot[1:]
            if str.lower(plot) not in ('eps', 'jpeg', 'jpg', 'pdf', 'pgf', 'png',
                            'ps', 'raw', 'rgba', 'svg', 'svgz', 'tif', 'tiff', 'webp'):
                plot = None
        else:
            plot = None
        
        if not plot is None:
            try:
                os.makedirs(Params.out_path)
            except FileExistsError:
                pass
            except PermissionError:
                print("No permission to make directory for saving plots")
                plot = None    
       
        if not plot is None:
            for k,v in {'elevation': self.elevation,
                        'shielding': self.shielding,
                        'quartz': self.quartz 
                       }.items():
                if not v is None:
                    try:
                        os.makedirs(os.path.join(Params.out_path, k))
                    except FileExistsError:
                        pass
                    except PermissionError:
                        print("No permission to make subdirectory for saving plots")
                        plot = None

                

        if unit not in Params.units.keys():
            raise ValueError("Invalid unit='{}': see riversand.print_units()"
                             .format(unit))
            
            
        result_cols = ['name', 'scaling', 'nuclide',
                       'E', 'delE-', 'delE+', 'NRMSE', 'Tavg', 'error']
        results = pd.DataFrame(columns=result_cols)
        results['name'] = self.samples['name']

        errR, errS, errC = self.validate(verbose=False, multi=True)
        if len(errR+errS+errC)>0:
            print("Cannot validate dataset for multi-catchment processing; "+
                  "use '.validate()' for details")
            return results
        
        if verbose:
            print("Processing multi-catchment dataset")
            if isinstance(bins, Number):
                print("Bin size : {} m".format(bins))
            else:
                print("Custom bins")
            print("Scaling method : {}".format(scaling))
            if shielding in {'sample', 'topo'}:
                print("Topographic shielding from {} data".format(shielding))
            elif isinstance(shielding, Number):
                print("Topographic shielding : {}".format(shielding))
            if not self.quartz is None:
                print("Correcting for quartz-free lithologies")
            if plot is None:
                print("Not saving plots")
            else:
                print("Saving plots as .{} in '{}'".format(plot, Params.out_path))
            print("")
            print("")

        # some printout formatting
        If = len(str(len(self.samples))) # fill length for id
        Nf = max([len(l) for l in self.samples['name'].values]) # fill length for sample name
        #Ef = 26 # fill length for erosion rate result        
        order_str = "{: >"+str(If)+"} {: <"+str(Nf)+"} : "
        #order_str_E2 = "{: >"+str(If)+"} {: <"+str(Nf)+"} : {: >"+str(Ef)+"}"  # erosion rate result
        #order_str_F1 = "{: >"+str(If)+"} {: <"+str(Nf)+"} : "  # erosion rate result
        #order_str_F2 = "{: >"+str(If)+"} {: <"+str(Nf)+"} : {}"  # erosion rate result

        self.restandardize(verbose=False)
        # modifies self.samples to restandardized values
        # in either case poly_E_results() uses restandardized values for calculation but does not touch self.samples
        
        for idx, sample_data in self.samples.iterrows(): # iterate over all samples
        
            E = np.nan
            delE = (np.nan, np.nan)
            NRMSE = np.nan
            error_code = None

            name = str(sample_data['name'])
            # nuclide and mineral from input or default value
            try:
                nuclide = sample_data['nuclide']+' '
            except:
                nuclide = Params.default_values['nuclide']+' '
            try:
                nuclide += sample_data['mineral']
            except:
                nuclide += Params.default_values['mineral']
                
            if verbose: # print index and sample name
                print(order_str.format(idx, name), end='')
            
            
            if sample_data['name'] not in self.valid_catchments: #set by self.validate()
                error_code = 'no catchment polygon'
            else:
                # get integer id of the catchment
                n = [n for n, c in enumerate(self.catchments.catchments) 
                     if str(c['properties'][self.cid])==name][0]
                
                try:
                    clips = self.clip_all_rasters(n) # function does not call rv.validate()
                    
                except OutOfBoundsError:
                    #    "clip_all_rasters() : catchment polygon out of bounds"
                    error_code = 'catchment out of bounds'
                except ValueError as e:
                    #"non-matching projection (epsg) of raster datasets"
                    #"shapefile projection (epsg) not matching raster datasets"
                    #"{} raster outside of valid range (0..1)".format(r)
                    error_code = e.args[0]
                    # the first two of these should never happen                    
                except RuntimeError:
                    #RuntimeError
                    #    "clip_all_rasters() : .catchments does not have n={} polygons".format(n)
                    error_code = 'indexing error; use process_single_catchment()' # this should not happen
                    
                else:
                    
                    if plot is not None:
                        for label in clips.keys():
                            if label in Raster.dtypes:
                                fig, ax = riversand.plot.plot_clipped_raster(
                                    clips, dtype=label)
                                if not fig is None:
                                    try:
                                        fig.savefig(os.path.join(
                                            Params.out_path,
                                            label, #subfolder
                                            str(sample_data.name)+'_'+name+'.'+plot))
                                    except PermissionError:
                                        print("No permission to save plot")
                                        #plot = None

                    topostats, summary = get_topostats(clips, bins=bins, centroid='from_clipped') # accepts iterable as bins
                    
                    if len(topostats)==0:
                        err = ['no quartz in catchment']
                        print("no quartz in catchment")
                        
                    else:
                        try:
                            # estimate minimum and maximum erosion rates
                            summary['elevation'] = summary['elevLo']
                            textline = get_textline(sample_data, summary, shielding=shielding)
                            E_Lo = get_E(textline) # erosion rates in cm/yr
                            
                            summary['elevation'] = summary['elevHi']
                            textline = get_textline(sample_data, summary, shielding=shielding)
                            E_Hi = get_E(textline) # erosion rates in cm/yr
                        except ValueError as e: # raised by get_textline() for invalid shielding
                            print(type(e).__name__, ":", e.args[0])
                            err = [e.args[0]]
                        
                        except RuntimeError as e: # raised by get_E()
                            error_code = e.args[0][10:] # e.g. "get_E() : sample appears to be saturated"
                        else:
                
                            # generate suitable erates
                            erates = guess_erates(E_Lo, E_Hi, scaling=scaling) # cm/yr
                            E, delE, NofE, RMSE, err = (
                                poly_E_results(sample_data, topostats, 
                                shielding=shielding, erates=erates, scaling=scaling, url=url))
                    
                            # catch out-of-bound errors
                            while 'minE too high' in err:
                                erates = guess_erates(3*erates[0]-2*erates[1], erates[-1]) # cm/yr
                                E, delE, NofE, RMSE, err = (
                                    poly_E_results(sample_data, topostats, 
                                    shielding=shielding, erates=erates, scaling=scaling, url=url))
                                
                            while 'maxE too low' in err:
                                erates = guess_erates(erates[0], 3*erates[-1]-2*erates[-2]) # cm/yr
                                E, delE, NofE, RMSE, err = (
                                    poly_E_results(sample_data, topostats, 
                                    shielding=shielding, erates=erates, scaling=scaling, url=url))
                            
                            NRMSE = RMSE/sample_data['N']
    
                            # printout formatting
                            Estr = "{:.1f}+/-{:.1f} {}".format(E*Params.units[unit],
                                                            delE[1]*Params.units[unit],
                                                            unit)
                            
                            if not plot is None:
                                if 'linear fit' in err:
                                    fig, ax = riversand.plot.plot_polyfit(
                                        E, delE, NofE, sample_data, unit=unit,
                                        linfit=True)
                                else:
                                    fig, ax = riversand.plot.plot_polyfit(
                                        E, delE, NofE, sample_data, unit=unit)
                                if not fig is None:
                                    try:
                                        fig.savefig(os.path.join(
                                            Params.out_path,
                                            str(sample_data.name)+'_'+name+'_'+scaling+'.'+plot))
                                    except PermissionError:
                                        print("No permission to save plot")
                                        #plot = None
                                        
                            if verbose: # print erosion rate result
                                print(Estr)
                            
                       
                        
            results.loc[idx, 'name'] = sample_data['name']   
            results.loc[idx, 'scaling'] = scaling
            results.loc[idx, 'nuclide'] = nuclide
                        
            results.loc[idx, 'E'] = E
            results.loc[idx, 'delE-'] = delE[0]
            results.loc[idx, 'delE+'] = delE[1]
            results.loc[idx, 'NRMSE'] = NRMSE
            
            if np.isnan(E):
                results.loc[idx, 'Tavg'] = E
            else:
                Tavg = 160/sample_data['density']/E
                #f = int(np.log10(Tavg))-2
                #Tavg = np.round(Tavg/10**f,)*10**f
                results.loc[idx, 'Tavg'] = int(Tavg)    
            
                
            if error_code:
                if verbose:
                    print(error_code)
                    
                results.loc[idx, 'error'] = error_code
                # 'catchment out of bounds'
                # 'no catchment polygon'
                # ...saturated...
                # ...concenration too low...
            else:
                if err==[]:
                    results.loc[idx, 'error'] = ''
                else:
                    results.loc[idx, 'error'] = "; ".join(err)
                ## possible error messages:
                # 'NRMSE = ...'
                # 'Root finding did not converge'
                # 'Cannot compute uncertainty'
                # 'linear fit'
                # 'Server cannot resolve erosion rates, dropping duplicates'
        
        #try: # if all went wrong there might not be a variable clips
        #     if 'quartz' not in clips.keys():
        #         results.drop(columns=['qtz'], inplace=True)
        #except:
        #    pass
        
        if verbose:
            print("\nProcessing finished.")
        return results
    
    
    def catchment_stats(self, verbose=True) -> pd.DataFrame:
        """
        Calculate topo statistics for all catchments in shapefile.
    
        Returns
        -------
        results : pd.DataFrame
            Table of results for each catchment. Keys are:
            - centr_lat, centr_long : centroid coordinates (in WGS84)
            - mean_elev, stdev_elev : mean and standard deviation of catchment elevations
            - median_elev : median elevation
            - relief : difference betweem highesta and lowest elevation
            - area : catchment area in km2 assuming geotiff resolution in metres!
            - mean_sf : mean shielding factor
            - error : "catchment out of bounds" if not within bounds of geotiff
            
        Raises
        ------
        ValueError
            self.elevation, self.catchments or self.samples are not defined.
            
        """
        
        from riversand.utils import projected_xy_to_longlat, get_xarray_centroid
        
        if self.elevation is None:
            raise ValueError("Missing elevation raster")
        if self.samples is None:
            raise ValueError("Missing sample data")           
        if self.catchments is None:
            raise ValueError("Missing catchment polygons")
            
        topo_cols = ['name', 'centr_lat', 'centr_long', 'mean_elev', 'stdev_elev',
                     'median_elev', 
                     'relief', 'area', 'mean_sf']
        results = pd.DataFrame(columns=topo_cols)
            
        errR, errS, errC = self.validate(verbose=False) # 'multi' determined from number of catchments
        if len(errR+errS+errC)>0:
            print("Cannot validate dataset; "+
                  "use '.validate()' for details")
            return results

        if verbose:
            print("Processing ", end="")
            end = ""
            
        for idx, c in enumerate(self.catchments.catchments):
            if self.cid:
                name = str(c['properties'][self.cid])
            else:
                name = ""#c['id']
            results.loc[idx, 'name'] = name
            if verbose:
                print(end+name, end="")
                end = ", "
                
            try:
                clips = self.clip_all_rasters(idx) # function does not call rv.validate()
            except OutOfBoundsError:
                results.loc[idx, 'error'] = 'catchment out of bounds'
            except ValueError:
                results.loc[idx, 'error'] = 'indexing error; this should never happen' # this should not happen
            else:
                results.loc[idx, 'error'] = ''
                
                try:
                    (cx, cy) = get_xarray_centroid(clips['elevation'])
                except RuntimeError: # raised by get_xarray_centroid() for empty xarray
                    (cx, cy) = (np.nan, np.nan)
                    
                (cx, cy) = projected_xy_to_longlat((cx, cy), clips['epsg'])
                
                area_per_pixel = float(
                    np.abs(clips['elevation'].transform[0] * 
                           clips['elevation'].transform[4])) # in km2
                area = np.count_nonzero(~np.isnan(clips['elevation']))*area_per_pixel*1e-6
                
                results.loc[idx, 'centr_lat'] = np.round(cy, 5)
                results.loc[idx, 'centr_long'] = np.round(cx, 5)
                
                if area > 0: # avoid RuntimeWarning e.args[0]=='Mean of empty slice';  
                    results.loc[idx, 'mean_elev'] = np.round(
                        float(np.nanmean(clips['elevation'])), 1)
                    results.loc[idx, 'stdev_elev'] = np.round(
                        float(np.nanstd(clips['elevation'])), 1)
                    results.loc[idx, 'median_elev'] = np.round(
                        float(np.nanmedian(clips['elevation'])), 1)
                    results.loc[idx, 'relief'] = np.round(
                        np.nanmax(clips['elevation'])-np.nanmin(clips['elevation']), 1)
                    if 'shielding' in clips.keys():
                        results.loc[idx, 'mean_sf'] = np.round(
                            np.nanmean(clips['shielding']), 5)
                else:
                    results.loc[idx, 'mean_elev'] = np.nan
                    results.loc[idx, 'stdev_elev'] = np.nan
                    results.loc[idx, 'median_elev'] = np.nan
                    results.loc[idx, 'relief'] = np.nan
                    if 'shielding' in clips.keys():
                        results.loc[idx, 'mean_sf'] = np.nan
                    
                results.loc[idx, 'area'] = np.round(area, 1)
                
        
        if verbose:
            print(" finished.")
        
        results = results.sort_values(by=['name'], ignore_index=True)
        
        return results
    
    
    def __repr__(self):
        s = []
        s += ["---------------"]
        if self.path_to_data is None:
            s += ["Input data folder: not set"]
        else:
            s += ["Input data folder: '{}'".format(self.path_to_data)]
        s += ["Output folder for plots: '{}'\n".format(Params.out_path)]
        
        if self.elevation or self.shielding or self.quartz:
            s += ["---------------"]
            s += ["Raster data:\n"]
        if self.elevation:
            s += [str(self.elevation)+"\n"]
        if self.shielding:
            s += [str(self.shielding)+"\n"]
        if self.quartz:
            s += [str(self.quartz)+"\n"]   
            
        if not self.samples is None:
            s += ["---------------"]
            s += ["Sample data:"]
            s += ["{:d} sample(s)\n".format(len(self.samples))]
            
        if not self.catchments is None:
            s += ["---------------"]
            s += ["Catchment polygons:\n"]
            s += [str(self.catchments)+"\n"]
            
        if self.epsg or self.res:
            s += ["==============="]
            s += ["Validated projection:"]
            if self.epsg:
                s += ["epsg   : {}".format(self.epsg)]
            if self.res:
                s += ["res    : {}".format(self.res)]
            s += [""]
            
        if self.valid_catchments is None:
            pass
        elif self.valid_catchments==[]:
            s += ["---------------"]
            s += ["Catchments / samples cross-validated:"]
            s += ["no matches found"]
        else:
            s += ["---------------"]
            s += ["Catchments / samples cross-validated:"]
            s += ["Found {} matches".format(len(self.valid_catchments))]
            s += [", ".join(self.valid_catchments)]

        return "\n".join(s)
    


    
    def clip_all_rasters(self, n:int=0) -> dict:
        """
        Clip all rasters 'elevation', 'shielding', 'quartz' to one polygon.
        
        n : int
            Number of the catchment.
            > polygon = self.catchment.catchments[n]['geometry']
            For single-catchment shapefile set n=0 or or n=None (default value).
        
        Returns
        -------
        clips : dict
            Keys are 'elevation', 'shielding', 'quartz' with xarrays of the
            corresponding clipped raster. Additional keys 'epsg', 'name'.
        
        Exceptions
        ----------
        RuntimeError
            "polygon indexing error (arg n)"
        ValueError
            "non-matching projection (epsg) of raster datasets"
            "shapefile projection (epsg) not matching raster datasets"
            "{} values are nodata".format(r))
            "{} values are outside of valid range (0..1)".format(r)
        OutOfBoundsError 
            "catchment polygon out of bounds"
            
        """
        
        
        from riversand.utils import clip_raster
                
        #if self.epsg is None:
        #    print("WARNING : dataset does not seem validated, run .validate()")
        #errR, errS, errC = self.validate(verbose=False)
        #if len(errR+errC)>0:
        #    raise ValueError("Cannnot validate data; use .validate() for details")
        
        try:
            polygon = self.catchments.catchments[n]['geometry'] # raises IndexError
            #try:
            #    c_name = str(self.catchments.catchments[n]['properties'][self.cid])
            #except:
            #    c_name = '' # for Exception message
        except IndexError:
            raise ValueError("polygon indexing error (arg n)"
                             ) from None
        try: # catchment name
            name = str(self.catchments.catchments[n].properties[self.cid])
        except:
            name = ''
            
        clips = {}
        
        rr = [r for r in [self.elevation, self.shielding, self.quartz]
             if r is not None]
        
        epsgs = set(r.epsg for r in [self.elevation, self.shielding, self.quartz]
             if r is not None)
        if len(epsgs)!=1:
            raise ValueError("non-matching projection (epsg) of raster datasets")
        if self.catchments.epsg not in epsgs:
            raise ValueError("shapefile projection (epsg) not matching raster datasets")
        
        for r in rr:
            label = r.dtype

            try:
                Z = clip_raster(polygon, r.src, label, name)
            except ValueError:
                raise OutOfBoundsError("catchment polygon out of bounds"
                                 ) from None
            clips[label] = Z
            
        for r in {'elevation', 'shielding', 'quartz'}:
            if r in clips.keys():
                if (sum(~np.isnan(clips[r].values.flatten()))==0): # all empty slice
                    raise ValueError("{} values are nodata".format(r))
            
        for r in {'shielding', 'quartz'}:
            if r in clips.keys():            
                if (np.nanmin(clips[r].values)<0 or np.nanmax(clips[r].values)>1):
                    raise ValueError("{} values are outside of valid range".format(r))
        
        # out of bounds is caught by validation of the clipped catchment
        if any(x is None for x in clips.values()):
            raise OutOfBoundsError("catchment polygon out of bounds"
                             ) from None
        
        clips['epsg'] = epsgs.pop() # validated above
        clips['name'] = name # also stored with each clip['elevation].attrs['name']
        return clips



def main():
    pass
    
        
if __name__ == "__main__":
    main()

    
