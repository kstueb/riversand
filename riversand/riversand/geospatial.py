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

import os, sys

import numpy as np
import pandas as pd
import collections
from numbers import Number

import rasterio
import rasterio.crs
import rasterio.mask

from pyproj.crs import CRS
import pyproj.exceptions

import fiona

import warnings
warnings.filterwarnings("error") # convert warnings to raise exceptions

# =============================================================================
# Raster datasets
# =============================================================================

def get_geotiff(fname:str) -> rasterio.DatasetReader:
    """ Get file handle to geotiff. """

    src = None
    
    if not os.path.isfile(fname):
        raise FileNotFoundError("{}: No such file or directory".format(fname))
    
    try:
        src = rasterio.open(fname, 'r')
        src.close()
    except rasterio.errors.NotGeoreferencedWarning as e:
        raise pyproj.exceptions.CRSError(e.args[0]) # 'Dataset has no geotransform, gcps, or rpcs. The identity matrix will be returned.'
    except:# rasterio.RasterioIOError as error:
        pass
    
    return src


class Raster():
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
        s += ["dtype : {}".format(self.dtype)]
        s += ["fname : {}".format(self.fname)]
        s += ["src   : {}".format(self.src)]
        s += ["epsg  : {}".format(self.epsg)]
        s += ["res   : {}".format(self.res)]
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
        """ Set Catchment.cid (fyi only, Riversand.cid is relevant) """
        
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
        s += ["fname : {}".format(self.fname)]
        s += ["src   : {}".format(self.src)]
        
        s += ["attrs : {}".format(self.attrs)]
        s += ["len   : {}".format(len(self.catchments))]
        s += ["epsg  : {}".format(self.epsg)]
        
        if self.cid:
            s += ["cid   : {}".format(self.cid)]
        
        return "\n".join(s)


# =============================================================================
# Riversand object
# =============================================================================
 

class Riversand():
    """
    A Riversand object contains all the data needed to calculate catchmentwide
    erosion rates. Methods:
    
    .add_rasters()  # add Raster objects (.elevation, .shielding, .quartz)
    .add_catchments()  # add Catchment object (.catchments)
    .add_samples()  # add sample data (.samples)
 
    .validate()   # validate projection and resolution of the geospatial data
        (.res, .crs, .epsg.)
        
    .process_single_catchment()  # calculate erosion rate for a single catchment
    .process_multi_catchment()  # calculate erosion rates for a shapefile with
        multiple catchments (use .set_cid() to set the shapefile attribute to
        be used as catchment name)
    
    """
    

    
    def __init__(self, path:str=None):
        self.path_to_data = None
        self.elevation = None # Raster object
        self.shielding = None # Raster object
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
            elif dtype=='quartz':
                self.quartz = Raster(fname, dtype)
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
        Add / replace self.samples with data from dict or file.
        Prints error messages.
        
        If reading from file exceptions are printed but data is still stored in
        self.samples if possible.
        """
        
        if isinstance(samples, dict):
            try:
                self.add_samples_from_dict(samples, add)
            except Exception as e:
                print("Error adding sample data from dictionary:\n   "+str(e))
            
        elif isinstance(samples, str):
            try:
                self.add_samples_from_file(samples)
            except FileNotFoundError as e:
                    print("{} : {}: No such file"
                          .format(type(e).__name__, e.args[0].split(':')[0]))
            except Exception as e:
                print("Error adding sample data from file:\n   "+str(e))
        else:
            raise TypeError("add_samples() missing mandatory argument "+
                            "(file name or dictionary)")

                
    def add_samples_from_file(self, fname:str):
        """
        Write sample data to self.samples.
        
        Raises ValueErrors "\n".join(errors) found by:
        errR, errS, errC = self.validate(verbose=False)
        
        """
        
        from riversand import params
                    
        if self.path_to_data is not None:
            fname = os.path.join(self.path_to_data, fname)
            
        df = pd.DataFrame()
        
        try:
            df = pd.read_csv(fname)
            
            # I believe this is only needed for csv files
            for c in df.columns:
                try:
                    df[c] = df[c].str.strip(' \n\t')
                except:
                    pass
        except:
            try:
                xls = pd.ExcelFile(fname)
            except ImportError as e:
                raise e
            except IOError as e:
                if e.errno==2:
                    raise FileNotFoundError("{}: No such file or directory"
                                            .format(e.filename))
                else:
                    raise e
            else:
                df = xls.parse(0)
    
        for c in df.columns:
            c_new = c.strip(' \n\t')
            df.rename(columns={c: c_new}, inplace=True)
        
        if df.empty:
            raise ValueError("No sample data in file")
        
        topo_keys = ['lat', 'long', 'elevation']
        mandatory_keys = [k for k in params.all_keys if k not in topo_keys]
        other_keys = [k for k in topo_keys if k in df.keys()] # additional 'lat', 'long', 'elevation'
        other_keys += [k for k in df.keys() if k not in params.all_keys] # additional other
        
        self.samples = pd.DataFrame(df, columns=mandatory_keys+other_keys) # sorted columns
        self.samples.dropna(axis=1, how='all', inplace=True) # remove columns that are all nan
        
        self.valid_catchments = None # set by self.validate()
        
        _, errS, _ = self.validate(verbose=False)
        if errS:
            raise ValueError(errS[0])
        
        return errS # list of one error message (sample data)
        

    def add_samples_from_dict(self, data:dict, add:bool=False):
        """
        Write sample data to self.samples or add data.
        
        """
        
        from riversand import params
        from riversand.utils import validate_nuclide, validate_sample
            
        if isinstance(data, (pd.DataFrame)):
            raise NotImplementedError("add_samples_from_dict() not yet implemented for pandas DataFrame")
            
        if not isinstance(data, (dict, pd.Series)):
            raise TypeError("add_samples_from_dict() 'data' must be dictionary")
        
        if add==False:
            self.samples = None
                
        if add==True and (
                (self.samples is None) or
                (not isinstance(self.samples, pd.DataFrame))):
            add = False # can only add to existing data
            
        topo_keys = ['lat', 'long', 'elevation']
        mandatory_keys = [k for k in params.all_keys if k not in topo_keys]  
        
        if 'shielding' not in data.keys(): #default shielding=1 if not defined 
            data['shielding'] = 1
        
        exc = None
        try:
            nuclide_vals = validate_nuclide(data)
            sample_vals = validate_sample(data, shielding=True)
        except ValueError as exc:
            raise ValueError("Invalid sample data: {}".format(exc)) from exc
        else:
            other_keys = [k for k in topo_keys if k in data.keys()] # additional 'lat', 'long', 'elevation'
            other_keys += [k for k in data.keys() if k not in params.all_keys] # additional other
            other_vals = dict((k, data[k]) for k in other_keys)
        
            
            if add==False: # new dataset or adding to empty Dataframe
                self.samples = pd.DataFrame(nuclide_vals|sample_vals|other_vals,
                                            columns=mandatory_keys+other_keys,
                                            index=[0])       
            else:
                self.samples.reset_index(inplace=True, drop=True)
                # add new data to the end but use column ordering off added (validated) data
                df = pd.DataFrame(nuclide_vals|sample_vals|other_vals,
                                  columns=mandatory_keys+other_keys,
                                  index=[len(self.samples)])
                self.samples = (pd.concat([df, self.samples])
                                .sort_index().reset_index(drop=True))        
                
            self.valid_catchments = None # set by self.validate()
        
    
    def get_valid_catchments(self):
        self.validate(verbose=False)
        if self.valid_catchments is None:
            return []
        else:
            return self.valid_catchments
    
        
    def validate(self, multi:bool=None, verbose=True):
        """
        Validate Riversand object and print error report.
                 
        """
        
        from riversand.utils import validate_topo, validate_nuclide, validate_sample
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
            if self.quartz is not None:
                proj_errs += int(self.quartz.epsg!=self.elevation.epsg)
                proj_errs += int(self.quartz.res!=self.elevation.res)
                #print("checking quartz "+str(proj_errs))
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
            
            if not outR==[]: # outR=None if not validated; outR=["err mssg"] if invalid 
                raise ValidationError("No valid raster data, cannot validate shapefile projection")
        
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
        Calculate catchmentwide erosion rates for a single-catchment shapefile.

        Parameters
        ----------
        bins : numeric, optional
            bin size in meters for elevation binning. The default is 500.
        scaling : str, optional
            scaling method 'St', 'Lm' or 'LSDn'. The default is 'LSDn'.
        shielding : str or numeric, optional
            shielding method, 'sample', 'topo' or numeric (0 to 1).
            The default is 1.
        unit : str, optional
            unit for plotting. The default is 'mm/yr'.
        plot : str, optional
            flag for savinf plots, 'jpg' or 'png'. The default is None.

        Returns
        -------
        results : pd.DataFrame
            Table of results for each sample.

        """
        
        from riversand import params
        if url is None: url = params.url
            
        #from riversand.utils import validate_topo, validate_nuclide, validate_sample
        from riversand.utils import eliminate_quartzfree, get_topostats
        from riversand.calc import get_textline, get_E, guess_erates
        from riversand.calc import poly_E_results, get_RMSE
        import riversand.plot
        
        if scaling not in params.scalingmethods:
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
            
        if unit not in params.units.keys():
            raise ValueError("Invalid unit='{}': see params.units for valid options"
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
            if plot=='png':
                print("Saving plots as .png in '{}'".format(params.out_path))
            if plot in ('jpg', 'auto'):
                print("Saving plots as .jpg in '{}'".format(params.out_path))
            print("")
        
        try:
            clips = self.clip_all_rasters() # function calls self.validate()
        except ValueError as e:
            # catchment is out of bounds; this shold never happen for
            # single-catchment datasets
            print("Catchment is out of bounds")
            return results
            
        else:
            if 'quartz' in clips.keys():
                clips = eliminate_quartzfree(clips) # prints report message
                
            if plot is not None:
                for label in clips.keys():
                    if label in Raster.dtypes:
                        riversand.plot.plot_clipped_raster(clips, c_name='',
                                                 dtype=label, fname=plot)
                        
            topostats, summary = get_topostats(clips, bins=bins, centroid='from_clipped') # accepts iterable as bins
            
            try: # mock textline to make sure that 'shielding' is set correctly
                summary['elevation'] = 100
                _ = get_textline(self.samples.iloc[0],
                                 summary,
                                 shielding=shielding)
            except ValueError as e: # raised by get_textline() for invalid shielding
                print(type(e).__name__, ":", e.args[0])
                return results
                
            for idx, sample_data in self.samples.iterrows(): # sample_data as pd.Series
                
                E = np.nan
                delE = (np.nan, np.nan)
                NRMSE = np.nan
            
                # sample name from input or default value
                try:
                    name = sample_data['name']
                except:
                    name = params.default_values['name']
                
                # nuclide and mineral from input or default value
                try:
                    nuclide = sample_data['nuclide']+' '
                except:
                    nuclide = params.default_values['nuclide']+' '
                try:
                    nuclide += sample_data['mineral']
                except:
                    nuclide += params.default_values['mineral']
                                    
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
                    
                    E, delE, NofE, err = (
                        poly_E_results(sample_data, topostats, 
                        shielding=shielding, erates=erates, scaling=scaling,url=url))
                    
                    # catch out-of-bound errors
                    while 'minE too high' in err:
                        erates = guess_erates(3*erates[0]-2*erates[1], erates[-1]) # cm/yr
                        E, delE, NofE, err = (
                            poly_E_results(sample_data, topostats, 
                            shielding=shielding, erates=erates, scaling=scaling,url=url))
                        
                    while 'maxE too low' in err:
                        erates = guess_erates(erates[0], 3*erates[-1]-2*erates[-2]) # cm/yr
                        E, delE, NofE, err = (
                            poly_E_results(sample_data, topostats, 
                            shielding=shielding, erates=erates, scaling=scaling,url=url))
                    
                    NRMSE = get_RMSE(NofE)/sample_data['N']
                    
                    if plot is not None:
                        if plot in {'jpg', 'png'}:
                            fullname = "{}_{}_{}.{}".format(idx, name, scaling, plot)
                        else:
                            fullname = "{}_{}_{}.jpg".format(idx, name, scaling)
                        riversand.plot.plot_polyfit(E, delE, NofE, sample_data,
                                                    unit=unit, fname=plot,
                                                    fullname=fullname)
                    
                finally:
                    results.loc[idx, 'name'] = name    
                    results.loc[idx, 'scaling'] = scaling
                    results.loc[idx, 'nuclide'] = nuclide
                    
                    results.loc[idx, 'E'] = E
                    results.loc[idx, 'delE-'] = delE[0]
                    results.loc[idx, 'delE+'] = delE[1]
                    results.loc[idx, 'NRMSE'] = NRMSE
                    results.loc[idx, 'Tavg'] = np.round(160/sample_data['density']/E,)
            
                    if err==[]:
                        results.loc[idx, 'error'] = ''
                    else:
                        results.loc[idx, 'error'] = "; ".join(err)
                    ## possible error messages:
                    # 'NRMSE = ...'
                    # 'Root finding did not converge'
                    # 'Cannot compute uncertainty'
                    # 'Server cannot resolve erosion rates, dropping duplicates'
        
        return results
    
        
    def process_multi_catchment(self,
            bins=500, scaling='LSDn', shielding=1, unit='mm/yr',
            plot=None,
            url=None, verbose=True) -> pd.DataFrame:
        """
        Calculate catchmentwide erosion rates for a multi-catchment shapefile.
        
        Parameters
        ----------
        bins : numeric, optional
            bin size in meters for elevation binning. The default is 500.
        scaling : str, optional
            scaling method 'St', 'Lm' or 'LSDn'. The default is 'LSDn'.
        shielding : str or numeric, optional
            shielding method, 'sample', 'topo' or numeric (0 to 1).
            The default is 1.
        unit : str, optional
            unit for plotting. The default is 'mm/yr'.
        plot : str, optional
            flag for savinf plots, 'jpg' or 'png'. The default is None.
    
        Returns
        -------
        results : pd.DataFrame
            Table of results for each sample.
            
        """
        
        from riversand import params
        if url is None: url = params.url
            
        #from riversand.utils import validate_topo, validate_nuclide, validate_sample
        from riversand.utils import eliminate_quartzfree, get_topostats
        from riversand.calc import get_textline, get_E, guess_erates
        from riversand.calc import poly_E_results, get_RMSE
        import riversand.plot
        
        if scaling not in params.scalingmethods:
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
            
        if unit not in params.units.keys():
            raise ValueError("Invalid unit='{}': see params.units for valid options"
                             .format(unit))
            
            
        result_cols = ['name', 'scaling', 'nuclide', 'qtz',
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
            if plot=='png':
                print("Saving plots as .png in '{}'".format(params.out_path))
            if plot in ('jpg', 'auto'):
                print("Saving plots as .jpg in '{}'".format(params.out_path))
            print("")
            print("")

        # some printout formatting
        If = len(str(len(self.samples))) # fill length for id
        Nf = max([len(l) for l in self.samples['name'].values]) # fill length for sample name
        Ef = 26 # fill length for erosion rate result        
        order_str = "{: >"+str(If)+"} {: <"+str(Nf)+"} : "
        #order_str_E2 = "{: >"+str(If)+"} {: <"+str(Nf)+"} : {: >"+str(Ef)+"}"  # erosion rate result
        order_str_F1 = "{: >"+str(If)+"} {: <"+str(Nf)+"} : "  # erosion rate result
        #order_str_F2 = "{: >"+str(If)+"} {: <"+str(Nf)+"} : {}"  # erosion rate result

        
        for idx, sample_data in self.samples.iterrows(): # iterate over all samples
                
            Qpc = 0 # quartz-free lithology removed; default value
            E = np.nan
            delE = (np.nan, np.nan)
            NRMSE = np.nan
            error_code = None

            name = sample_data['name']
            # nuclide and mineral from input or default value
            try:
                nuclide = sample_data['nuclide']+' '
            except:
                nuclide = params.default_values['nuclide']+' '
            try:
                nuclide += sample_data['mineral']
            except:
                nuclide += params.default_values['mineral']
                
            if verbose: # print index and sample name
                print(order_str.format(idx, name), end='')
            
            if sample_data['name'] not in self.valid_catchments: #set by self.validate()
                error_code = 'no catchment polygon'
            else:
                # get integer id of the catchment
                n = [n for n, c in enumerate(self.catchments.catchments) 
                     if c['properties'][self.cid]==name][0]
                try:
                    clips = self.clip_all_rasters(n) # function calls rv.validate()
                except ValueError as e:
                    #print(e) # out of bounds ; is checked by .validate() and should not happen for single-catchment 
                    error_code = 'catchment out of bounds'
                else:
                    if 'quartz' in clips.keys():
                        clips, Qpc = eliminate_quartzfree(clips, verbose=False, Qpc=True)
                
                    if plot is not None:
                        for label in clips.keys():
                            if label in Raster.dtypes:
                                riversand.plot.plot_clipped_raster(
                                    clips, c_name=name,
                                    dtype=label, fname=plot)
                                
                    topostats, summary = get_topostats(clips, bins=bins, centroid='from_clipped') # accepts iterable as bins
                    
                    
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
                        error_code = e.args[0][10:]
                    else:
            
                        # generate suitable erates
                        erates = guess_erates(E_Lo, E_Hi, scaling=scaling) # cm/yr
                        
                        E, delE, NofE, err = (
                            poly_E_results(sample_data, topostats, 
                            shielding=shielding, erates=erates, scaling=scaling,url=url))
                
                        # catch out-of-bound errors
                        while 'minE too high' in err:
                            erates = guess_erates(3*erates[0]-2*erates[1], erates[-1]) # cm/yr
                            E, delE, NofE, err = (
                                poly_E_results(sample_data, topostats, 
                                shielding=shielding, erates=erates, scaling=scaling,url=url))
                            
                        while 'maxE too low' in err:
                            erates = guess_erates(erates[0], 3*erates[-1]-2*erates[-2]) # cm/yr
                            E, delE, NofE, err = (
                                poly_E_results(sample_data, topostats, 
                                shielding=shielding, erates=erates, scaling=scaling,url=url))
                        
                        NRMSE = get_RMSE(NofE)/sample_data['N']

                        # printout formatting
                        Estr = "{:.1f}+/-{:.1f} {}".format(E*params.units[unit],
                                                        delE[1]*params.units[unit],
                                                        unit)
                        
                        if plot is not None:
                            if plot in {'jpg', 'png'}:
                                fullname = "{}_{}_{}.{}".format(idx, name, scaling, plot)
                            else:
                                fullname = "{}_{}_{}.jpg".format(idx, name, scaling)
                            riversand.plot.plot_polyfit(E, delE, NofE, sample_data,
                                                        unit=unit, fname=plot,
                                                        fullname=fullname)
                        if verbose: # print erosion rate result
                            print(Estr)
                            
                       
                        
            results.loc[idx, 'name'] = sample_data['name']   
            results.loc[idx, 'scaling'] = scaling
            results.loc[idx, 'nuclide'] = nuclide
            results.loc[idx, 'qtz'] = np.round(100-Qpc,1)
                        
            results.loc[idx, 'E'] = E
            results.loc[idx, 'delE-'] = delE[0]
            results.loc[idx, 'delE+'] = delE[1]
            results.loc[idx, 'NRMSE'] = NRMSE
            results.loc[idx, 'Tavg'] = np.round(160/sample_data['density']/E,)
                
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
                # 'Server cannot resolve erosion rates, dropping duplicates'
        
        try: # if all went wrong there might not be a variable clips
            if 'quartz' not in clips.keys():
                results.drop(columns=['qtz'], inplace=True)
        except:
            pass
        
        return results
    
    
    def catchment_stats(self,
            bins=100, verbose=True) -> pd.DataFrame:
        """
        Calculate topo statistics for a multi-catchment shapefile.
        
        Parameters
        ----------
        bins : numeric, optional
            bin size in meters for elevation binning. The default is 500.
    
        Returns
        -------
        results : pd.DataFrame
            Table of results for each sample.
            
        """
        
        from riversand.utils import eliminate_quartzfree, get_topostats
        
        if self.elevation is None:
            raise ValueError("Missing elevation raster")
        if self.samples is None:
            raise ValueError("Missing sample data")           
        if self.catchments is None:
            raise ValueError("Missing catchment polygons")
            
        # # multi / single is determined by the number of polygons but cid must be set to identify the polygons
        # if self.cid is None:
        #     raise ValueError("No catchment identifier defined; use .set_cid()")
        # if self.cid not in self.catchments.attrs: # self.set_cid() does not allow to set an invalid cid
        #     raise ValueError("Invalid catchment identifier cid='{}'; use .set_cid()"
        #                                   .format(self.cid))
            
        topo_cols = ['name', 'centr_lat', 'centr_long', 'mean_elev',
                     'relief', 'area', 'mean_sf', 'qtz_pc']
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
                name = c['properties'][self.cid]
            else:
                name = ""
            results.loc[idx, 'name'] = name
            if verbose:
                print(end+name, end="")
                end = ", "
                
            try:
                clips = self.clip_all_rasters(idx) # function calls rv.validate()
            except ValueError as e:
                results.loc[idx, 'error'] = 'catchment out of bounds'
            else:
                if 'quartz' in clips.keys():
                    clips, Qpc = eliminate_quartzfree(clips, verbose=False, Qpc=True)
                    results.loc[idx, 'qtz_pc'] = np.round(100-Qpc,1)
                topostats, summary = get_topostats(clips, bins=bins, centroid='from_clipped') # accepts iterable as bins
                
                results.loc[idx, 'area'] = summary['areakm2']
                results.loc[idx, 'relief'] = np.nanmax(clips['elevation'])-np.nanmin(clips['elevation'])
                if 'shielding' in clips.keys():
                    results.loc[idx, 'mean_sf'] = np.nanmean(clips['shielding'])
                results.loc[idx, 'centr_lat'] = summary['lat']
                results.loc[idx, 'centr_long'] = summary['long']
                results.loc[idx, 'mean_elev'] = np.nanmean(clips['elevation'])
        
        if verbose:
            print(" finished.")
        
        results = results.sort_values(by=['name'], ignore_index=True)
        
        return results
    
    
    def __repr__(self):
        s = []
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
            s += ["---------------"]
            s += ["Validated projection:"]
            if self.epsg:
                s += ["epsg  : {}".format(self.epsg)]
            if self.res:
                s += ["res   : {}".format(self.res)]
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
    


    def clip_all_rasters(self, n:int=None) -> dict:
        """
        Clip all rasters 'elevation', 'shielding', 'quartz' to a single polygon.
        
        n : number of the catchment
            polygon = rv.catchment.catchments[n]['geometry']
            for single-catchment shapefile set n=0 or or n=None (default value)
        
        """
        
        from riversand.utils import clip_raster
        
        if n is None:
            n=0
        
        #if self.epsg is None:
        #    print("WARNING : dataset does not seem validated, run .validate()")
        errR, errS, errC = self.validate(verbose=False)
        if len(errR+errC)>0:
            raise ValueError("Cannnot validate 'rv'; use 'rv.validate()' for details")
        
        try:
            polygon = self.catchments.catchments[n]['geometry'] # raises IndexError
            #try:
            #    c_name = self.catchments.catchments[n]['properties'][self.cid]
            #except:
            #    c_name = '' # for Exception message
        except IndexError as e:
            raise ValueError("clip_all_rasters() : 'catchments' does not have "+
                             "n={} polygons".format(n))
        
        clips = {}
        
        rr = [r for r in [self.elevation, self.shielding, self.quartz]
             if r is not None]
        
        for r in rr:
            label = r.dtype
            try:    
                Z = clip_raster(polygon, r.src, label)
            except ValueError as e:
                raise ValueError("clip_all_rasters() : catchment polygon out of bounds")
            clips[label] = Z
            
        # out of bounds is caught by validation of the clipped catchment
        if any(x is None for x in clips.values()):
            raise ValueError("clip_all_rasters() : catchment polygon out of bounds")
        
        clips['epsg'] = self.epsg # validated at the beginning of the function
        return clips



def main():
    pass
    
        
if __name__ == "__main__":
    main()

    
