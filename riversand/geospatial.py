#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 01:15:49 2022

***** geospatial.py ***********************************************************

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
    
The current version of the online calculator is obtained from calc.py
version = calc.get_version()


@author: Konstanze Stübner, kstueb@gmail.com

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
    except:# rasterio.RasterioIOError as error:
        pass
    
    return src
    
    
def get_epsg(rio_crs:rasterio.crs.CRS) -> int:
    """
    Get epsg code of raster coordinate reference system.
    > get_epsg(src.crs)
    
    """
    
    epsg = None
    try: 
        crs = CRS.from_user_input(rio_crs) 
        #crs = CRS.from_epsg(crs.to_epsg()) # this trick adds the area_of_use to raster_crs
        epsg = crs.to_epsg()
    except pyproj.exceptions.CRSError as error:
        print(error)

    if epsg is None:
        try:
            wkt = rio_crs.wkt
            crs = CRS.from_wkt(wkt)
            #crs = CRS.from_epsg(crs.to_epsg()) # this trick adds the area_of_use to raster_crs
            epsg = crs.to_epsg()
        except AttributeError:
            pass
        except pyproj.exceptions.CRSError as error:
            print(error)
            
    # crs may not be specified in src, eg if the raster is Matlab generated without the Mapping toolbox   
    return epsg


class Raster():
    dtypes = ('elevation', 'shielding', 'quartz') # valid Raster types
    
    def __init__(self, fname:str, dtype:str):
        
        dtype = dtype.lower()
        if dtype not in self.dtypes:
            raise TypeError("Invalid Raster dtype '{}' (must be 'elevation', 'shielding' or 'quartz')".format(dtype))
        
        if not isinstance(fname, str):
            raise TypeError("Invalid Raster fname (must be string)")
            
        src = get_geotiff(fname) # closed rasterio.DatasetReader
        if src:
            self.dtype = dtype
            self.fname = fname
            self.src = src
            self.epsg = get_epsg(src.crs) # epsg code, int
            self.res = src.res # resolution, tuple
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
    dtypes = ('single', 'multi', '') # valid Catchment types
    
    def __init__(self, fname:str, dtype:str=''):
        
        dtype = dtype.lower()
        if dtype not in self.dtypes:
            raise TypeError("Invalid Catchment dtype '{}' (must be 'single' or 'multi')"
                            .format(dtype))

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
                    
                crs = None
                try: # get crs from wkt
                    crs = CRS.from_wkt(src.crs_wkt)
                except: # get crs from 'init'
                    crs = src.crs
                    try:
                        crs = CRS.from_epsg(crs['init'].split(":")[1])
                    except:
                        print("Cannot determine shapefile projection")
                    
            if dtype=='':
                if len(catchments)==1:
                    self.dtype = 'single'
                else:
                    self.dtype = 'multi'
            else:
                self.dtype = dtype
                
            self.fname = fname
            self.src = src
            
            self.attrs = list(attrs)
            self.crs = crs
            self.catchments = catchments
            
            
            self.cid = None

        else:
            raise IOError("Cannot read shapefile {}".format(fname))
            
        
    def set_dtype(self, dtype:str):
        if dtype not in self.dtypes:
            raise TypeError("Invalid Catchment dtype '{}' (must be 'single' or 'multi')".format(dtype))
        self.dtype = dtype
        
        
    def set_cid(self, cid:str='id'):
        """ Set Catchment.cid; this parameter is fyi only. """
        
        if cid is None:
            self.cid = None
        else:
            if self.attrs is None:
                print("'{}' is not a valid catchment identifier;".format(cid))
                print("   shapefile has no fields")
                self.cid = None
            if cid in self.attrs:
                self.cid = cid
            else:
                print("'{}' is not a valid catchment identifier;".format(cid))
                print("   shapefile has fields: {}".format(", ".join(self.attrs)))
                self.cid = None


    def get_names(self, cid:str=None) -> list:
        """
        Get sorted list of all catchment names from attribute field 'cid'.
        Duplicates are included, missing values as 'None'.
        
        """
        
        if cid is None:
            try:
                cid = self.cid
            except:
                raise ValueError("No catchment identifier set; use function .set_cid()")
        
        if cid not in self.attrs:
            raise ValueError("Invalid catchment identifier '{}';\n"
                                 .format(cid) +
                                 "   attribute fields are: {}"
                                 .format(", ".join(self.attrs)))
         
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
        
        c_names = self.get_names(cid)        
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
        s += ["epsg  : {}".format(self.crs.to_epsg())]
        
        if self.cid:
            s += ["cid   : {}".format(self.cid)]
        
        return "\n".join(s)


# =============================================================================
# Riversand object
# =============================================================================
 

class Riversand():
    
    def __init__(self):
        self.elevation = None # Raster object
        self.shielding = None # Raster object
        self.quartz = None # Raster object
        
        self.res = None # project resolution
        self.epsg = None # project projection
        self.cid = None # catchment identifier for multi-catchment datsets
        
        self.catchments = None # Catchment object

        self.samples = None # pd.DataFrame of sample data
        
        self.valid_catchments = None # set by self.validate('catchment')
        # not validated if within raster bounds or if sample data is valid
        # set to None by .add_catchments(), .set_cid(),
        #     .add_samples_from_file(), .add_samples_from_dict()
            
    def add_raster(self, fname:str, path:str='', dtype:str='elevation'):
        """ Add raster dataset to the project. """
        if isinstance(fname, Raster):
            raise NotImplementedError("add_raster() not yet implemented for Raster object")
            
        fullname = os.path.join(path, fname)
        
        dtype = dtype.lower()
        if dtype=='elevation':
            self.elevation = Raster(fullname, dtype)
        elif dtype=='shielding':
            self.shielding = Raster(fullname, dtype)
        elif dtype=='quartz':
            self.quartz = Raster(fullname, dtype)
        else:
            Raster(fullname, dtype) # raises exceptions
        
        self.epsg = None
        self.res = None
        
        
    def add_catchments(self, fname:str, path:str='', dtype:str='single'):
        """ Add catchment shapefile to the project. """
        
        if isinstance(fname, Catchment):
            raise NotImplementedError("add_catchments() not yet implemented for Catchment object")
            
        fullname = os.path.join(path, fname)
        dtype = dtype.lower()
        self.catchments = Catchment(fullname, dtype)
        self.cid = None
        
        self.valid_catchments = None # set by self.validate()
        
        
    def set_cid(self, cid:str='id'):
        """ Set attribute 'cid'. """
        
        if self.catchments is None:
            raise ValueError("No catchments defined")
        if self.catchments.attrs is None:
            print("'{}' is not a valid catchment identifier;".format(cid))
            print("   shapefile has no fields")
            self.cid = None
            self.catchments.set_cid(None)
        if cid in self.catchments.attrs:
            self.cid = cid
            self.catchments.set_cid(cid)
        else:
            print("'{}' is not a valid catchment identifier;".format(cid))
            print("   shapefile has fields: {}".format(", ".join(self.catchments.attrs)))
            self.cid = None
            self.catchments.set_cid(None)
            
        self.valid_catchments = None # set by self.validate()
        
        
    def add_samples(self, **kwargs):
        """
        Add / replace self.samples with data from dict or file.
        Prints error messages.
        
        If reading from file exceptions are printed but data is still stored in
        self.samples if possible.
        """
        
        if 'data' in kwargs.keys():
            try:
                self.add_samples_from_dict(**kwargs)
            except Exception as e:
                print("Error adding sample data from dictionary:\n   "+str(e))
            
        elif 'fname' in kwargs.keys():
            try:
                self.add_samples_from_file(**kwargs)
            except Exception as e:
                print("Error adding sample data from file:\n   "+str(e))
        else:
            raise TypeError("add_samples() missing required keyword argument 'data' or 'fname'")

                
    def add_samples_from_file(self, **kwargs):
        """
        Write sample data to self.samples.
        kwargs: fname, path
        
        Raises ValueErrors "\n".join(errors) found by:
        errors = self.validate('sample', verbose=False)
        
        """
        
        from riversand import params
        
        if 'fname' in kwargs.keys():
            fname = kwargs['fname']
        else:
            raise TypeError("add_samples_from_file() missing required keyword argument 'fname'")
            
        if 'path' in kwargs.keys():
            path = kwargs['path']
        else:
            path = ''
                    
        if not isinstance(fname, str):
            raise TypeError("add_samples_from_file() 'fname' must be string")
        if not isinstance(path, str):
            raise TypeError("add_samples_from_file() 'path' must be string")
            
        df = pd.DataFrame()
        
        try:
            df = pd.read_csv(os.path.join(path, fname))
            
            # I believe this is only needed for csv files
            for c in df.columns:
                try:
                    df[c] = df[c].str.strip(' \n\t')
                except:
                    pass
        except:
            try:
                xls = pd.ExcelFile(os.path.join(path, fname))
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
        
        errors = self.validate('sample', verbose=False)
        if errors:
            raise ValueError("\n".join(errors))
        
        return errors
        

    def add_samples_from_dict(self, **kwargs):
        """
        Write sample data to self.samples or add data.
        kwargs: data, add (default add=False)
        
        """
        
        from riversand import params
        from riversand.utils import validate_nuclide, validate_sample
        
        if 'data' in kwargs.keys():
            data = kwargs['data']
        else:
            raise TypeError("add_samples_from_dict() missing keyword argument 'data'")
            
        if isinstance(data, (pd.DataFrame)):
            raise NotImplementedError("add_samples_from_dict() not yet implemented for pandas DataFrame")
            
        if not isinstance(data, (dict, pd.Series)):
            raise TypeError("add_samples_from_dict() 'data' must be dictionary")
        
        if 'add' in kwargs.keys() and not (kwargs['add']==False):
            add = True
            if (self.samples is None) or not (isinstance(self.samples, pd.DataFrame)):
                add = False # can only add to existing pd.DataFrame
        else:
            add = False
            self.samples = None
                
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
        
            
    def validate(self, *args, **kwargs):
        """
        Validate Riversand object and print error report.
                          
        kwargs :
            'verbose'=False : returns list of errors instead of printing report.
            'dtype'=str : forces catchment dtype of 'single' or 'multi.
                Raises ValueError if 'single' for a multi-polygon shapefile;
                this is the only error this function should ever raise.
        
        """
        
        from riversand.utils import validate_topo, validate_nuclide, validate_sample
        from riversand.utils import feature_in_raster

        if 'verbose' in kwargs.keys():
            verbose = kwargs['verbose']
        else:
            verbose = True
            
        if 'dtype' in kwargs.keys():
            c_dtype = kwargs['dtype'] # 'single' or 'multi'
        else:
            c_dtype = None
            
        args = [a[:6].lower() for a in args] # lowercase and clip to first 6 chars
        args = ['raster' if a=='r' else a for a in args] # accept 'r', 'c', 's' as input
        args = ['catchm' if a=='c' else a for a in args]
        args = ['sample' if a=='s' else a for a in args]

        if len(args)==0:
            args = ['raster', 'catchm', 'sample']
        
        
        
        errR = None
        errS = None
        errC = None
        
        # validate raster datasets
        if 'raster' in (arg.lower() for arg in args):
            errR = []
            
            if not self.elevation:
                errR += ["No elevation raster defined"]
            else:
                proj_errs = 0
                if self.shielding:
                    if self.shielding.epsg != self.elevation.epsg:
                        proj_errs += 1
                    if self.shielding.res != self.elevation.res:
                        proj_errs += 1
                if self.quartz:
                    if self.quartz.epsg != self.elevation.epsg:
                        proj_errs += 1
                    if self.quartz.res != self.elevation.res:
                        proj_errs += 1
                if proj_errs>0:
                    errR += ["Conflicting projections in raster data"]
                    
            if errR:
                self.epsg = None
                self.res = None
            else:
                self.epsg = self.elevation.epsg
                self.res = self.elevation.res

        # validate sample datasets
        if 'sample' in (arg.lower() for arg in args):
            errS = []
            
            if ((self.samples is None) or
                not isinstance(self.samples, pd.DataFrame) or
                self.samples.empty):
                errS += ["No sample data defined"]
            else:
                
                tmp = self.samples.reset_index(drop=True)
                
                # (1) Validate that mandatory columns 'N','delN' exist
                missing_columns = []
                for mc in ['N', 'delN']:
                    if mc not in tmp.keys():
                        missing_columns += [mc]
                if len(missing_columns)>0:
                    errS += ["Missing column(s): {}".format(", ".join(missing_columns))]

                # (2) Check row-wise that validation for online calculator will pass
                if errS==[]:
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
                        errS += ["Invalid / missing data in lines: {}"
                                 .format(", ".join(errors_in_line))]
        
        # validate catchment
        if 'catchm' in (arg.lower() for arg in args):
            errC = []
            
            if ((self.catchments is None) or
                not isinstance(self.catchments, Catchment)):
                errC += ["No catchment data defined"]
            else:
                # validate projection
                if errR is None:
                    errR = self.validate('raster', verbose=False)
                if len(errR)>0:
                    errC += ["Raster data invalid, cannot validate shapefile projection"]
                else:
                    C_epsg = self.catchments.crs.to_epsg()
                    if C_epsg!=self.epsg:
                        errC += ["Shapefile projection (epsg={}) "
                                 .format(C_epsg) +
                                 "does not match raster projection (epsg={})"
                                 .format(self.epsg)]
                
                if c_dtype=='single':
                    self.catchments.set_dtype('single')
                if c_dtype=='multi':
                    self.catchments.set_dtype('multi')
                    
                nC = len(self.catchments.catchments)
                if nC==0:
                    errC += ["No polygons in shapefile"]
                    
                if (nC>1) and c_dtype=='single':
                    raise ValueError("Cannot process as single-catchment"+
                                     "dataset with {} polygons in shapefile"
                                     .format(nC))
                    
                if (nC>1) and (self.catchments.dtype=='single'):
                    self.catchments.set_dtype('multi')
                
                # validate multi-catchment
                if self.catchments.dtype=='multi':
                    if self.cid is None:
                        errC += ["Multi-catchment shapefile; use .set_cid() "+
                                 "to set catchment identifier"]
                    elif self.cid not in self.catchments.attrs:
                        errC += ["Multi-catchment shapefile with invalid cid='{}'\n"
                                 .format(self.cid) +
                                 "   use .set_cid()"]
                    else:                        
                        # valid catchment names, no duplicates, not null, sorted
                        c_names = (self.catchments
                                   .get_valid_names(self.cid))
                        try:
                            s_names = list(self.samples['name'].values)
                        except:
                            s_names = [] # empty list if no column 'name'
                            
                        self.valid_catchments = [c for c in c_names if c in s_names]
  
                # validate single-catchment
                if self.catchments.dtype=='single':
                    out_of_bounds = 0
                    rr = [r.src for r in [self.elevation, self.shielding, self.quartz]
                         if r is not None]
                    for r in rr:
                        try:
                            fir = feature_in_raster(
                                self.catchments.catchments[0]['geometry'],
                                r)
                            if fir==False:
                                out_of_bounds += 1
                        except:
                            pass
                    if out_of_bounds>0:
                        errC += ["Catchment polygon out of bounds of {} raster dataset(s)"
                                 .format(out_of_bounds)]
                        
            
            if self.samples is None: # avoid empty list if self.samples is undefined
                self.valid_catchments = None
                
                
                    
        # print output or return list of errors
        if verbose:            
            if errR is not None:
                if errR:
                    print("Errors in raster data:")
                    for e in errR:
                        print("   {}".format(e))
                else:
                    print("Raster data validated.")

            if errS is not None:
                if errS:
                    print("Errors in sample data:")
                    for e in errS:
                        print("   {}".format(e))
                else:
                    print("Sample data validated.")
                    
            if errC is not None:
                if errC:
                    print("Errors in catchment data:")
                    for e in errC:
                        print("   {}".format(e))
                else:
                    print("Catchment data validated.\n")
                
                if self.valid_catchments is not None:
                    print("Catchments / samples cross-validated:")
                    if len(self.valid_catchments)==0:
                        print("   No matches found")
                    else:
                        print("   Found {} matche(s)".format(len(self.valid_catchments)))
            return None
        else:
            errs = []
            if errR: errs += errR
            if errS: errs += errS
            if errC: errs += errC
            return errs
        
    
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
            raise ValueError("Invalid scaling '{}' (must be 'St', 'Lm' or 'LSDn')")
        
        if isinstance(shielding, str):
            if shielding not in {'topo', 'sample'}:
                raise ValueError("Invalid shielding (must be 'topo', 'sample' or numeric)")
        if isinstance(shielding, Number):
            if not (0<=shielding<=1):
                raise ValueError("Invalid shielding (must be 0..1)")
            
        if unit not in params.units.keys():
            raise ValueError("Invalid unit '{}' (see params.units for valid options)"
                             .format(unit))
            
        result_cols = ['name', 'scaling', 'nuclide',
                       'E', 'delE-', 'delE+', 'NRMSE', 'error']
        results = pd.DataFrame(columns=result_cols)
        
        errors = self.validate(c_type='single', verbose=False)
        if len(errors)>0:
            print("Cannot validate 'rv' for single-catchment processing; "+
                             "use 'rv.validate()' for details")
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
                                                 label=label, fname=plot)
                        
            topostats, summary = get_topostats(clips, bins=bins, centroid='from_clipped') # accepts iterable as bins
        
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
                except RuntimeError as e:
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
            raise ValueError("Invalid scaling '{}' (must be 'St', 'Lm' or 'LSDn')")
            
        if isinstance(shielding, str):
            if shielding not in {'topo', 'sample'}:
                raise ValueError("Invalid shielding (must be 'topo', 'sample' or numeric)")
        if isinstance(shielding, Number):
            if not (0<=shielding<=1):
                raise ValueError("Invalid shielding (must be 0..1)")
            
        if unit not in params.units.keys():
            raise ValueError("Invalid unit '{}' (see params.units for valid options)"
                             .format(unit))
            
        result_cols = ['name', 'scaling', 'nuclide', 'qtz',
                       'E', 'delE-', 'delE+', 'NRMSE', 'error']
        results = pd.DataFrame(columns=result_cols)
        results['name'] = self.samples['name']

        errors = self.validate(c_type='multi', verbose=False)
        if len(errors)>0:
            print("Cannot validate 'rv' for multi-catchment processing; "+
                  "use 'rv.validate()' for details")
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
                                    label=label, fname=plot)
                                
                    topostats, summary = get_topostats(clips, bins=bins, centroid='from_clipped') # accepts iterable as bins
                        
                    try:
                        # estimate minimum and maximum erosion rates
                        summary['elevation'] = summary['elevLo']
                        textline = get_textline(sample_data, summary, shielding=shielding)
                        E_Lo = get_E(textline) # erosion rates in cm/yr
                        
                        summary['elevation'] = summary['elevHi']
                        textline = get_textline(sample_data, summary, shielding=shielding)
                        E_Hi = get_E(textline) # erosion rates in cm/yr
                    except RuntimeError as e:
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
            results.loc[idx, 'qtz'] = 100-Qpc
                        
            results.loc[idx, 'E'] = E
            results.loc[idx, 'delE-'] = delE[0]
            results.loc[idx, 'delE+'] = delE[1]
            results.loc[idx, 'NRMSE'] = NRMSE
                
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
        
        if 'quartz' not in clips.keys():
            results.drop(columns=['qtz'], inplace=True)
        
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
            
        if self.epsg or self.res:
            s += ["Validated projection:"]
            if self.epsg:
                s += ["epsg  : {}".format(self.epsg)]
            if self.res:
                s += ["res   : {}".format(self.res)]
            s += [""]
            
        if not self.samples is None:
            s += ["---------------"]
            s += ["Sample data:"]
            s += ["{:d} sample(s)\n".format(len(self.samples))]
            
        if not self.catchments is None:
            s += ["---------------"]
            s += ["Catchment polygons:\n"]
            s += [str(self.catchments)+"\n"]
            
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
        
        val_errs = self.validate('raster', 'catchment', verbose=False)
        if len(val_errs)>0:
            raise ValueError("Cannnot validate 'rv'; use 'rv.validate()' for details")
        
        try:
            polygon = self.catchments.catchments[n]['geometry'] # raises IndexError
            try:
                c_name = self.catchments.catchments[n]['properties'][self.cid]
            except:
                c_name = '' # for Exception message
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
    """ Testing geospatial.py """
    
    import sys, os
    
    rootdir = "/home/user/Dokumente/11_catchment_erates_calculator/2022 CatchCalc_v0/riversand-1.0"
    sys.path.append(rootdir)
    path = os.path.join(rootdir, "tests/test_data")
    
    import riversand
    
        
if __name__ == "__main__":
    main()

    
