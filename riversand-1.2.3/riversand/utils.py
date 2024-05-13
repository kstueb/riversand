#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

*******************************************************************************
utils.py  :  utilities for the riversand package 

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

Validation functions for online calculator:
    validate_topo()
    validate_nuclide()
    validate_sample()
    - Return validated and complete dictionaries
    - Fill in default values from params.Params if non-mandatory keys are missing
    - Raise ValueError if values are invalid / mandatory values are missing 


restandardize()

    
Geospatial processing functions:
    feature_in_raster()
    clip_raster()
    eliminate_quartzfree()
    get_xarray_centroid()
    get_polygon_centroid()
    projected_xy_to_longlat()
    get_bins()
    get_topostats()
    

import_data()

"""

import numpy as np
import pandas as pd
pd.options.mode.copy_on_write = True

import re
from numbers import Number
from copy import deepcopy

import rasterio
import rasterio.crs
import rasterio.mask
from rasterio import MemoryFile

import xarray as xr

from riversand.params import Params

# =============================================================================
# Validation for online calculator
# =============================================================================

def valid_namestring(name:str) -> bool:
    name_regex = "^[a-zA-Z0-9_-]*$"
    if len(name)>32:
        return False
    if not re.match(name_regex, name):
        return False
    return True


def validate_topo(item, shielding:bool=False) -> dict:
    """
    Generate dict of topo data.
    
    Parameters
    ----------
    item : dict or pd.Series
        Single sample.
        Mandatory keys 'lat', 'long', 'elevation'. Optional key 'shielding'.
    
    Raises
    ------
    ValueError
        Mandatory item missing or ill-defined.
    
    """
           
    out = dict()
    
    # mandatory keys:
    if ('lat' not in item.keys() or
        item['lat'] is None):
        raise ValueError("Latitude is not defined")
        
    if ('long' not in item.keys() or
        item['long'] is None):
        raise ValueError("Longitude is not defined")
        
    if ('elevation' not in item.keys() or
        item['elevation'] is None):
        raise ValueError("Elevation is not defined")
        
    try:
        out['lat'] = float(item['lat'])
    except: # ValueError (e.g. str) / TypeError (e.g. list)
        out['lat'] = np.nan
    if np.isnan(out['lat']):
        raise ValueError("Latitude is not a number")
    if not (-90 <= out['lat'] <= 90):
        raise ValueError("Latitude is out of bounds (-90..+90)")
    
    try:
        out['long'] = float(item['long'])
    except: # ValueError (e.g. str) / TypeError (e.g. list)
        out['long'] = np.nan
    if np.isnan(out['long']):
        raise ValueError("Longitude is not a number")
    if not (-180 <= out['long'] <= 180):
        raise ValueError("Longitude is out of bounds (-180..+180)")
    
    try:
        out['elevation'] = float(item['elevation'])
    except: # ValueError (e.g. str) / TypeError (e.g. list)
        out['elevation'] = np.nan
    if np.isnan(out['elevation']):
        raise ValueError("Elevation is not a number")
    if out['elevation'] < (-500):
        raise ValueError("Elevation is less than -500")
        
    # optional keys:
    if shielding:
        if ('shielding' not in item.keys() or
            item['shielding'] is None):
            raise ValueError("Shielding is not defined")
        try:
            out['shielding'] = float(item['shielding'])
        except: # ValueError (e.g. str) / TypeError (e.g. list)
            out['shielding'] = np.nan
        if np.isnan(out['shielding']):
            raise ValueError("Shielding is not a number")
        if not (0 <= out['shielding'] <= 1):
            raise ValueError("Shielding is out of bounds (0..1)")
            
    return out


def validate_nuclide(item) -> dict:
    """
    Generate dict of nuclide data.

    Parameters
    ----------
    item : dict or pd.Series
        Single sample.
        Mandatory keys 'N', 'delN'. Default values for 'name', 'nuclide',
        'mineral', 'standardization'.
    
    Raises
    ------
    ValueError
        Mandatory item missing or ill-defined.
    
    """
        
    out = dict()
    
    # mandatory keys:
    if ('N' not in item.keys() or
        item['N'] is None):
        raise ValueError("Concentration N is not defined")
        
    if ('delN' not in item.keys() or
        item['delN'] is None):
        raise ValueError("Uncertainty delN is not defined")
        
        
    try:
        out['N'] = round(float(item['N'])) # accepts string "1.2E+6"
    except: # ValueError (e.g. str) / TypeError (e.g. list)
        raise ValueError("Concentration N is not a number") from None
    if out['N'] <= 0:
        raise ValueError("Concentration N is 0 or less")
        
    try:
        out['delN'] = round(float(item['delN']))
    except: # ValueError (e.g. str) / TypeError (e.g. list)
        raise ValueError("Uncertainty delN is not a number") from None
    if out['delN'] <= 0:
        raise ValueError("Uncertainty delN is 0 or less")
        
    # optional keys with default values:
    if 'name' in item.keys():
        #name_regex = "^[a-zA-Z0-9_-]*$"
        try:
            out['name'] = str(item['name'])
        except:
            out['name'] = '??' # illegal character
        
        if not valid_namestring(out['name']):
            raise ValueError("Illegal characters in sample name or >32 characters")
    else:
        out['name'] = Params.default_values['name'] # default value
    
    if 'nuclide' in item.keys():
        try:
            out['nuclide'] = str(item['nuclide'])
        except:
            out['nuclide'] = '??'
        if out['nuclide'] not in {'Be-10', 'Al-26'}:
            raise ValueError("Illegal nuclide")
    else:
        out['nuclide'] =Params.default_values['nuclide'] # default value
        if 'standardization' in item.keys():
            raise ValueError("Must specify nuclide if standardization is specified")
        
    if 'mineral' in item.keys():
        try:
            out['mineral'] = str(item['mineral'])
        except:
            out['mineral'] = '??'
        if out['mineral'] != 'quartz':
            raise ValueError("Illegal mineral")
    else:
        out['mineral'] = Params.default_values['mineral'] # default value

    if out['nuclide'] == 'Be-10':
        if 'standardization' in item.keys():
            try:
                out['standardization'] = str(item['standardization'])
            except:
                out['standardization'] = '??'
        else:
            out['standardization'] = Params._default_standards[out['nuclide']]
        if out['standardization'] not in Params.Be_stds.keys():
            raise ValueError("Illegal Be standardization")

    if out['nuclide'] == 'Al-26':
        if 'standardization' in item.keys():
            try:
                out['standardization'] = str(item['standardization'])
            except:
                out['standardization'] = '??'
        else:
            out['standardization'] = Params._default_standards[out['nuclide']]
        if out['standardization'] not in Params.Al_stds.keys():
            raise ValueError("Illegal Al standardization")    

    return out


def validate_sample(item, shielding:bool=False) -> dict:
    """
    Generate dict of sample data.
    
    Parameters
    ----------
    item : dict or pd.Series
        Single sample.
        Default values for 'press_flag', 'density', 'thickness', 'erate',
        'year'. Optional key 'shielding' .
    
    Raises
    ------
    ValueError
        Mandatory item missing or ill-defined.
    
    """
        
    out = dict()

    # optional keys with default values:
    if 'press_flag' in item.keys():
        try:
            out['press_flag'] = str(item['press_flag'])
        except:
            out['press_flag'] = '??'
        if out['press_flag'] not in {'std','ant'}:
            raise ValueError("Illegal pressure flag")
    else:
        out['press_flag'] = Params.default_values['press_flag'] # default value
            
    if 'density' in item.keys():
        try:
            out['density'] = float(item['density'])
        except:
            out['density'] = np.nan
        if np.isnan(out['density']):
            raise ValueError("Density is not a number")
        if out['density'] <= 0:
            raise ValueError("Density is 0 or less")
    else:
        out['density'] = Params.default_values['density'] # default value
        
    if 'thickness' in item.keys():
        try:
            out['thickness'] = float(item['thickness'])
        except:
            out['thickness'] = np.nan
        if np.isnan(out['thickness']):
            raise ValueError("Thickness is not a number")
        if out['thickness'] < 0:
            raise ValueError("Thickness is less than 0")
    else:
        out['thickness'] = Params.default_values['thickness'] # default value
        
    if 'erate' in item.keys():
        try:
            out['erate'] = float(item['erate'])
        except:
            out['erate'] = np.nan
        if np.isnan(out['erate']):
            raise ValueError("Erosion rate is not a number")
        if out['erate'] < 0:
            raise ValueError("Erosion rate is less than 0")
    else:
        out['erate'] = Params.default_values['erate'] # default value

    if 'year' in item.keys():
        try:
            out['year'] = round(float(item['year']))
        except:
            raise ValueError("Collection year is not a number") from None
        if not (1700 <= out['year'] <= 2100):
            raise ValueError("Collection year is out of bounds (1700..2100)")
    else:
        out['year'] = Params.default_values['year'] # default value
        
    # optional keys:
    if shielding:
        if ('shielding' not in item.keys() or
            item['shielding'] is None):
            raise ValueError("Shielding is not defined")
        try:
            out['shielding'] = float(item['shielding'])
        except:
            out['shielding'] = np.nan
        if np.isnan(out['shielding']):
            raise ValueError("Shielding is not a number")
        if not (0 <= out['shielding'] <= 1):
            raise ValueError("Shielding is out of bounds (0..1)")

    return out


def restandardize(df) -> pd.DataFrame:
    """
    Converts any pd.DataFrame with columns N, delN, nuclide, standardization
    to default standardization.
     
    Raises
    ------
    KeyError
        Mandatory column missing
        
    """

    if not all([k in df.columns for k in {'N', 'delN', 'nuclide', 'standardization'}]):
        raise KeyError("restandardize() needs columns 'N', 'delN', 'nuclide', 'standardization'. "+
                       "Try using import_data() to recast data to required format.")
    
    df['_N'], df['_delN'], df['_standardization'] = df['N'], df['delN'], df['standardization']
    
    df['N'], df['delN'], df['standardization'] = np.nan, np.nan, np.nan
    
    nucls = ['Be-10', 'Al-26']
    stdss = [pd.Series(Params.Be_stds), pd.Series(Params.Al_stds)]
    for nu, st in zip(nucls, stdss): 
        sel = df['nuclide']==nu
        
        df.loc[sel, 'N'] = (
            df.loc[sel,'_standardization'].map(st)
            * df.loc[sel,'_N']
            )
        df.loc[sel, 'delN'] = (
            df.loc[sel,'_standardization'].map(st)
            * df.loc[sel,'_delN']
            )
        df.loc[sel*df['N'].notnull(), 'standardization'] = (
            Params._default_standards[nu]
            )
    df.drop(columns=['_N', '_delN', '_standardization'], inplace=True) 
    return df



# =============================================================================
# Geospatial processing
# =============================================================================
    
def feature_in_raster(polygon:dict, src:rasterio.DatasetReader) -> bool:
    """
    Check whether feature is within the bounds of a raster.
    
    Parameters
    ----------
    polygon : catchment.catchments[0]['geometry']
    
    """
    
    if 'coordinates' not in polygon.keys():
        raise TypeError("feature_in_raster() argument 'polygon' not a valid polygon dictionary")
    if not isinstance(src, rasterio.io.DatasetReader):
        raise TypeError("feature_in_raster() argument 'src' not a valid rasterio.DatasetReader")

    (left, bott, right, top) = rasterio.features.bounds(polygon)
    bounds = src.bounds
    return ((bounds[0] < left) and (bounds[1] < bott) and
            (bounds[2] > right) and (bounds[3] > top))


def clip_raster(polygon:dict, src:rasterio.DatasetReader,
                label:str) -> xr.DataArray:
    """
    Clip raster with catchment polygon and return as xr.DataArray
    with attributes attrs['transform'] and attrs['label'].
    
    Parameters
    ----------
    polygon : dict
        > polygon = catchment.catchments[0]['geometry']
    label : str
        'elevation', 'shielding', 'quartz'
    
    Raises
    ------
    TypeError
    ValueError
        "clip_raster() : catchment polygon out of bounds"
        
    """
    
    if 'coordinates' not in polygon.keys():
        raise TypeError("clip_raster() argument 'polygon' not a valid polygon dictionary")
    if not isinstance(src, rasterio.io.DatasetReader):
        raise TypeError("clip_raster() argument 'src' not a valid rasterio.DatasetReader")

    if src.closed:
        src = rasterio.open(src.name, 'r')

    if not feature_in_raster(polygon, src):
        raise ValueError("clip_raster() : catchment polygon out of bounds")

    # read the raster data with rasterio and clip the catchment
    try:
        Z, Z_transform = rasterio.mask.mask(src, [polygon], nodata=np.nan,
                                            crop=True)
        
    except TypeError as e: # recast raster to float to facilitate nodata=nan
        data = src.read()
        profile = src.profile
        src.close()

        profile.update(dtype='float32', driver='GTiff')
        with MemoryFile() as memfile:
            with memfile.open(**profile) as dst: # open as DatasetWriter
                dst.write(data)
                del data

            src = memfile.open()  # reopen as DatasetReader
            
        Z, Z_transform = rasterio.mask.mask(src, [polygon], nodata=np.nan,
                                            crop=True)
        
    Z_meta = src.meta
    #Z_bounds = src.bounds
    Z_meta.update({'driver': 'GTiff',
           'height': Z.shape[1],
           'width': Z.shape[2],
           'transform': Z_transform})
    
    # Convert raster to xr.DataArray
    xx = (np.arange(Z_meta['width'])*Z_transform[0] +
          Z_transform[2] + 0.5*Z_transform[0])
    yy = (np.arange(Z_meta['height'])*Z_transform[4] +
          Z_transform[5] + 0.5*Z_transform[4])
    Z = xr.DataArray(Z.squeeze(), dims=('y','x'), coords={
                              'x' : xx,
                              'y' : yy },
                    )
    src.close()
    
    Z.attrs['transform'] = Z_transform
    Z.attrs['label'] = label
    
    return Z


def eliminate_quartzfree(clips:dict, verbose=True, Qpc=False) -> dict:
    """
    Removes the quartz-free areas of a catchment determined from the 'quartz'
    raster from all xr.DataArrays in clips.
    
    This function does not modify the values of clips.
    
    clips : dict
        Keys 'quartz', 'elevation', 'shielding'.
        The key 'epsg' set by clip_all_rasters() is ignored.
    
    """
    
    if 'quartz' not in clips.keys():
        raise ValueError("eliminate_quartzfree() argument 'clips' missing "+
                         "required key 'quartz'")
    
    XX = deepcopy(clips)
    
    quartzfree_pc = 0
    
    for dtype in ['elevation', 'shielding']:        
        if dtype in XX.keys():
            A0 = int(XX[dtype].count())
            XX[dtype] = XX[dtype].where(XX['quartz']==1)
            A1 = int(XX[dtype].count())
            quartzfree_pc = 100*(A0-A1)/A0

    if verbose:
        print("Removed {:.1f} % of the catchment as quartz-free"
              .format(quartzfree_pc))
    
    if Qpc:
        return XX, quartzfree_pc # how much quartz-free (percent) has been removed
    else:
        return XX


def get_xarray_centroid(X:xr.DataArray) -> tuple:
    """
    Parameters
    ----------
    X : xr.DataArray
        Projected raster dataset clipped to the catchment polygon.
        > X = clips['elevation']
    
    Returns
    -------
    (centr_x, centr_y) : tuple of floats
        Projected longitude, latitude.
    
    """

    A = xr.DataArray((X.values>0), coords=X.coords)
    centr_x = float(np.sum(A.sum('y')/np.sum(A.sum('y')) * A.x))
    centr_y = float(np.sum(A.sum('x')/np.sum(A.sum('x')) * A.y))

    return (centr_x, centr_y)

    

def get_polygon_centroid(polygon:dict) -> tuple:
    """
    Parameters
    ----------
    polygon : dict
        > polygon = catchment.catchments[0]['geometry']
    
    Returns
    -------
    (centr_x, centr_y) : tuple of floats
        Projected longitude, latitude.
    
    """
    
    if 'coordinates' not in polygon.keys():
        raise TypeError("get_polygon_centroid() argument 'polygon' not a valid polygon dictionary")
        
    from shapely.geometry import Polygon
        
    centr_x = Polygon(polygon['coordinates'][0]).centroid.x
    centr_y = Polygon(polygon['coordinates'][0]).centroid.y

    return (centr_x, centr_y)

    
def projected_xy_to_longlat(xy:tuple, epsg:int) -> tuple:
    """  
    Parameters
    ----------
    xy : tuple
        Projected coordinates
    epsg : int
        epsg code
    
    Returns
    -------
    (long, lat) : tuple of floats
        Unprojected longitude, latitude.
    
    """
    
    from pyproj import Transformer
    # note that 4326 is lat/lon, whereas utm is x/y
    # use option always_xy=True to avoid trouble
    transformer = Transformer.from_crs(epsg, 4326, always_xy=True)
    #transformer = Transformer.from_crs("EPSG:4326", "EPSG:26917")

    (long, lat) =  transformer.transform(xy[0], xy[1])
    return (long, lat)


def get_bins(Z:xr.DataArray, binsize:int=100) -> np.ndarray:
    """
    Get pretty elevation bins.
    
    Parameters
    ----------
    Z : ndarray or xr.DataArray
        Values to be binned.
    
    """

    return np.arange(start=np.floor(Z.min()/binsize)*binsize,
                      step=binsize,
                      stop=np.ceil(Z.max()/binsize)*binsize + binsize
                     )


def get_topostats(clips, bins, centroid='from_clipped',
                  polygon=None,
                  epsg=None,
                  validate_transforms=True,
                  ) -> (pd.DataFrame, dict):
    """
    Compute elevation statistics. 
    
    Parameters
    ----------
    clips : dict of xr.DataArrays
        Keys are raster dtypes ('elevation', 'shielding, ...).
        Values have attributes .transform (Affine) and .label (dtype).
        Additional key 'epsg' has the projection.
        > clips = geospatial.clip_all_rasters(rv)
            
    centroid : str or tuple
        "from_clipped" or (long, lat) in projected coordinates.
    
    Returns
    -------
    df : pd.DataFrame
        Topo statistics.
    summary : dict
        Percentile elevations (35, 50, 65%) and centroid coordinates.
        Keys are 'elevLo', 'elev50', 'elevHi', 'lat', 'long', 'areakm2'.
        
    """

    from riversand.geospatial import Raster
    
    if isinstance(epsg, int):
        pass
    elif 'epsg' in clips.keys(): # clips['epsg']=None for non-validated 
        epsg = clips['epsg']
    if epsg is None:
        raise ValueError("epsg cannot be determined from 'clips'; "+
                         "specify as keyword argument")
    
    if isinstance(centroid, str): # allow mistakes
        centroid = centroid.lower()
        if centroid in 'from_clipped':
            centroid = 'from_clipped'            
    elif (isinstance(centroid, tuple) and len(centroid)==2):
        pass
    else:
        raise TypeError("centroid must be string 'from_clipped' "+
                        "or tuple (long, lat)")
    
    # Z is xr.DataArray of elevation raster
    Z = clips['elevation']
    # XX are all valid raster datasets; specifically excluding clips['epsg']
    XX = dict([(key,val) for key,val in clips.items() if key in Raster.dtypes])
    
    # Perform some validation of the xr.DataArrays.
    # should be redundant if clips = clip_all_rasters()
    if len(clips)>1:
        shape = Z.shape #clips['elevation']
        
        for X in XX.values(): 
            if not X.shape==shape:
                raise ValueError("get_topostats() argument 'clips' with "+
                                 "non-matching shape/resolution")

        if validate_transforms:
            transform = clips['elevation'].transform
            for X in XX.values():
                if not X.transform==transform:
                    raise ValueError("get_topostats() argument 'clips' with "+
                                     "non-matching transforms")

    
    # Compute bins and centroid coordinates
    if isinstance(bins, Number):
        bins = get_bins(Z, bins)
    
    if centroid=='from_clipped':
        (centr_x, centr_y) = projected_xy_to_longlat(
                                              get_xarray_centroid(Z), epsg)
    elif centroid=='from_polygon':
        (centr_x, centr_y) = projected_xy_to_longlat(
                                       get_polygon_centroid(polygon), epsg)
    elif isinstance(centroid, tuple):
        centr_x = float(centr_x)
        centr_y = float(centr_y)

    area_per_pixel = np.abs(clips['elevation'].transform[0] * 
                            clips['elevation'].transform[4])/1000000 # in km2
    
    # Z denotes elevation raster
    # X denotes any other raster
    df = pd.DataFrame()
    i = 0
    for this_bin in zip(bins[:-1], bins[1:]):
        X_bins = {} #dict of rasters for the current bin
        Z_bin = xr.apply_ufunc(np.digitize, Z, this_bin)
        for key, val in XX.items():
            X_bins[key] = val.where(Z_bin==1) # select the appropriate bin
    
        if np.count_nonzero(~np.isnan(Z_bin))>0:
            temp = pd.DataFrame({'bin': this_bin[0], 'bin1': this_bin[1]}|
                                {key: np.nanmean(val) for key, val in X_bins.items()}|
                                {'area': np.count_nonzero(~np.isnan(X_bins['elevation']))*area_per_pixel,
                                 'lat': centr_y, 'long': centr_x}, index = [i])
            df = pd.concat([df, temp], ignore_index=True, sort=False)
            i += 1
            
    df['wt'] = df['area']/sum(df['area'])
    
    #df.insert(0, 'bins', list(zip(df.bin0, df.bin1)))
    #df.drop(columns=['bin0', 'bin1'], inplace=True)
    df.drop(columns=['bin1'], inplace=True)
    
    # no harm having a quartz columns
    if 'quartz' in df.columns:
        df.drop(columns='quartz', inplace=True)
    
    summary = {
        'elevLo' : np.nanpercentile(Z, 35),
        'elev50' : np.nanpercentile(Z, 50.),
        'elevHi' : np.nanpercentile(Z, 65),
        'lat' : centr_y,
        'long' : centr_x,
        'areakm2' : sum(df['area'])
        }
    for k, v in clips.items():
        if k not in ('elevation', 'epsg'):
            summary[k] = np.nanmean(v)
    summary['epsg'] = epsg # avoids that epsg is storead as float
    
    return df, summary


def import_data(df:pd.DataFrame) -> pd.DataFrame:
    """
    Recast pandas DataFrame df for exchange with online calculator. Prints error
    messages; prints warnings if default values are filled in for missing columns.
    Returns an empty DataFrame if problems cannot be resolved.
    
    Parameters
    ----------
    df : pd.DataFrame (see below)
        > df = pd.read_excel('./folder/file.ods') # read from excel spreadsheet
        > df = import_data(df)
                
    Returns
    -------
    out : pd.DataFrame
        columns: {'name', 'lat', 'long', 'elevation', 'press_flag',
                  'thickness', 'density', 'erate', 'shielding', 'year',
                  'nuclide', 'mineral', 'N', 'delN', 'standardization'}

    
    Input DataFrame df must have the following data:
    Sample site information: column names are case-insensitive; several aliases are accepted:
    - 'name'
    - 'lat', 'latitude' :             -90 to +90 decimal degrees
    - 'long', 'longitude', 'lon' :    -180 to +180 decimal degrees
    - 'elevation', 'elev', 'pressure', 'press' :  m or hPa, depending on 'press_flag'
    - 'press_flag':                   'std', 'ant' or 'pre' (standard or antarctica elevation model or pressure in hPa)
    - 'thickness', 'thick' :          sample thickness in cm
    - 'density', 'rho' :              density in g/cm3
    - 'shielding' :                   shielding factor between 0.0 and 1.0
    - 'erate', 'erosion_rate' :       erosion rate in cm/yr
    - 'year' :                        year of sampling betweem 1700 and 2100
    
    Nuclide information: either separate columns for Be and Al data or
    all data in N and nuclide specified in separate column:
    - 'N' :                           nuclide concentration in atoms/g
    - 'delN' (or 'dN', 'N_uncert') :  uncertainty on nuclide concentration in atoms/g
    - 'nuclide' :                     either 'Be-10' or 'Al-26'
    - 'mineral' :                     must be 'quartz'
    - 'standardization', 'standard', 'std' :  see online information
    
    - 'Be' (or 'Be-10', '10Be')
    - 'Be_uncert' (or 'Be_uncertainty' etc)
    - 'Be_standardization' (or 'Be_standard', 'Be_std' etc)
    and/or equivalent for Al
    
    Note that for nuclide information only 'uncertainty', 'standardization' etc
    are case-insensitive. Be, Al and N are case-sensitive.
    
    There are default values if columns 'nuclide', 'mineral' or 'standardization' are missing.
    There are also default values for some of the sample site information; see print_defaults().
    
    """
    # It is recommended to use lowercase letters for all column names except for those referring to nuclide
    # data ('N', 'delN' or 'Be', 'Be_uncert', 'Be_standardization' etc.).
    # All other column names are converted to lowercase, therefore two column names 'name' and 'Name' are
    # considered identical and will cause an error because it is unclear which column to use for sample names.
    
    # missing values raise a ValueError "Non-numeric or missing value in ..."
    # nuclide concentrations <= 0 are silently excluded
    # uncertainties <= 0 raise ValueError "Be uncertainty must be >0."

    dfS, dfN = import_data2(df)
    # raises TypeError if input is not a pd.DataFrame
    # prints messages and returns empty dfS, dfN if input is faulty
    
    out = pd.merge(dfS, dfN, on='name', how='right') # duplicates sample info for duplicate analyses
    out = out[ ['name'] + [ col for col in out.columns if col != 'name' ] ] # move 'name' to beginning
    return out


def import_data2(df:pd.DataFrame) -> (pd.DataFrame, pd.DataFrame):
    """
    Recast pandas DataFrame df for exchange with online calculator.
    Returns separate DataFrames for sample and nuclide information.
    See also import_data().
    
    Parameters
    ----------
    df : pd.DataFrame (see below)
        > df = pd.read_excel('./folder/file.ods') # read from excel spreadsheet
        > dfS, dfN = import_data2(df)
                
    Returns
    -------
    dfS : pd.DataFrame
        columns: {'name', 'lat', 'long', 'elevation', 'press_flag',
                  'thickness', 'density', 'erate', 'shielding', 'year'}
        dfS is empty if input data is faulty
    dfN : pd.DataFrame
        columns: {'name', 'nuclide', 'mineral', 'N', 'delN', 'standardization'}
        dfN is empty if input data is faulty
    
    Raises
    ------
    TypeError
        arg is not pd.DataFrames.
        
    """
    
    
    # create pd.DataFrame dfS with sample site information (each sample exactly once)
    # create pd.DataFrame dfN with nuclide information (each sample of dfS at least once)

    # sample site info:
    dfS = pd.DataFrame(columns = ['name', 'lat', 'long', 'elevation',
                                  'press_flag', 'thickness', 'density',
                                  'shielding', 'erate', 'year'])
    dfN = pd.DataFrame(columns = ['name', 'nuclide', 'mineral',
                                  'N', 'delN', 'standardization'])
    
    if isinstance(df, pd.DataFrame):
        temp = df
    else:
        raise TypeError("import_data() argument 'df' must be pandas DataFrame. Try:\n"+
                        "df = pandas.read_excel('file.xlsx')\n"+
                        "df = riversand.import_data(df)")
        
    # fold cases
    keep_cases =  ['Be', 'Be-10', '10Be', 
                   'Al', 'Al-26', '26Al',
                   'N', 'dN', 'delN', 'N_uncert']
    # also uncertainties can be 'Be_uncert' or 'Be_uncert*' (case-insensitive)
    # standardizations can be 'Be_standard' or 'Be_standard*' or 'Be_std' (case-insensitive)
    # note that 'standardisation', 'standardization', 'standard', 'std' are accepted for 'N'
    renames = {}
    for col in temp.columns:
        if col[:2] in ['N_']:
            if col[2:8].casefold()=='uncert':
                renames[col] = col[:2]+'uncert' # uncertainties 'N_uncert...'
        elif col[:3] in ['Be_', 'Al_']:
            if col[3:9].casefold()=='uncert': 
                renames[col] = col[:3]+'uncert' # uncertainties 'Be_uncert...', 'Al_uncert..'
            if col[3:11].casefold()=='standard':
                renames[col] = col[:3]+'standardization' # standardizations 'Be_standard..', 'Al_standard.'
            if col[3:].casefold()=='std':
                renames[col] = col[:3]+'standardization' # standardizations 'Be_std', 'Al_std'
        elif col not in keep_cases:
            renames[col] = col.casefold() # casefold all other column names
    temp = temp.rename(columns=renames)
    if len(temp.columns.unique()) < len(temp.columns):
        print("Non-unique column names.")
        return dfS, dfN
    
            
    # Check that all SAMPLE SITE data are provided; rename columns if necessary
    # name
    if 'name' not in temp.columns:
        print("Missing column 'name'.")
        return dfS, dfN
    else:
        temp['name'] = temp['name'].astype(str)
        if not all([valid_namestring(nm) for nm in temp['name'].to_numpy().tolist()]):
            print("Invalid sample name.")
            return dfS, dfN

    # latitude
    keys = ['latitude', 'lat']
    rename = {k:'lat' for k in keys if k in temp.columns}
    if len(rename)==1:
        temp = temp.rename(columns=rename)
    elif len(rename)==0:
        print("Missing column 'latitude'.")
        return dfS, dfN
    elif len(rename) > 1:
        print("Multiple columns may be interpreted as 'latitude'.")
        return dfS, dfN
    temp['lat'] = temp['lat'].fillna('nan') # force ValueError for any missing value
    try:
        temp['lat'] = pd.to_numeric(temp['lat'])
    except:
        print("Non-numeric or missing value in 'latitude'.")
        return dfS, dfN
    if ((temp.loc[:,'lat']<(-90)).any() & (temp.loc[:,'lat']>90).any()):
        print("Latitude is out of bounds (-90..+90).")
        return dfS, dfN
    
    # longitude
    keys = ['longitude', 'long', 'lon']
    rename = {k:'long' for k in keys if k in temp.columns}
    if len(rename)==1:
        temp = temp.rename(columns=rename)
    elif len(rename)==0:
        print("Missing column 'longitude'.")
        return dfS, dfN
    elif len(rename) > 1:
        print("Multiple columns may be interpreted as 'longitude'.")
        return dfS, dfN
    temp['long'] = temp['long'].fillna('nan') # force ValueError for any missing value
    try:
        temp['long'] = pd.to_numeric(temp['long'])
    except:
        print("Non-numeric or missing value in 'longitude'.")
        return dfS, dfN
    if ((temp.loc[:,'long']<(-180)).any() | (temp.loc[:,'long']>180).any()):
        print("Longitude is out of bounds (-180..+180).")
        return dfS, dfN
    
    # elevation
    keys = ['elevation', 'elev', 'pressure', 'press']
    rename = {k:'elevation' for k in keys if k in temp.columns}
    if len(rename)==1:
        temp = temp.rename(columns=rename)
    elif len(rename)==0:
        print("Missing column 'elevation' or 'pressure'.")
        return dfS, dfN
    elif len(rename) > 1:
        print("Multiple columns may be interpreted as 'elevation' or 'pressure'.")
        return dfS, dfN
    temp['elevation'] = temp['elevation'].fillna('nan') # force ValueError for any missing value
    try:
        temp['elevation'] = pd.to_numeric(temp['elevation'])
    except:
        print("Non-numeric or missing value in 'elevation' or 'pressure'.")
        return dfS, dfN
    if (temp.loc[:,'elevation']<(-500)).any():
        print("Elevation (pressure) is less than -500.")
        return dfS, dfN
    
    #press_flag
    if 'press_flag' in temp.columns:
        if ~temp.loc[:,'press_flag'].isin(['std', 'ant', 'pre']).all():
            print("Invalid 'press_flag'; must be 'std', 'ant' or 'pre'.")
            return dfS, dfN
    else:
        print("No column 'press_flag', assuming default value '{}'."
             .format(Params.default_values['press_flag']))
        temp.loc[:,'press_flag'] = Params.default_values['press_flag']
    
    keys = ['thickness', 'thick']
    rename = {k:'thickness' for k in keys if k in temp.columns}
    if len(rename)==1:
        temp = temp.rename(columns=rename)
    elif len(rename)==0:
        print("No column 'thickness', assuming default value {} cm."
              .format(Params.default_values['thickness']))
        temp.loc[:,'thickness'] = Params.default_values['thickness']
    elif len(rename) > 1:
        print("Multiple columns may be interpreted as 'thickness'.")
        return dfS, dfN
    temp['thickness'] = temp['thickness'].fillna('nan') # force ValueError for any missing value
    try:
        temp['thickness'] = pd.to_numeric(temp['thickness'])
    except:
        print("Non-numeric or missing value in 'thickness'.")
        return dfS, dfN
    if (temp.loc[:,'thickness']<0).any():
        print("Sample thickness is less than 0.")
        return dfS, dfN
    
    keys = ['density', 'rho']
    rename = {k:'density' for k in keys if k in temp.columns}
    if len(rename)==1:
        temp = temp.rename(columns=rename)
    elif len(rename)==0:
        print("No column 'density', assuming default value {} g/cm3."
              .format(Params.default_values['density']))
        temp.loc[:,'density'] = Params.default_values['density']
    elif len(rename) > 1:
        print("Multiple columns may be interpreted as 'density'.")
        return dfS, dfN
    temp['density'] = temp['density'].fillna('nan') # force ValueError for any missing value
    try:
        temp['density'] = pd.to_numeric(temp['density'])
    except:
        print("Non-numeric or missing value in 'density'.")
        return dfS, dfN
    if (temp.loc[:,'density']<=0).any():
        print("Density is 0 or less.")
        return dfS, dfN
    
    if 'shielding' in temp.columns:
        temp['shielding'] = temp['shielding'].fillna('nan') # force ValueError for any missing value
        try:
            temp['shielding'] = pd.to_numeric(temp['shielding'])
        except:
            print("Non-numeric or missing value in 'shielding'.")
            return dfS, dfN
    else:
        print("No column 'shielding', assuming default value {:.1f}."
              .format(Params.default_values['shielding']))
        temp.loc[:,'shielding'] = Params.default_values['shielding']
    if ((temp.loc[:,'shielding']<0).any() | (temp.loc[:,'shielding']>1).any()):
        print("Shielding is out of bounds (0..1).")
        return dfS, dfN

    keys = ['erate', 'erosion_rate']
    rename = {k:'erate' for k in keys if k in temp.columns}
    if len(rename)==1:
        temp = temp.rename(columns=rename)
    elif len(rename)==0:
        print("No column 'erate', assuming default value {} cm/yr."
              .format(Params.default_values['erate']))
        temp.loc[:,'erate'] = Params.default_values['erate']
    elif len(rename) > 1:
        print("Multiple columns may be interpreted as erosion rate ('erate').")
        return dfS, dfN
    temp['erate'] = temp['erate'].fillna('nan') # force ValueError for any missing value
    try:
        temp['erate'] = pd.to_numeric(temp['erate'])
    except:
        print("Non-numeric or missing value in 'erate'.")
        return dfS, dfN
    if (temp.loc[:,'erate']<0).any():
        print("Erosion rate is less than 0.")
        return dfS, dfN
    
    if 'year' in temp.columns:
        temp['year'] = temp['year'].fillna('nan') # force ValueError for any missing value
        try:
            temp['year'] = pd.to_numeric(temp['year'], downcast='integer')
        except:
            print("Non-numeric or missing value in 'year'.")
            return dfS, dfN
    else:
        print("No column 'year', assuming default value {}."
              .format(Params.default_values['year']))
        temp.loc[:,'year'] = Params.default_values['year']
    if ((temp.loc[:,'year']<1700).any() | (temp.loc[:,'year']>2100).any()):
        print("Year of collection is out of bounds (1700..2100).")
        return dfS, dfN
    
    dfS = temp[dfS.keys()] # copy sample site columns to dfS
    dfS = dfS[~dfS.duplicated()] # drop rows with identical sample site information
    if dfS['name'].duplicated().any(): # raise error if duplicate names remain
        print("Duplicate samples with non-matching site information (samples: {})."
              .format(', '.join(set(dfS.loc[dfS['name'].duplicated(), 'name']))))
        return dfS, dfN
    
    
    
    
    # if there is 'N' and 'delN' then 'nuclide' defaults to 'Be-10'; 'mineral' is optional
    # or else look for 'Be' and 'Be_uncert' and 'Al' and 'Al_uncert'

    def validate_standardization(temp):
        # raise exception if any standardization is invalid
        col = temp.loc[:,'nuclide']=='Be-10'
        #invalid_std = [v for v in temp.loc[col,'standardization'].values if v not in Params.Be_stds.keys()] #list of invalid standard names
        invalid_std = [i for i,v in temp.loc[col,'standardization'].items() if v not in Params.Be_stds.keys()] #list of indices with invalid standards
        if len(invalid_std)>0:
            raise ValueError("Invalid standardization for Be-10 data (lines: {})."
                            .format(', '.join(str(v) for v in invalid_std)))
        col = temp.loc[:,'nuclide']=='Al-26'
        #invalid_std = [v for v in temp.loc[col,'standardization'].values if v not in Params.Al_stds.keys()]
        invalid_std = [i for i,v in temp.loc[col,'standardization'].items() if v not in Params.Al_stds.keys()]
        if len(invalid_std)>0:
            raise ValueError("Invalid standardization for Al-26 data (lines: {})."
                            .format(', '.join(str(v) for v in invalid_std)))
            
    # Check that NUCLIDE data are provided; rename columns

    # rename nuclide columns to 'Be' and 'Al'; raise exception if 'N' as well as 'Be','Al' are specified
    keys = ['Be', 'Be-10', '10Be']
    rename = {k:'Be' for k in keys if k in temp.columns}
    if len(rename)>1:
        print("Multiple columns may be interpreted as Be concentration.")
        return dfS, dfN
    elif (len(rename)==1) and ('N' in temp.columns):
        print("Multiple columns ('N','Be') may be interpreted as Be concentration.")
        return dfS, dfN
    elif len(rename)==1:
        temp = temp.rename(columns=rename)
        
    keys = ['Be_standardization', 'Be_standardisation', 'Be_standard', 'Be_std']
    rename = {k:'Be_standardization' for k in keys if k in temp.columns}
    if len(rename)>1:
        print("Multiple columns may be interpreted as 'Be_standardization'.")
        return dfS, dfN
    elif len(rename)==1:
        temp = temp.rename(columns=rename)

    keys = ['Al', 'Al-26', '26Al']
    rename = {k:'Al' for k in keys if k in temp.columns}
    if len(rename)>1:
        print("Multiple columns may be interpreted as Al concentration.")
        return dfS, dfN
    elif (len(rename)==1) and ('N' in temp.columns):
        print("Multiple columns ('N','Al') may be interpreted as Al concentration.")
        return dfS, dfN
    elif len(rename)==1:
        temp = temp.rename(columns=rename)
        
    keys = ['Al_standardization', 'Al_standardisation', 'Al_standard', 'Al_std']
    rename = {k:'Al_standardization' for k in keys if k in temp.columns}
    if len(rename) > 1:
        print("Multiple columns may be interpreted as 'Al_standardization'.")
        return dfS, dfN
    elif len(rename)==1:
        temp = temp.rename(columns=rename)


    if 'N' in temp.columns:
        df1 = pd.DataFrame()
        try:
            #temp['N'] = pd.to_numeric(temp['N'])
            df1 = temp.loc[pd.to_numeric(temp.loc[:, 'N'])>0,:]
        except:
            print("Non-numeric value in nuclide concentration 'N'.")
            return dfS, dfN
        
        # N, delN
        keys = ['delN', 'dN', 'N_uncert']
        rename = {k:'delN' for k in keys if k in temp.columns}
        if len(rename)==1:
            df1 = df1.rename(columns=rename)
        elif len(rename)==0:
            print("Missing uncertainty 'delN' for concentration 'N'.")
            return dfS, dfN
        elif len(rename) > 1:
            print("Multiple columns may be interpreted as uncertainty 'delN'.")
            return dfS, dfN
        df1['delN'] = df1['delN'].fillna('nan') # force ValueError for any missing value
        try:
            df1['delN'] = pd.to_numeric(df1['delN'])
        except:
            print("Non-numeric or missing value in nuclide uncertainty 'delN'.")
            return dfS, dfN
        if (df1.loc[:,'delN']<=0).any():
            print("Nuclide uncertainty must be >0.")
            return dfS, dfN
        
        # nuclide
        if 'nuclide' in df1.columns:
            if ~df1.loc[:,'nuclide'].isin(['Be-10','Al-26']).all():
                print("Invalid nuclide; must be 'Be-10' or 'Al-26'.")
                return dfS, dfN
        else:
            print("No column 'nuclide', assuming default '{}'."
                  .format(Params.default_values['nuclide']))
            df1.loc[:,'nuclide'] = Params.default_values['nuclide']

        # mineral
        if 'mineral' in df1.columns:
            if ~df1.loc[:,'mineral'].isin(['quartz']).all():
                print("Invalid mineral; must be 'quartz'.")
                return dfS, dfN
        else:
            #print("No column 'mineral', assuming default mineral '{}'."
            #  .format(Params.default_values['mineral']))
            df1.loc[:,'mineral'] = Params.default_values['mineral']

        # standardization
        keys = ['standardisation', 'standardization', 'standard', 'std']
        rename = {k:'standardization' for k in keys if k in df1.columns}
        if len(rename)==1:
            df1 = df1.rename(columns=rename)
            try:
                validate_standardization(df1)
            except ValueError as e:
                print(e)
                return dfS, dfN
        elif len(rename)==0:
            print("No column 'standardization', assuming default values ({})."
                  .format(', '.join(["{}: {}".format(k,v) for k,v in Params._default_standards.items()])))
            for k,v in Params._default_standards.items():
                col = (df1.loc[:,'mineral']=='quartz') & (df1.loc[:,'nuclide']==k)
                df1.loc[col,'standardization'] = v
        elif len(rename) > 1:
            print("Multiple columns may be interpreted as 'standardization'.")
            return dfS, dfN
        dfN = df1.loc[:, dfN.columns] # df1 has only nuclide-relevant columns and rows
        
        
    elif ('Be' in temp.columns) or ('Al' in temp.columns): # 'mineral' and 'nuclide' are ignored
        dfBe = pd.DataFrame()
        dfAl = pd.DataFrame()
        
        if 'Be' in temp.columns:
            try:
                dfBe = temp.loc[pd.to_numeric(temp.loc[:, 'Be'])>0,:]
            except:
                print("Non-numeric value in Be concentration.")
                return dfS, dfN
            dfBe.loc[:,'mineral'] = 'quartz'
            dfBe.loc[:,'nuclide'] = 'Be-10'
            dfBe.loc[:,'N'] = dfBe.loc[:,'Be']
            try:
                dfBe.loc[:,'delN'] = dfBe.loc[:,'Be_uncert']
            except:
                print("Missing uncertainty 'Be_uncert'.")
                return dfS, dfN
            dfBe['delN'] = dfBe['delN'].fillna('nan') # force ValueError for any missing value
            try:
                dfBe['delN'] = pd.to_numeric(dfBe['delN'])
            except:
                print("Non-numeric or missing value in Be uncertainty.")
                return dfS, dfN
            if (dfBe.loc[:,'delN']<=0).any():
                print("Be uncertainty must be >0.")
                return dfS, dfN
            try:
                dfBe.loc[:,'standardization'] = dfBe.loc[:,'Be_standardization']
            except:
                print("No column 'Be_standardization', assuming default values (Be-10: {})."
                      .format(Params._default_standards['Be-10']))
                dfBe.loc[:,'standardization'] = Params._default_standards['Be-10']
            try:
                validate_standardization(dfBe)
            except ValueError as e:
                print(e)
                return dfS, dfN
            dfBe.index = dfBe.index*2
            
        if 'Al' in temp.columns:
            try:
                dfAl = temp.loc[pd.to_numeric(temp.loc[:, 'Al'])>0,:]
            except:
                print("Non-numeric value in Al concentration.")
                return dfS, dfN
            dfAl.loc[:,'mineral'] = 'quartz'
            dfAl.loc[:,'nuclide'] = 'Al-26'
            dfAl.loc[:,'N'] = dfAl.loc[:,'Al']
            try:
                dfAl.loc[:,'delN'] = dfAl.loc[:,'Al_uncert']
            except:
                print("missing uncertainty 'Al_uncert'.")
                return dfS, dfN
            dfAl['delN'] = dfAl['delN'].fillna('nan') # force ValueError for any missing value
            try:
                dfAl['delN'] = pd.to_numeric(dfAl['delN'])
            except:
                print("Non-numeric or missing value in Al uncertainty.")
                return dfS, dfN
            if (dfAl.loc[:,'delN']<=0).any():
                print("Al uncertainty must be >0.")
                return dfS, dfN
            try:
                dfAl.loc[:,'standardization'] = dfAl.loc[:,'Al_standardization']
            except:
                print("No column 'Al_standardization', assuming default values (Al-26: {})."
                      .format(Params._default_standards['Al-26']))
                dfAl.loc[:,'standardization'] = Params._default_standards['Al-26']
            try:
                validate_standardization(dfAl)
            except ValueError as e:
                print(e)
                return dfS, dfN
            dfAl.index = dfAl.index*2 + 1
            
        
        df1 = pd.concat([dfBe, dfAl], axis=0)#, ignore_index=True)
        dfN = df1.loc[:, dfN.columns] # dfN has only nuclide-relevant columns and rows
        dfN = dfN.sort_index()
    
    else:
        print("Missing columns for nuclide concentrations: either 'N','delN','standardization' "+
              "or 'Be','Be_uncert','Be_std','Al','Al_uncert','Al_std'")
        return dfS, dfN
    
    
    checksum = set(dfN.loc[:,'name']) - set(dfS.loc[:,'name'])
    if len(checksum)>0:
        print("No sample sites info for samples: {} (this error should never occur!!)"
              .format(','.join(sorted(checksum))))
        dfN = dfN[dfN.name.isin(dfS.name)] # drop rows
    checksum = set(dfS.loc[:,'name']) - set(dfN.loc[:,'name'])
    if len(checksum)>0:
        print("No nuclide data for samples: {}."
              .format(', '.join(sorted(checksum))))    
        dfS = dfS[dfS.name.isin(dfN.name)] # drop rows
        
    return dfS, dfN
    
    
    
def get_textblock(df) -> str:
    """
    Get textblock for online calculator. Relies on import_data() for validation.
    Returns an empty string if import_data() encounters problems with the input
    and returns an empty DataFrame.
    
    Parameters
    ----------
    df : pd.DataFrame
        > df = pd.read_excel('./folder/file.ods')   # read from excel spreadsheet
        > df = import_data(df)                      # recast for online calculator / validate
        > textblock = get_textblock(df)             # convert data to string
        or
        > df = pd.read_excel('./folder/file.ods')   # read from excel spreadsheet
        > textblock = get_textblock(df)
                
    Returns
    -------
    textblock : str
        Textblock for online calculator v.3
        
    Raises
    ------
    ValueError
        Raised by validate_topo(), validate_sample(), validate_nuclide()
        e.g. if mandatory item is missing or ill-defined.
        This should never happen because get_textblock() validates the input
        using import_data2().

    """

    
    dfS, dfN = import_data2(df)
        
    dfS = dfS.to_dict(orient='records') # list of dicts, sample site information
    dfN = dfN.to_dict(orient='records') # list of dicts, nuclide information
    
    #import_data2() returns empty df's if data validation failed
    if len(dfS)==0: return ''
    if len(dfN)==0: return ''
    
    # the validation should not be necessary; data is validated by import_data2()
    def lineS(d:dict): # sample and topo data for each row in dfS
        try:
            out1 = validate_topo(d, shielding=True) # 'lat', 'long', 'elevation', 'shielding'
            out3 = validate_sample(d) # 'press_flag', 'density', 'thickness', 'erate', 'year'
        except ValueError as e:
            raise e
        out = {**out1, **out3}    
            
        textline = "{} {:.5f} {:.5f} {:.3f} {} {} {} {:.5f} {:.5f} {};".format(
            d['name'], out['lat'], out['long'], out['elevation'], out['press_flag'],
            out['thickness'], out['density'], out['shielding'], out['erate'], out['year'])
        return textline
    
    # the validation should not be necessary; data is validated by import_data2()
    def lineN(d:dict): # nuclide data for each row in dfN
        try:
            out = validate_nuclide(d) # 'name', 'N', 'delN', 'nuclide', 'mineral', 'standardization'
        except ValueError as e:
            raise e
            
        textline = "{} {} {} {} {} {};".format(
            out['name'], out['nuclide'], out['mineral'], out['N'], out['delN'], out['standardization'])
        return textline
    
    textblockS = [lineS(d) for d in dfS]
    textblockN = [lineN(d) for d in dfN]
    textblockS = ' '.join(textblockS)
    textblockN = ' '.join(textblockN)
    return textblockS + ' ' + textblockN