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
    - Fill in default values from params.py if non-mandatory keys are missing
    - Raise ValueError if values are invalid / mandatory values are missing 

Geospatial processing functions:
    feature_in_raster()
    clip_raster()
    eliminate_quartzfree()
    get_xarray_centroid()
    get_polygon_centroid()
    projected_xy_to_longlat()
    get_bins()
    get_topostats()
    

"""

import numpy as np
import pandas as pd

import re
from numbers import Number
from copy import deepcopy

import rasterio
import rasterio.crs
import rasterio.mask
from rasterio import MemoryFile

import xarray as xr

# =============================================================================
# Validation for online calculator
# =============================================================================

def validate_topo(item, shielding:bool=False) -> dict:
    """
    Generate dict of topo data.
    mandatory keys 'lat', 'long', 'elevation'
    optional keys 'shielding'
    
    item : dict or pd.Series
    
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
    except:
        out['lat'] = np.nan
    if np.isnan(out['lat']):
        raise ValueError("Latitude is not a number")
    if not (-90 <= out['lat'] <= 90):
        raise ValueError("Latitude is out of bounds (-90..+90)")
    
    try:
        out['long'] = float(item['long'])
    except:
        out['long'] = np.nan
    if np.isnan(out['long']):
        raise ValueError("Longitude is not a number")
    if not (-180 <= out['long'] <= 180):
        raise ValueError("Longitude is out of bounds (-180..+180)")
    
    try:
        out['elevation'] = float(item['elevation'])
    except:
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
        except:
            out['shielding'] = np.nan
        if np.isnan(out['shielding']):
            raise ValueError("Shielding is not a number")
        if not (0 <= out['shielding'] <= 1):
            raise ValueError("Shielding is out of bounds (0..1)")
            
    return out


def validate_nuclide(item) -> dict:
    """
    Generate dict of nuclide data.
    mandatory keys 'N', 'delN'
    default values for 'name', 'nuclide', 'mineral', 'standardization'

    item : dict or pd.Series
    
    """
    
    from riversand import params
    
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
    except:
        raise ValueError("Concentration N is not a number") 
    if out['N'] <= 0:
        raise ValueError("Concentration N is 0 or less")
        
    try:
        out['delN'] = round(float(item['delN']))
    except:
        raise ValueError("Uncertainty delN is not a number") 
    if out['delN'] <= 0:
        raise ValueError("Uncertainty delN is 0 or less")
        
    # optional keys with default values:
    if 'name' in item.keys():
        name_regex = "^[a-zA-Z0-9_-]*$"
        try:
            out['name'] = str(item['name'])
        except:
            out['name'] = '??' # illegal character
        
        if len(out['name'])>32:
            raise ValueError("Sample name more than 32 characters")
        if not re.match(name_regex, out['name']):
            raise ValueError("Illegal characters in sample name")
    else:
        out['name'] = params.default_values['name'] # default value
    
    if 'nuclide' in item.keys():
        try:
            out['nuclide'] = str(item['nuclide'])
        except:
            out['nuclide'] = '??'
        if out['nuclide'] not in {'Be-10', 'Al-26'}:
            raise ValueError("Illegal nuclide")
    else:
        out['nuclide'] = params.default_values['nuclide'] # default value
        
    if 'mineral' in item.keys():
        try:
            out['mineral'] = str(item['mineral'])
        except:
            out['mineral'] = '??'
        if out['mineral'] != 'quartz':
            raise ValueError("Illegal mineral")
    else:
        out['mineral'] = params.default_values['mineral'] # default value

    standardization = params.default_standards[out['nuclide']]
    if 'standardization' in item.keys():
        try:
            out['standardization'] = str(item['standardization'])
        except:
            out['standardization'] = '??'
        if out['standardization'] != standardization:
            raise ValueError("Illegal standardization")
    else:
        out['standardization'] = standardization

    return out


def validate_sample(item, shielding:bool=False) -> dict:
    """
    Generate dict of sample data.
    default values for 'press_flag', 'density', 'thickness', 'erate', 'year'
    optional keys 'shielding'
    
    item : dict or pd.Series
    
    """
    
    from riversand import params
    
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
        out['press_flag'] = params.default_values['press_flag'] # default value
            
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
        out['density'] = params.default_values['density'] # default value
        
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
        out['thickness'] = params.default_values['thickness'] # default value
        
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
        out['erate'] = params.default_values['erate'] # default value

    if 'year' in item.keys():
        try:
            out['year'] = round(float(item['year']))
        except:
            raise ValueError("Collection year is not a number")
        if not (1700 <= out['year'] <= 2100):
            raise ValueError("Collection year is out of bounds (1700..2100)")
    else:
        out['year'] = params.default_values['year'] # default value
        
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


# =============================================================================
# Geospatial processing
# =============================================================================
    
def feature_in_raster(polygon:dict, src:rasterio.DatasetReader) -> bool:
    """
    Check whether feature is within the bounds of a raster.
    
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
    
    Raises ValueError if polygon is out of bounds.
    
    polygon : catchment.catchments[0]['geometry']
    label : 'elevation', 'shielding', 'quartz'
    
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
    
    clips : dict with keys 'quartz', 'elevation', 'shielding'.
            the key 'epsg' set by clip_all_rasters() is ignored
    
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
    X : projected raster dataset clipped to the catchment polygon.
        X = clips['elevation']
    
    Returns tuple of floats, projected longitude, latitude.
    
    """

    A = xr.DataArray((X.values>0), coords=X.coords)
    centr_x = float(np.sum(A.sum('y')/np.sum(A.sum('y')) * A.x))
    centr_y = float(np.sum(A.sum('x')/np.sum(A.sum('x')) * A.y))

    return (centr_x, centr_y)

    

def get_polygon_centroid(polygon:dict) -> tuple:
    """
    polygon : catchment.catchments[0]['geometry']
    
    Returns tuple of floats, projected longitude, latitude.
    
    """
    
    if 'coordinates' not in polygon.keys():
        raise TypeError("get_polygon_centroid() argument 'polygon' not a valid polygon dictionary")
        
    from shapely.geometry import Polygon
        
    centr_x = Polygon(polygon['coordinates'][0]).centroid.x
    centr_y = Polygon(polygon['coordinates'][0]).centroid.y

    return (centr_x, centr_y)

    
def projected_xy_to_longlat(xy:tuple, epsg:int) -> tuple:
    """    
    xy : tuple, projected coordinates
    epsg : integer epsg code
    
    Returns tuple of floats, unprojected longitude, latitude.
    
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
    
    Z : ndarray or xr.DataArray of the values to be binned.
    
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
    
    clips : dict of xr.DataArrays.
            keys are raster dtypes ('elevation', 'shielding, ...)
            values have attributes .transform (Affine) and .label (dtype)
            additional key 'epsg' has the projection
            clips = geospatial.clip_all_rasters(rv)
            
    centroid : "from_clipped" or tuple (long, lat) in projected coordinates.
    
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