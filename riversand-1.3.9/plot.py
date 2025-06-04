"""

*******************************************************************************
plot.py  :  plotting functions

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

Plotting functions for the Riversand project

plot_raster(R) - plot rasterio (Raster or Riversand object)
plot_clipped_raster(clips) - plot clipped catchment (dict of xr.DataArray)
plot_polyfit(E, delE, NofE, sample_data)

"""


from scipy import optimize
import numpy as np

import rasterio
from rasterio.plot import show

import matplotlib.pyplot as plt

from riversand.params import Params

grey_maps = {'elevation'  : 'gray',
             'shielding'  : 'gray',
             'snow'       : 'gray',
             'vegetation' : 'gray',
             'quartz'     : 'binary'
             }

colormaps = {'elevation'  : 'terrain',
             'shielding'  : 'Reds_r',
             'snow'       : 'Blues_r',
             'vegetation' : 'Greens_r',
             'quartz'     : 'binary'
             }


def plot_raster(R, dtype='elevation',
                cmap=None, vmin=None, vmax=None):
    """
    Plot raster datasets. Returns matplotlib.pyplot Figure and Axis objects,
    which can be modified and/or saved.

    > fig, ax = riversand.plot_raster(rv.shielding)
    > fig.savefig('filename.png')
    
    Parameters
    ----------
    R : Riversand object (R=rv) or Raster object (R=rv.elevation, rv.shielding,...)
    
    dtype : str, optional
        If R=rv, specifies, which raster to plot ('elevation', 'shielding',
        'snow', 'vegetation', 'quartz'). The default is 'elevation'.

    cmap : str, optional
        Registered colormap name used to map scalar data to colors.
        The default is None (automatic colormap).
        
    vmin, vmax : float, optional
        Definition of the data range covered by the colormap.
        The default is None (complete value range; vmax=1 for shielding).
        
    Returns
    -------
    fig : Figure
    ax : Axes

    """
    from riversand.geospatial import Riversand, Raster

    fig, ax = plt.subplots()

    # called as plot_raster(rv, dtype='shielding'):
    if isinstance(R, Riversand):
        if dtype=='elevation':
            R = R.elevation
        elif dtype=='shielding':
            R = R.shielding
        elif dtype=='snow':
            R = R.snow
        elif dtype=='vegetation':
            R = R.vegetation
        elif dtype=='quartz':
            R = R.quartz
        else:
            print("Invalid dtype '{}'".format(dtype))
            plt.close()
            return fig, ax
        
    if R is None:
        print("No raster '{}' available".format(dtype))
        plt.close()
        return fig, ax
    
    alpha = 1
    # called as plot_raster(rv.shielding):
    if isinstance(R, Raster):
        dtype = R.dtype
        if dtype=='elevation':
            with rasterio.open(R.fname, 'r') as src:
                ar = src.read()
                if vmin is None: vmin = np.nanmin(ar)
                if vmax is None: vmax = np.nanmax(ar)
                #try:
                #    vmin = np.nanmin(ar[ar!=src.nodata])
                #except:
                #    vmin = np.nanmin(ar)
                
        if dtype in {'shielding', 'snow', 'vegetation'}:
            with rasterio.open(R.fname, 'r') as src:
                ar = src.read()
                if vmin is None: vmin = max([np.nanmin(ar), 0])
                if vmax is None: vmax = max([np.nanmax(ar), 1])
                
        if dtype=='quartz':
            with rasterio.open(R.fname, 'r') as src:
                ar = src.read()
                if vmin is None: vmin = 0
                if vmax is None: vmax = 1
                alpha = 0.5

        if cmap is None:
            cmap = colormaps[dtype] # set to default colormaps
        ax = show(ar, transform=src.transform, cmap=cmap,
                  vmin=vmin, vmax=vmax, alpha=alpha,
                  ax=ax)
    
    plt.xticks(fontsize=6)#, ax.set_xlabel('x', fontsize=12)
    plt.yticks(fontsize=6)#, ax.set_ylabel('y', fontsize=12)
    ax.set_aspect('equal')
    
    im = ax.get_images()[0]
    fig.colorbar(im, shrink=.7, label=dtype)
        
    plt.close()
    
    return fig, ax
    


def plot_clipped_raster(clips, dtype='elevation',
                        cmap=None, vmin=None, vmax=None):
    """
    Plot the clipped catchment. Returns matplotlib.pyplot Figure and Axis
    objects, which can be modified and/or saved.

    > fig, ax = riversand.plot_clipped_raster(clips, dtype='shielding')
    > fig.savefig('filename.png')
    
    Parameters
    ----------
    clips : dict of xr.DaraArray
        Keys are 'elevation', 'shielding', 'snow', 'vegetation', 'quartz'
        as well as 'name' (sample name) and 'epsg'.
   
    dtype : str, optional
        Specifies, which raster to plot ('elevation', 'shielding',
        'snow', 'vegetation', 'quartz'). The default is 'elevation'.
   
    cmap : str, optional
        Registered colormap name used to map scalar data to colors.
        The default is None (automatic colormap).
        
    vmin, vmax : float, optional
        Definition of the data range covered by the colormap.
        The default is None (complete value range; vmax=1 for shielding).
        
    Returns
    -------
    fig : Figure
    ax : Axes

    """
    
    if dtype not in clips.keys():
        raise ValueError("No raster '{}' in clips".format(dtype))
        
    fig, ax = plt.subplots()
    
    alpha = 1
    X = clips[dtype]
    if dtype=='elevation':
        if vmin is None: vmin = np.nanmin(X)
        if vmax is None: vmax = np.nanmax(X)
    
    if dtype in {'shielding', 'snow', 'vegetation'}:
        if vmin is None: vmin = max([np.nanmin(X), 0])
        if vmax is None: vmax = max([np.nanmax(X), 1])
    
    if dtype=='quartz':
        if vmin is None: vmin = 0
        if vmax is None: vmax = 1
        alpha = 0.5
    
    if cmap is None:
        cmap = colormaps[dtype] # set to default colormaps
    X.squeeze().plot.imshow(ax=ax, cmap=cmap, vmin=vmin, vmax=vmax, alpha=alpha,
                            cbar_kwargs={'label' : dtype, 'shrink' : .7})
    
        
    #outline.plot(edgecolor='red', facecolor='none', ax=ax[i])
    #outline.centroid.plot(ax=ax[i], color='red')
    plt.xticks(fontsize=6)#, ax.set_xlabel('x', fontsize=12)
    plt.yticks(fontsize=6)#, ax.set_ylabel('y', fontsize=12)
    ax.set_aspect('equal')
    ax.set_title(dtype)
    
    plt.close()
    
    return fig, ax


def plot_polyfit(E, delE, NofE, sample_data, unit='mm/yr',
                 **kwargs):
    """
    Plot nuclide concentration vs. erosion rate, data and polynomial fit.

    Parameters
    ----------
    E : number
        Calculated erosion rate in cm/yr.
        
    delE : tuple
        Uncertainty un the calculated erosion rate (delE-, delE+) in cm/yr.
    
    NofE : pd.Series
        Nuclide concentration as function of erosion rate.
        NofE.name indicates the scaling method.
        
        > E, delE, NofE, RMSE, err = riversand.poly_E_results(
                sample_data, topostats, shielding, erates, scaling)
            
    sample_data : dict, pd.Series or pd.DataFrame
        Sample and nuclide information of a single sample. Sample data will be
        re-standardized to 07KNSTD (Be) or KNSTD (Al)
        
        > sample_data = rv.samples.iloc[3] # sample at table row number 3
        > sample_data = rv.samples.loc[rv.samples['name']=='DB05'] # sample with the name 'DB05'
        Note that the second option only works if there is exactly one sample with the name 'DB05'

    unit : str, optional
        unit for plotting erosion rates, for example 'mm/yr' or 'cm/Ma'.
        The default is 'mm/yr'.
        
    Use optional kwarg 'linfit' to plot a linear function.
    
    Returns
    -------
    fig : Figure
    ax : Axes

    """
    
    if 'linfit' in kwargs.keys():
        from riversand.calc import NofE_linfitfunc as func
    else:
        from riversand import NofE_fitfunc as func
    from riversand.utils import restandardize
    
    units = Params.units
    
    if unit not in units.keys():
        print("'{}' is not a valid unit; default to 'mm/yr'".format(unit))
        unit = 'mm/yr'
    
    try:
        sample_data = restandardize(sample_data)
    except:
        pass#raise Warning("Cannot restandardize sample data")
            
        

    fig, ax = plt.subplots()
    
    try:
        ax.plot(NofE.index*units[unit], NofE, 'ok')
    except:
        pass#raise Warning("Cannot plot NofE")
        
    try:
        popt, pcov = optimize.curve_fit(func, NofE.index, NofE)
        x = np.linspace(NofE.index[0], NofE.index[-1], 100)
        y = func(x, *popt)
        ax.plot(x*units[unit], y, '-r')
    except:
        pass#raise Warning("Cannot plot a curve fit to NofE")
        
    if not np.isnan(E) and not np.isnan(delE).any():
        plt.errorbar(E*units[unit], sample_data['N'],
                 yerr=sample_data['delN'],
                 xerr=[[e*units[unit]] for e in delE],
                 capsize=5, fmt='or')
    elif not np.isnan(E):
        plt.errorbar(E*units[unit], sample_data['N'],
                 yerr=sample_data['delN'],
                 #xerr=[[e*units[unit]] for e in delE],
                 capsize=5, fmt='or')
    else:
        pass#raise Warning("Cannot plot a curve fit to NofE")

    if 'linfit' in kwargs.keys():
        ax.text(0.95, 0.95, 'linear fit!',
                transform=ax.transAxes, fontsize=14,
                verticalalignment='top', horizontalalignment='right',
                backgroundcolor='white', color='r')#, bbox=props)
        
    plt.xlabel('Erosion rate E, {}'.format(unit))
    plt.ylabel('Nuclide concentration N, at/g')
        
    plt.title("{}: E = {:.2f} +{:.2f}/-{:.2f} {}".format(
        NofE.name, E*units[unit],
        delE[1]*units[unit], delE[0]*units[unit], unit))
        
    plt.grid()
    
    plt.close()
    
    return fig, ax
    
