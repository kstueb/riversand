"""
Created on Sun Feb 19 14:17:07 2023

***** plot.py ***********************************************************

Plotting functions for the Riversand project

plot_raster(R) - plot rasterio (Raster or Riversand object)
plot_clipped_raster(clips) - plot clipped catchment (dict of xr.DataArray)
plot_polyfit(E, delE, NofE, sample_data)


@author: Konstanze Stübner, kstueb@gmail.com

"""

import os
from scipy import optimize
import numpy as np

import rasterio
from rasterio.plot import show

import matplotlib.pyplot as plt

def plot_raster(R, dtype='elevation', fname='show'):
    """
    Plot raster datasets.

    Parameters
    ----------
    R : Riversand object (R=rv) or Raster object (R=rv.elevation, rv.shielding,...)
    dtype : str, optional
        If R=rv, specify which raster to plot, e.g., 'elevation', 'shielding',
        'quartz'. The default is 'elevation'.
    fname : str, optional
        'show' : show the figure inline, do not save.
        'auto' : save figure as .jpg with auto-generated file name.
        'jpg'  : same as 'auto'.
        'png'  : same as 'auto' but save as .png.
        Custom filename including extension, e.g. 'my_image.tif'
        The default is 'auto'.

    Returns
    -------
    None.

    """
    from riversand import params
    from riversand.geospatial import Riversand, Raster

    fig, ax = plt.subplots()
    
    if isinstance(R, Raster):
        label = R.dtype
        with rasterio.open(R.fname, 'r') as src:
            ax = show(src.read(), transform=src.transform, cmap='gray', ax=ax)
    
    elif isinstance(R, Riversand):
        label = dtype
        if dtype=='elevation':
            with rasterio.open(R.elevation.fname, 'r') as src:
                ax = show(src.read(), transform=src.transform, cmap='gray', ax=ax)
        if dtype=='shielding':
            with rasterio.open(R.shielding.fname, 'r') as src:
                ax = show(src.read(), transform=src.transform, cmap='gray', ax=ax)
        if dtype=='quartz':
            with rasterio.open(R.quartz.fname, 'r') as src:
                ax = show(src.read(), transform=src.transform, cmap='gray', alpha=0.5, ax=ax)
    else:
        plt.close()
        return
    
   
    plt.xticks(fontsize=10), ax.set_xlabel('x', fontsize=12)
    plt.yticks(fontsize=10), ax.set_ylabel('y', fontsize=12)
    
    
    if fname=='show':
        plt.show()
        return
    
    fullname = "_{}{}".format(label[0].upper(), label[1:].lower())
    if fname in {'jpg', 'png'}:
        fullname = "{}.{}".format(fullname, fname)
    else:     
        fullname = "{}.jpg".format(fullname)
            
    plt.savefig(os.path.join(params.out_path, fullname), bbox_inches='tight', transparent=True)
    plt.close()
    
    
def plot_clipped_raster(clips, c_name='', label='elevation', fname='auto'):
    """
    Plot the clipped catchment.

    Parameters
    ----------
    clips : dict of xr.DaraArray with keys 'elevation', 'shielding', 'quartz'
        
    c_name : str, optional
        catchment name for saving the figure. The default is ''.
   
    label : str, optional
        selects, which raster to plot ('elevation', 'shielding', 'quartz').
        The default is 'elevation'.
   
    fname : str, optional
        'show' : show the figure inline, do not save.
        'auto' : save figure as .jpg with auto-generated file name.
        'jpg'  : same as 'auto'.
        'png'  : same as 'auto' but save as .png.
        Custom filename including extension, e.g. 'my_image.tif'
        The default is 'auto'.

    Returns
    -------
    None.

    """
    from riversand import params
    
    if label not in clips.keys():
        raise ValueError("No '{}' raster in clips".format(label))
        
    X = clips[label]
    fig, ax = plt.subplots()
    X.squeeze().plot.imshow(ax=ax, cmap='gray', cbar_kwargs={'label' : label})
    #outline.plot(edgecolor='red', facecolor='none', ax=ax[i])
    #outline.centroid.plot(ax=ax[i], color='red')
    plt.xticks(fontsize=10), ax.set_xlabel('x', fontsize=12)
    plt.yticks(fontsize=10), ax.set_ylabel('y', fontsize=12)
    ax.set_aspect('equal')
    ax.set_title(label)
    
    if fname=='show':
        plt.show()
        return
    
    if c_name=='':
        fullname = label
    else:
        fullname = "{}_{}".format(c_name, label)
        
    if fname in {'jpg', 'png'}:
        fullname = "{}.{}".format(fullname, fname)
    else:     
        fullname = "{}.jpg".format(fullname)
            
    plt.savefig(os.path.join(params.out_path, fullname), bbox_inches='tight', transparent=True)
    plt.close()
    

def plot_polyfit(E, delE, NofE, sample_data,
                 unit='mm/yr', fname='auto',
                 **kwargs):
    """
    Plot nuclide concentration vs. erosion rate, data and polynomial fit.

    Parameters
    ----------
    E : number
        calculated erosion rate in cm/yr.
        
    delE : tuple
        uncertainty un the calculated erosion rate (delE-, delE+) in cm/yr.
    
    NofE : pd.Series
        nuclide concentration as function of erosion rate.
        NofE.name indicates the scaling method.
        
        for exmaple:
            > NofE_all = riversand.get_NofE(sample_data, topostats, shileding, erates)
            > NofE = NofE_all['LSDn']
        or:
            > E, delE, NofE, err = riversand.calc.poly_E_results(
                sample_data, topostats, shielding, scaling, erates)
            
    sample_data : dict, pd.Series or pd.DataFrame
        sample and nuclide information of a single sample.
        
        for example:
            > sample_data = rv.samples.iloc[3] # sample at table row number 3
        or:
            > sample_data = rv.samples.loc[rv.samples['name']=='DB05'] # sample with the name 'DB05'
        Note that the second option only works if there is exactly one sample with the name 'DB05'

    unit : str, optional
        unit for plotting erosion rates, for example 'mm/yr' or 'cm/Ma'.
        The default is 'mm/yr'.
        
    fname : str, optional
        'show' : show the figure inline, do not save.
        'auto' : save figure as .jpg with auto-generated file name.
        'jpg'  : same as 'auto'.
        'png'  : same as 'auto' but save as .png.
        Custom filename including extension, e.g. 'my_image.tif'
        The default is 'auto'.

    kwargs :
        catchment name 'c_name':str
        scaling method 'scaling':str
            auto filename: "{}_{}.jpg".format(c_name, scaling)
        full file name 'fullname':str
            e.g., "DB05_LSDn.jpg"
        
    Returns
    -------
    None.

    """
    from riversand import params
    from riversand.params import units
    from riversand import NofE_fitfunc as func
    if unit not in units.keys():
        print("'{}' is not a valid unit; default to 'mm/yr'".format(unit))
        unit = 'mm/yr'
    
    if 'scaling' in kwargs.keys():
        scaling = kwargs['scaling']
    else:
        try:
            scaling = NofE.name
        except:
            scaling = ''
            
    if 'c_name' in kwargs.keys():
        c_name = kwargs['c_name']
    else:
        try:
            c_name = sample_data['name']
        except:
            c_name = 'Test'
        
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

        
    plt.xlabel('Erosion rate E, {}'.format(unit))
    plt.ylabel('Nuclide concentration N, at/g')
        
    plt.title("{}: E = {:.2f} +{:.2f}/-{:.2f} {}".format(
        NofE.name, E*units[unit],
        delE[1]*units[unit], delE[0]*units[unit], unit))
        
    plt.grid()
    
    if fname=='show':
        plt.show()
        return

    if 'fullname' in kwargs.keys():
        fname = kwargs['fullname']
    elif fname in {'jpg', 'png'}:
        fname = "{}_{}.{}".format(c_name, scaling, fname)
    else:     
        fname = "{}_{}.jpg".format(c_name, scaling)
    plt.savefig(os.path.join(os.path.join(params.out_path, fname)), bbox_inches='tight', transparent=True)
    plt.close()
    
