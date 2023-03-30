#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

*******************************************************************************
calc.py  :  url requests to stoneage.hzdr.de and calculation of erosion rate 

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

General functions:
- get_textline() : Cast 'sample' and 'topo' data to textstring for the online
    calculator. Input dict ot pd.Series. If possible, missing values are filled
    with default values from params.py.
- get_erates_from_server() and get_NofE_from_server : Get pd.DataFrame,
    diagnosis and version for any string sent to the online calculator.

Specialized function:
- get_E() : dict w/ keys 'St', 'Lm', 'LSDn' for a single sample-nuclide pair.
- get_NofE_FullTable() : Full table of 'NpredXX' for 'topostats' and suite of
    'erates'. Results are weighted by 'wt' ('XX_wt').
- get_NofE() : Wrapper for _FullTable(): groups data by 'E_cmyr'.

The work horse:
- poly_E_results() : -> (E, delE, NofE, error).
    Returns error=['minE too high'] or ['maxE too low']
    or a list of one or several other errors.

"""

import numpy as np
import pandas as pd
from numbers import Number

from copy import deepcopy
from scipy import optimize

import urllib.parse
import urllib.request


def get_textline(sample:pd.Series, topo:pd.Series, shielding) -> str:
    """
    Cast sample and topo data to textstring for the online calculator. Input
    data are validated by utils.validate_*().

    sample : pandas Series or dict
        > sample_data = S.data.iloc[n]
        Missing data are replaced by default values, faulty data raise ValueError.
    topo : pandas Series or dict 
        'lat', 'long', 'elevation' (and 'shielding' if shielding='topo').
        > topo_data = topotable.iloc[n]
        or:
        > summary['elevation'] = summary['elevLo']
        > topo_data = summary
    shielding : str or numeric
        'sample', 'topo' or a numeric shielding factor. Defines whether
        shielding is taken from sample or raster data.
    
    Returns
    -------
    textline : str
        String corresponding to one location and one set of sample data (Be-10
        or Al-26).
        
    Raises
    ------
    TypeError : input data format is invalid
    ValueError : data cannot be validated
        
    """
    
    sample = deepcopy(sample)
    topo = deepcopy(topo)
    
    from riversand.utils import validate_topo, validate_nuclide, validate_sample
        
    "'standardization' is reset; faulty entries have no effect"
        
    if isinstance(sample, pd.DataFrame):
        if len(sample)!=1:
            raise TypeError("get_textline() argument 'sample' must be length 1")
        sample = sample.iloc[0] # recast to pd.Series
            
        if isinstance(sample, pd.Series):
            sample = sample.to_dict()
        
        if not isinstance(sample, dict):
            raise TypeError("get_textline() argument 'sample' must be dict or pandas Series")
    
    if isinstance(topo, pd.DataFrame):
        if len(topo)!=1:
            raise TypeError("get_textline() argument 'topo' must be length 1")
        topo = topo.iloc[0] # recast to pd.Series
            
        if isinstance(topo, pd.Series):
            topo = topo.to_dict()
        
        if not isinstance(topo, dict):
            raise TypeError("get_textline() argument 'topo' must be dict or pandas Series")

   
    # 'standardization' is set to 07KNSTD (Be-10) or KNSTD (Al-26)        
    try:
        #TODO: maybe issue a warning if incorrect standardization is defined
        sample.pop('standardization') #set to default by validate_nuclide())
    except:
        pass
    
    if shielding=='sample':
        try:
            sf = float(sample['shielding'])
        except:
            sf = np.nan
        finally:
            if np.isnan(sf):
                raise ValueError("Invalid shielding='sample': no shielding "+
                                 "defined in sample data")
    if shielding=='topo':
        try:
            sf = float(topo['shielding'])
        except:
            sf = np.nan
        finally:
            if np.isnan(sf):
                raise ValueError("Invalid shielding='topo': no shielding "+
                                 "defined in topo data")
            
    if isinstance(shielding, Number):
        sample['shielding'] = shielding
        shielding = 'sample'
    
    if shielding=='sample':
        try:
            out1 = validate_topo(topo)
            out2 = validate_nuclide(sample)
            out3 = validate_sample(sample, shielding=True)
        except ValueError as e:
            raise e
            
    elif shielding=='topo':
        try:
            out1 = validate_topo(topo, shielding=True)
            out2 = validate_nuclide(sample)
            out3 = validate_sample(sample)
        except ValueError as e:
            raise e
    else:
        raise ValueError("Invalid shielding: must be 'topo', 'sample' or numeric")
    
    out = {**out1, **out2, **out3}
    textline = "{} {:.5f} {:.5f} {:.3f} {} {} {} {:.5f} {:.5f} {} ; {} {} {} {} {} {} ;".format(
        out['name'], out['lat'], out['long'], out['elevation'], out['press_flag'],
        out['thickness'], out['density'], out['shielding'], out['erate'], out['year'],
        out['name'], out['nuclide'], out['mineral'], out['N'], out['delN'], out['standardization'])
    return textline    



def get_erates_from_server(
        textline:str,
        url:str=None
        ) -> (pd.DataFrame, str, dict):
    """
    Get erosion rates in g/cm2/yr from online calculator.
    
    textline : str
        textline returned by get_textline()
        > textline = get_textline(sample, topo, shielding='sample')
        or
        > out_summary['Elevation'] = out_summary['elevLo']
        > textline = get_textline(sample, out_summary, shielding='topo')
        May also be any valid textblock of multiple samples / nuclide / erates.
        
    Returns
    -------
    dfE : pd.DataFrame
        erosion rate results; one row per sample and nuclide.
        keys are: 'name', 'nuclide', 'density',
                  'E_St', 'delE_int_St', 'delE_ext_St',
                  'E_Lm', 'delE_int_Lm', 'delE_ext_Lm',
                  'E_LSDn', 'delE_int_LSDn', 'delE_ext_LSDn' 
    diagnostics : str
        "No response from {}".format(url)
        "No diagnostics" as returned by the server
    version : dict
        returned by the server.
        
    """
    
    dfE = pd.DataFrame(columns=['name', 'nuclide', 'density',
              'E_St', 'delE_int_St', 'delE_ext_St',
              'E_Lm', 'delE_int_Lm', 'delE_ext_Lm',
              'E_LSDn', 'delE_int_LSDn', 'delE_ext_LSDn'])
    version = {}
    
    from riversand import xml
    
    from riversand import params
    if url is None: url = params.url
        
    if not isinstance(textline, str):
        raise TypeError("get_erates_from_server() 'textline' must be string")
    if len(textline)==0:
        raise TypeError("'textline' empty string")
        
    # request erosion rate from server
    request_vars = {"mlmfile" : "erosion_input_v3",
                    "reportType" : "XML",
                    "resultType" : "long",
                    "summary" : "no",
                    "plotFlag" : "no",
                    "text_block" : textline }
    form_data = urllib.parse.urlencode(request_vars).encode('ascii')
    try:        
        result = urllib.request.urlopen(url, form_data)
    except urllib.error.URLError:
        diagnostics = "No response from {}".format(url)
    else:    
        result_xml = result.read()
        dfE, diagnostics, version = xml.read_erates(result_xml)
    
    return dfE, diagnostics, version
    

def get_E(textline:str, url:str=None) -> dict:
    """
    Get erosion rates in cm/yr from online calculator for a 
    *single* sample and nuclide.

    Returns
    -------
    ret_dict : dict
        keys 'St', 'Lm', 'LSDn'; erosion rate in cm/yr
    
    Raises RuntimeError if diagnostics is anything but "No diagnostics". 
       
    """

    from riversand import params
    if url is None: url = params.url
    
    dfE = pd.DataFrame()
    ret_dict = {'St' : None, 'Lm' : None, 'LSDn' : None}
        
    # erates in g/cm2/yr
    dfE, diagnostics, version = get_erates_from_server(
        textline, url=url)
    
    if "saturated" in diagnostics:
        raise RuntimeError("get_E() : sample appears to be saturated")
    elif "Empty string returned from the server" in diagnostics:
        raise RuntimeError("get_E() : nuclide concentration too low for calculation")
    elif "No response from" in diagnostics:
        raise RuntimeError("get_E() : {}".format(diagnostics))
        
        
    if len(dfE)!=1: 
        raise RuntimeError("get_E() : invalid results, probably from invalid "+
                           "argument 'textline'; must be a single sample and nuclide")
        
    try:
        ret_dict = {'St' : float(dfE.loc[0,'E_St']) / float(dfE.loc[0,'density']),
                    'Lm' : float(dfE.loc[0,'E_Lm']) / float(dfE.loc[0,'density']),
                    'LSDn' : float(dfE.loc[0,'E_LSDn']) / float(dfE.loc[0,'density']),
                    #'diagnostics' : diagnostics
                    }
    except:
        diagnostics = "Cannot compute E [g/cm2] from server output"
    
    if diagnostics!="No diagnostics":
        raise RuntimeError("get_E() : {}".format(diagnostics))
    return ret_dict

    
def get_NofE_from_server(
        textline:str,
        url:str=None
        ) -> (pd.DataFrame, str, dict):
    """
    Get predicted nuclide concentrations in at/g from online calculator.
    
    textline : str
        any valid textblock of multiple samples / nuclide / erates.
        
    Returns
    -------
    dfN : pd.DataFrame
        nuclide concentration results; one row per sample and nuclide.
        keys are: 'name', 'nuclide', 'E_cmyr',
                  'NpredSt', 'NpredLm', 'NpredLSDn' 
    diagnostics : str
        "No response from {}".format(url)
        "No diagnostics" as returned by the server
    version : dict
        returned by the server.
        
    """
    
    dfN = pd.DataFrame(columns = ['name', 'nuclide', 'E_cmyr',
                                  'NpredSt', 'NpredLm', 'NpredLSDn'])
    version = {}
    
    from riversand import xml
    
    from riversand import params
    if url is None: url = params.url
        
    if not isinstance(textline, str):
        raise TypeError("get_NofE_from_server() 'textline' must be string")
    if len(textline)==0:
        raise TypeError("'textline' is empty string")
        
    # request erosion rate from server
    request_vars = {"mlmfile" : "NofE_input_v3", "text_block" : textline}
    form_data = urllib.parse.urlencode(request_vars).encode('ascii')
    
    try:        
        result = urllib.request.urlopen(url, form_data)
    except urllib.error.URLError:
        diagnostics = "No response from {}".format(url)
    else:    
        result_xml = result.read() # extract result
        dfN, diagnostics, version = xml.read_NofE(result_xml)
    
    return dfN, diagnostics, version
    
    
    
def get_NofE_FullTable(
        sample_data:pd.Series,
        topostats: pd.DataFrame,
        shielding:str,
        erates:np.ndarray,
        url:str=None
    ) -> (pd.DataFrame, str, dict):
    """
    Get predicted nuclide concentrations for given topographic statistics
    for a single sample and a suite of erates.
    - lat, long, elevation from topostats
    - N, delN, press_flag, density, etc. from sample_data
    Sample names are auto-generated (one for each topostats entry).
    Sample thickness is 0.

    Parameters
    ----------
    sample_data : pd.Series or dict; single sample.
    topostats : pd.DataFrame; elevation-binned topo data.
    shielding : 'topo', 'sample' or numeric.
    erates : iterable or numeric.
    
    Returns
    ----------
    dfA : pd.DataFrame
        full table of nuclide predictions for each elevation bin and each erate
        keys are: 'name', 'nuclide', 'E_cmyr',
                  'NpredSt', 'NpredLm', 'NpredLSDn', # at/g as returned by the server
                  'St_wt', 'Lm_wt', 'LSDn_wt'        # at/g weighted by topostats['wt']
    diagnostics : str
    version : dict
        returned from the server.
        
    """
    
    from riversand import params
    if url is None: url = params.url
    
    if isinstance(erates, Number):
        erates = [erates]
    
    if len(erates)==0:
        raise TypeError("get_NofE_FullTable() argument 'erates' must have at least one value")
        
    if isinstance(sample_data, pd.DataFrame) and len(sample_data)!=1:
        raise TypeError("get_NofE_FullTable() argument 'sample_data' must "+
                         "be dict, pandas Series or single-row pandas Dataframe")
        
    # if topostats is dict or pd.Series cast to DataFrame
    # set 'wt'=1.0 if not defined
    # Note the .map at the end of this function! topostats must be pd.DataFrame with column 'wt'
    if isinstance(topostats, dict):
        temp = pd.DataFrame(topostats, index=[0])
        if 'wt' not in topostats.keys():
            temp['wt'] = 1.
        topostats = temp
        
    if isinstance(topostats, pd.Series): # e.g. if a single row of topostats is passed 
        temp = pd.DataFrame(data=[topostats])
        if 'wt' not in temp.keys():
            temp['wt'] = 1.
        topostats = temp
        
    # create text_block from all erates and all topotable entries
    tempS = deepcopy(sample_data)
    newline = '.'
    text_block = ''
    for e, erate in enumerate(erates):
        for i, tempT in topostats.iterrows():
            name_str = ("B_{}_{}".format(i,e))
            tempS['name'] = name_str
            tempS['erate'] = erate
            try:
                newline = get_textline(tempS,
                                       tempT, #topostats.iloc[i],
                                       shielding=shielding)
            except ValueError as e:
                raise ValueError("Cannot get valid input for the online calculator:\n   "+
                                 str(e))
            text_block = text_block + newline + ' '
             
    dfA, diagnostics, version = get_NofE_from_server(text_block, url=url)
    
    # apply weighting factor 'wt' to calculate 'St', 'Lm', 'LSDn'
    if len(dfA)>0:
        # generate a column 'topo' in order to map the weights from dfT.wt to the samples
        dfA[['n', 'topo', 'erate']] = dfA['name'].str.split("_", expand=True)
        dfA.drop(columns=['n', 'erate'], inplace=True)
        
        # apply weighting factor to the exposure ages and get weighted sum of exposure ages
        for x,Nx in zip(params.scalingmethods, params.XML_scalingmethods):
            # x='St', 'Nx'='NpredSt' etc:
            dfA[x] = dfA['topo'].astype(int).map(topostats['wt']) * dfA[Nx].astype(float)
        
        dfA['E_cmyr'] = dfA['E_cmyr'].apply(pd.to_numeric)
        dfA.drop(columns=['topo'], inplace=True)
        
    return dfA, diagnostics, version


def get_NofE(
        sample_data:pd.Series,
        topostats:pd.DataFrame,
        shielding:str,
        erates:np.ndarray,
        url:str=None
    ) -> pd.DataFrame:
    """
    Wrapper for get_NofE_FullTable:
    Predicted nuclide concentrations for a given catchment topography
    for a single sample and a suite of erosion rates.
    
    Parameters
    ----------
    sample_data : pd.Series or dict; single sample.
    topostats : pd.DataFrame; elevation-binned topo data.
    shielding : 'topo', 'sample' or numeric.
    erates : iterable or scalar; (suite of) erosion rates in cm/yr.
    
    Returns
    ----------
    NofE : pd.DataFrame
        Predicted nuclide concentrations for each erosion rate in erates
        and each scaling method
        keys are : 'St', 'Lm', 'LSDn'
        index name is 'E_cmyr'.
    
    Raises RuntimeError if diagnostics is anything but "No diagnostics".
        
    """

    from riversand.utils import validate_topo, validate_nuclide, validate_sample

    if not isinstance(topostats, pd.DataFrame):
        raise TypeError("get_NofE() argument 'topostats' must be pandas DataFrame")

    if isinstance(sample_data, pd.DataFrame):
        if len(sample_data)==1: # recast single-line DataFrame to Series
            sample_data = sample_data.iloc[0]
                    
    if not isinstance(sample_data, (pd.Series, dict)):
        raise TypeError("get_NofE() argument 'sample_data' must be dict, "+
                        "pandas Series or single-row pandas Dataframe")
    
    # validate sample_data and a line of topodata to avoid Exceptions raise by
    # get_NofE_FullTable() (ValueError "Cannot get valid input...")
    try:
        _ = validate_topo(topostats.iloc[0])
    except:
        raise ValueError("Invalid values in 'topostats' or empty table")
    try:
        _ = validate_sample(sample_data)
    except:
        raise ValueError("Invalid values in 'sample_data'")
    try:
        _ = validate_nuclide(sample_data)
    except:
        raise ValueError("Invalid values in 'sample_data'")
        
    if (shielding=='topo') and ('shielding' not in topostats.keys()):
        raise ValueError("shielding='topo' but no shielding factor defined in 'topostats'")
            
    if (shielding=='sample'):
        try:
            _ = validate_sample(sample_data, shielding=True)
        except:
            raise ValueError("shielding='sample' but no shielding factor defined in 'sample_data'")
        
    dfA, diagnostics, version = get_NofE_FullTable(
        sample_data, topostats, shielding, erates, url=url)
    
    # saturated samples or low nuclide concentration are never an issue with NofE:
    if "No response from" in diagnostics:
        raise RuntimeError("get_NofE() : {}".format(diagnostics))
                        
    if len(dfA)==0: # not sure why this would ever happen
        raise RuntimeError("unexpected error in get_NofE() : empty data table")
    
    # 'sample_name' and 'erosion_rate_cm_yr' are tags returned by the server
    if dfA.groupby('E_cmyr').count().loc[:,'name'].nunique()==1:
        # if grouping by erosion rate is done correctly each group should have the exact same number of entries
        NofE = dfA.groupby('E_cmyr').sum(numeric_only=True)
        
    else:
        # rounding errors in the erosion rates sent to the server may lead to problems
        temp = dfA['name'].str.split('_', expand=True)
        dfA['temp'] = temp[2]
        temp = dfA.groupby('temp').sum()
        tempE = dfA.groupby('temp').mean()
        temp['E_cmyr'] = tempE['E_cmyr']
        NofE = temp.set_index('E_cmyr')
    return NofE


def guess_erates(*args, **kwargs) -> np.ndarray:
    """
    Generate a suite of initial erosion rates in cm/yr.

    One positional argument : erates = 0.5*E ... 2*E
    Two positional arguments : erates = E1 ... E2
    Positional arguments can be dict of erosion rates with keys 'St', 'Lm',
        'LSDn' if 'scaling' is specified as keyword argument.        
    'N' defines the number of output erosion rates; defaults to N=6

    """
    
    if 'N' in kwargs.keys(): N = kwargs['N']
    else: N = 6

    if 'scaling' in kwargs.keys(): scaling = kwargs['scaling']
    else: scaling = None
        
    from riversand import params
    
    if len(args)==2:
        E1 = args[0]
        E2 = args[1]
        
        if isinstance(E1, Number) and isinstance(E2, Number):
            pass
        elif isinstance(E1, dict) and isinstance(E2, dict):
            try:
                E1 = E1[scaling]
                E2 = E2[scaling]
            except:
                raise ValueError("guess_erates() invalid keyword argument "+
                                 "scaling '{}'"
                                 .format(scaling))
        else:
            raise ValueError("guess_erates() invalid arguments")
       
    elif len(args)==1:
        E = args[0]
        
        if isinstance(E, Number):
            pass
        elif isinstance(E, dict):
            try:
                E = E[scaling]
            except:
                raise ValueError("guess_erates() invalid keyword argument "+
                                 "scaling '{}'"
                                 .format(scaling))
        else:
            raise ValueError("guess_erates() invalid arguments")
        E1 = E/2
        E2 = 2*E
    else:
        raise ValueError("guess_erates() invalid arguments")
        
    minE = min(E1, E2)
    maxE = max(E1, E2)
    
    minE = np.log(max([minE, 0.001*maxE]))
    maxE = np.log(maxE)
    return np.exp(np.linspace(minE, maxE, int(N))) # cm/yr


def NofE_fitfunc(x, a, b, c):
    return  a/x**2 + b/x + c


def poly_E_results(
        sample_data:pd.Series,
        topostats:pd.DataFrame,
        shielding:str,
        erates:np.ndarray,
        scaling:str,
        url:str=None,
        strict=True
    ) -> (Number, tuple, pd.DataFrame, list):
    """
    Calculate catchmentwide erosion rate E and uncertainty delE for a single
    sample based on 'topostats' and a suite of initial erosion rates.
    
    sample_data : pd.Series or dict; single sample.
    topostats : pd.DataFrame; elevation-binned topo data.
    shielding : 'topo', 'sample' or numeric.
    erates : iterable or scalar; (suite of) erosion rates in cm/yr.
    scaling : str; scaling method 'St', 'Lm' or 'LSDn'.
    
    Returns
    ----------
    E : float; erosion rate in cm/yr.
    delE : tuple of floats, uncertainty (-delE, +delE).
    NofE : pd.Series with keys defined by 'scaling' and index 'E_cmyr'.
    error : list or error strings. In case of 'minE too high' or 'maxE too low'
        the list has only this error!
    
    """

    from riversand import params
    if url is None: url = params.url
    
    if scaling not in params.scalingmethods:
        raise ValueError("Invalid scaling '{}' (must be 'St', 'Lm' or 'LSDn')")
    # other requirements are validated by get_NofE() and raised as TypeError or ValueError
    
    if isinstance(erates, Number):
        erates = [erates]
        
    if strict:
        if len(erates)<4:
            raise RuntimeError("poly_E_results() : argument 'erates' should "+
                               "have at least 4 values; use argument "+
                               "'strict=False' to override")
    # return nan and error string in case of error
    y = None
    E_poly = np.nan
    delE = (np.nan, np.nan)
    
    error = []
    
    erates[:] = np.sort(erates)
    
    try:
        NofE = get_NofE(sample_data, topostats, shielding, erates, url=url)
    except TypeError as e: # input data format
        raise e
    except ValueError as e: # input data
        raise e
    except RuntimeError as e: # no server response
        raise e
        
    if not NofE.index.is_unique:
        # rounding errors in the erates sent to the server may result in non-unique indices;
        # drop duplicates
        NofE = NofE[~NofE.index.duplicated(keep='first')]
        error += ["Server cannot resolve erosion rates, dropping duplicates"]
        #raise Warning("Server cannot resolve the provided erosion rates, dropping duplicates")
    
    Nmeas = float(sample_data['N'])
    
    x = NofE.index
    y = NofE[scaling]
    
    if Nmeas >= y.iloc[0]:
        error = ['minE too high']
        #"ERROR: Erosion rate is outside of specified bracket (<{} cm/yr)".format(x[0]))
        
    elif Nmeas <= y.iloc[-1]:
        error = ['maxE too low']
        #"ERROR: Erosion rate is outside of specified bracket (>{} cm/yr)".format(x[-1]))
        
    else:
        popt, pcov = optimize.curve_fit(NofE_fitfunc, x, y)
        a,b,c = popt
        sol = optimize.root_scalar(NofE_fitfunc, args=(a, b, c-Nmeas), method='toms748', bracket=[x[0], x[-1]])
            
        if sol.converged:
            E_poly = sol.root            

        else:
            error += ["Root finding did not converge"]
            #"ERROR: No solution found"

        try:
            delE = uncert_on_E(E_poly, y, sample_data)
        except:
            error += ["Cannot compute uncertainty"]
        
        #RMSE = np.sqrt(np.sum((y - NofE_fitfunc(x, *popt))**2)/len(y))
        RMSE = get_RMSE(pd.Series(y, index=x))
        if (RMSE / Nmeas)>1e-3:
            error += ["NRMSE = {:.2e} suggests a poor fit of the polynomial".format(RMSE / Nmeas)]
            #raise Warning("NRMSE = {:.2e} suggests a poor fit of the polynomial!".format(RMSE / Nmeas))
            #print("NRMSE={:.2e}".format(RMSE/Nmeas))
                         
    return E_poly, delE, NofE[scaling], error



def uncert_on_E(E:float, NofE:pd.Series, sample_data:pd.Series) -> tuple:
    """ Uncertainty on E from uncertainty on N and curve fit to NofE. """
    
    if isinstance(sample_data, pd.DataFrame):
        if len(sample_data)==1:
            sample_data = deepcopy(sample_data.iloc[0])
        else:
            raise TypeError("uncert_on_E() argument 'sample_data' must be "+
                            "dict, pandas Series or single-row pandas Dataframe")
    
    N1 = sample_data['N']+sample_data['delN']
    N2 = sample_data['N']-sample_data['delN']
    
    x = NofE.index
    popt, pcov = optimize.curve_fit(NofE_fitfunc, x, NofE)
    a,b,c = popt

    # find values x1 and x2 that bracket the erosion rates E1, E2 (E+/-delE)
    x1, y1 = min(NofE.index), max(NofE)
    x2, y2 = max(NofE.index), min(NofE)
    
    while y1<N1: # stepwise 10% decrease of x1
        x1 = 0.9*x1
        y1 = NofE_fitfunc(x1, a, b, c)
    
    while y2>N2: # stepwise 10% increase of x2
        x2 = 1.1*x2
        y2 = NofE_fitfunc(x2, a, b, c)
    
    # find erosion rates corresponding to N+delN and N-delN
    sol = optimize.root_scalar(NofE_fitfunc, args=(a, b, c-N1), method='toms748', bracket=[x1, E])
    if sol.converged:
        E1 = sol.root
        delE1 = E-E1
        
    sol = optimize.root_scalar(NofE_fitfunc, args=(a, b, c-N2), method='toms748', bracket=[E, x2])
    if sol.converged:
        E2 = sol.root
        delE2 = E2-E

    return (delE1, delE2)



def get_RMSE(NofE):   
    """ Root mean squared error. """
    popt, pcov = optimize.curve_fit(NofE_fitfunc, NofE.index, NofE)
    return np.sqrt(np.sum((NofE - NofE_fitfunc(NofE.index, *popt))**2)/len(NofE))


def get_version(url:str=None) -> dict:
    
    from riversand import params
    if url is None: url = params.url
    textline = get_textline({'N':1e5, 'delN':1e3}, {'lat':0, 'long':0, 'elevation':0}, shielding=1)
    
    df, diagn, version = get_erates_from_server(textline, url=url)
    
    if "No response from" in diagn:
        raise RuntimeError("get_version() : {}".format(diagn))
        
    return version

def main():
    pass

if __name__ == "__main__":
    main()