#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

*******************************************************************************
xml.py  :  decoding of the xml server output 

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

Conversion of xml-strings returned by the server to
a pandas DataFrame, a diagnostics string and a version dictionary.

> request_vars = {"mlmfile" : "erosion_input_v3",
                  "reportType" : "XML",
                  "resultType" : "long",
                  "summary" : "no",
                  "plotFlag" : "no",
                  "text_block" : textline }
> form_data = urllib.parse.urlencode(request_vars).encode('ascii')
> result = urllib.request.urlopen(url, form_data)
> result_xml = result.read()

read_erates(result_xml)
read_NofE(result_xml)
read_ages(result_xml) # not implemented

"""

import pandas as pd
import xml.etree.ElementTree as et


def read_erates(result_xml:str) -> (pd.DataFrame, str, dict):
    """    
    Returns
    -------
    df : pandas DataFrame with one row per sample and nuclide.
        Units are density: g/cm3, erosion rate: g/cm2/yr
    diagnostics : str
        diagnostics string.
    version : dict
        version info returned from the server.
    """
    
    version = {}
    cols = ['sample_name', 'nuclide', 'sample_density',
        'E_gcm2_St',   'delE_Int_gcm2_St',   'delE_Ext_gcm2_St',
        'E_gcm2_Lm',   'delE_Int_gcm2_Lm',   'delE_Ext_gcm2_Lm',
        'E_gcm2_LSDn', 'delE_Int_gcm2_LSDn', 'delE_Ext_gcm2_LSDn']
    df = pd.DataFrame(columns=cols)
    
    out = df.rename(columns={
        'sample_name' : 'name', 'sample_density': 'density',
        'E_gcm2_St':   'E_St',   'delE_Int_gcm2_St':   'delE_int_St',   'delE_Ext_gcm2_St':   'delE_ext_St',
        'E_gcm2_Lm':   'E_Lm',   'delE_Int_gcm2_Lm':   'delE_int_Lm',   'delE_Ext_gcm2_Lm':   'delE_ext_Lm',
        'E_gcm2_LSDn': 'E_LSDn', 'delE_Int_gcm2_LSDn': 'delE_int_LSDn', 'delE_Ext_gcm2_LSDn': 'delE_ext_LSDn',
        })
    
    if len(result_xml)==0:
        ret_str = "Empty string returned from the server; this may happen for a very small nuclide concentration"
        return out, ret_str, {}
    
    try:
        root = et.fromstring(result_xml)
    except Exception as err:
        print(err)
        return out, err, {}

    tags = [i.tag for i in root]
    
    iD = [i for i, tag in enumerate(tags) if tag=="diagnostics"]
    iV = [i for i, tag in enumerate(tags) if tag=="version"]
    assert len(iD)==1 # exactly 1 Element <diagnostics>
    assert len(iV)<=1 # max 1 Element <version>; may be 0 for invalid data

    iE = [i for i, tag in enumerate(tags) if tag=="erosionRateResult"]
    
    # diagnostics:
    if len(iD)==1:
        dignosticsElement = root[iD[0]]
        diagnostics = dignosticsElement.text
        
        if diagnostics=="No diagnostics":
            pass
        elif 'saturated' in diagnostics:
            pass
        else:
            # This is a bug in the server output: if there are any issues with
            # the data the server returns some diagnostics but no version info,
            # and the root tag is "exposureAgeResult" even if erosion rate was
            # requested. The data is all jumbled zeroes
            return out, diagnostics, {}
        
    # version dict:
    if len(iV)==1:
        versionElement = root[iV[0]]
        for item in versionElement:
            version[item.tag] = item.text
   
    i = 0
    for erateResult in [root[i] for i in iE]:
        assert erateResult.tag=='erosionRateResult'
        sample_name = None
        sample_density = None
        
        for item in erateResult:
            assert item.tag in {'sample_name', 'sample_density', 'nuclide_result'}
            
            if item.tag=='sample_name':
                sample_name = item.text
            
            if item.tag=='sample_density':
                sample_density = item.text
            
            if item.tag=='nuclide_result':
                this_nuclide = {'sample_name': sample_name,
                                'sample_density': sample_density
                               }
                delInt = True
                for E in item:                    
                    if 'delE' in E.tag: # sort out internal and external errors
                        
                        if delInt==True:
                            Etag = E.tag[:4]+"_Int"+E.tag[4:]
                            delInt = False
                        else:
                            Etag = E.tag[:4]+"_Ext"+E.tag[4:]
                            delInt = True
                        this_nuclide[Etag] = E.text # <nuclide_result> to dict
                    else:
                        this_nuclide[E.tag] = E.text # <nuclide_result> to dict
    
                df.loc[i, :] = this_nuclide
                i += 1
                    
    df.rename(columns={
        'sample_name' : 'name', 'sample_density': 'density',
        'E_gcm2_St':   'E_St',   'delE_Int_gcm2_St':   'delE_int_St',   'delE_Ext_gcm2_St':   'delE_ext_St',
        'E_gcm2_Lm':   'E_Lm',   'delE_Int_gcm2_Lm':   'delE_int_Lm',   'delE_Ext_gcm2_Lm':   'delE_ext_Lm',
        'E_gcm2_LSDn': 'E_LSDn', 'delE_Int_gcm2_LSDn': 'delE_int_LSDn', 'delE_Ext_gcm2_LSDn': 'delE_ext_LSDn',
        },
        inplace=True)
    
    return df, diagnostics, version


def read_NofE(result_xml:str) -> (pd.DataFrame, str, dict):
    """    
    Returns
    -------
    df : pandas DataFrame with one row per sample and nuclide.
        Units are erosion rate: g/cm2/yr, nuclide concentrations atoms/gram
    diagnostics : str
        diagnostics string.
    version : dict
        version info returned from the server.
    """
    
    version = {}
    cols = ['sample_name', 'nuclide', 'erosion_rate_cm_yr',
            'NpredSt', 'NpredLm', 'NpredLSDn']
    df = pd.DataFrame(columns=cols)
    
    out = df.rename(columns={
        'sample_name' : 'name', 'erosion_rate_cm_yr': 'E_cmyr',
        })
    
    if len(result_xml)==0:
        ret_str = "Empty string returned from the server; this may happen for a very small nuclide concentration"
        return out, ret_str, {}
    
    try:
        root = et.fromstring(result_xml)
    except Exception as err:
        print(err)
        return out, err, {}

    tags = [i.tag for i in root]
    
    iD = [i for i, tag in enumerate(tags) if tag=="diagnostics"]
    iV = [i for i, tag in enumerate(tags) if tag=="version"]
    assert len(iD)==1 # exactly 1 Element <diagnostics>
    assert len(iV)<=1 # max 1 Element <version>; may be 0 for invalid data

    iE = [i for i, tag in enumerate(tags) if tag=="NofEResult"]
    
    # diagnostics:
    if len(iD)==1:
        dignosticsElement = root[iD[0]]
        diagnostics = dignosticsElement.text
        
        if diagnostics=="No diagnostics":
            pass
        elif 'saturated' in diagnostics:
            pass
        else:
            # This is a bug in the server output: if there are any issues with
            # the data the server returns some diagnostics but no version info,
            # and the root tag is "exposureAgeResult" even if erosion rate was
            # requested. The data is all jumbled zeroes
            return out, diagnostics, {}
        
    # version dict:
    if len(iV)==1:
        versionElement = root[iV[0]]
        for item in versionElement:
            version[item.tag] = item.text
   
    i = 0
    for NofEResult in [root[i] for i in iE]:
        assert NofEResult.tag=='NofEResult'
        sample_name = None
        erosion_rate_cm_yr = None
        
        for item in NofEResult:
            assert item.tag in {'sample_name', 'erosion_rate_cm_yr', 'nuclide_result'}
            
            if item.tag=='sample_name':
                sample_name = item.text
            
            if item.tag=='erosion_rate_cm_yr':
                erosion_rate_cm_yr = item.text
            
            if item.tag=='nuclide_result':
                this_nuclide = {'sample_name': sample_name,
                                'erosion_rate_cm_yr': erosion_rate_cm_yr
                               }
                for E in item:
                    this_nuclide[E.tag] = E.text # <nuclide_result> to dict
                        
                df.loc[i, :] = this_nuclide
                i += 1
                    
    df.rename(columns={
        'sample_name' : 'name', 'erosion_rate_cm_yr': 'E_cmyr',
        },
        inplace=True)
    
    return df, diagnostics, version