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

ages_from_xml(result_xml)
erates_from_xml(result_xml)

obsolete read_ages(result_xml)
obsolete read_erates(result_xml)
NofE_from_xml(result_xml)

"""

import pandas as pd
import xml.etree.ElementTree as et


def ages_from_xml(result_xml:str, norm=True) -> (pd.DataFrame, str, dict):
    """
    Translate <calcs_v3_age_data> xml string to pandas DataFrame.
    norm=False exludes the normalized nuclide concentrations.

    Returns
    -------
    df : pandas DataFrame with one row per sample and nuclide.
        Units are years.
    diagnostics : str
        diagnostics string.
    version : dict
        version info returned from the server.
    """
    
    """
    xml (v 2024-05-10)
    root : <calcs_v3_age_data>
    children :
        <exposureAgeResult>   # one for each sample site
        <diagnostics>
        <version>
        
    <exposureAgeResult> :
        item <sample_name>   # one for each sample site
        + following 18 items; sequence repeats if more than one analysis per sample site
        item <t10quartz_St>
        item <delt10quartz_int_St>
        item <delt10quartz_ext_St>
        item <Nnorm10quartz_St>
        item <delNnorm10quartz_int_St>
        item <delNnorm10quartz_ext_St>
        item <t10quartz_Lm>
        item <delt10quartz_int_Lm>
        item <delt10quartz_ext_Lm>
        item <Nnorm10quartz_Lm>
        item <delNnorm10quartz_int_Lm>
        item <delNnorm10quartz_ext_Lm>
        item <t10quartz_LSDn>
        item <delt10quartz_int_LSDn>
        item <delt10quartz_ext_LSDn>
        item <Nnorm10quartz_LSDn>
        item <delNnorm10quartz_int_LSDn>
        item <delNnorm10quartz_ext_LSDn>
    """
    
    version = {}
    # tags returned by the server
    Be10cols = ['t10quartz_St', 'delt10quartz_int_St', 'delt10quartz_ext_St',
                'Nnorm10quartz_St', 'delNnorm10quartz_int_St', 'delNnorm10quartz_ext_St',
                't10quartz_Lm', 'delt10quartz_int_Lm', 'delt10quartz_ext_Lm',
                'Nnorm10quartz_Lm', 'delNnorm10quartz_int_Lm', 'delNnorm10quartz_ext_Lm',
                't10quartz_LSDn', 'delt10quartz_int_LSDn', 'delt10quartz_ext_LSDn',
                'Nnorm10quartz_LSDn', 'delNnorm10quartz_int_LSDn', 'delNnorm10quartz_ext_LSDn']
    
    # empty dataframe formatted to match the output in case of an exception
    out = pd.DataFrame(columns=['sample_name', 'nuclide']+[c.replace('10quartz', '') for c in Be10cols])
    out = out.rename(columns={'sample_name': 'sample'})

    if len(result_xml)==0:
        ret_str = "Empty string returned from the server. This may happen if at least one sample has an extremely low nuclide concentration."
        return out, ret_str, version

    try:
        root = et.fromstring(result_xml)
    except Exception as err:
        print(err)
        return out, err, version

    tags = [i.tag for i in root] # ['exposureAgeResult', 'exposureAgeResult', ..., 'diagnostics', 'version']
    # note that a tag 'exposureAgeResult' may have more than one result

    iD = [i for i, tag in enumerate(tags) if tag=="diagnostics"]
    iV = [i for i, tag in enumerate(tags) if tag=="version"]
    assert len(iD)==1 # exactly 1 element <diagnostics>
    assert len(iV)<=1 # max 1 element <version>; for erosion rates this may be 0 for invalid data

    # version dict:
    if len(iV)==1:
        versionElement = root[iV[0]]
        for item in versionElement:
            version[item.tag] = item.text
            
    # diagnostics:
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
        pass
        return out, diagnostics, version
        
    ## START xml conversion
    nuclide_dict = {'10quartz': 'Be-10', '26quartz': 'Al-26'}
    df = pd.DataFrame(columns=['sample_name', 'nuclide']) # tag 'sample_name' returned by the server

    i = -1 # index in result DataFrame
    for element in root: # iterate over 'exposureAgeResult's; each may have more than one analysis
        if element.tag=='exposureAgeResult': # skip the 'dignostics' and 'version' elements
            for item in element:
                #assert item.tag in {'sample_name', 't10quartz_St', 'delt10quartz_int_St', etc}
                
                if item.tag=='sample_name':
                    sample_name = item.text # name for one or several analyses
                    
                elif item.tag in {'t10quartz_St', 't26quartz_St'}:
                    i += 1 # new analysis i.e. next line in results
                    nuclide = item.tag[1:9] # nuclide string '10quartz', '26quartz'
                    df.loc[i, 'sample_name'] = sample_name
                    df.loc[i, 'nuclide'] = nuclide_dict[nuclide] # translate to 'Be-10', 'Al-26'
                    col = item.tag.replace(nuclide, '') # remove the 'xxquartz' substring
                    df.loc[i, col] = item.text
                else:
                    col = item.tag.replace(nuclide, '') # remove the 'xxquartz' substring
                    df.loc[i, col] = item.text
    
    df = df.rename(columns={'sample_name': 'sample'})
    
    #df.nuclide = df.nuclide.map(pd.Series(repl)) # replace '10quartz' by 'Be-10'
    if norm:
        return df, diagnostics, version
    else: # short table excluding the normalized values
        return df.iloc[:,[0,1,2,3,4,8,9,10,14,15,16]], diagnostics, version
    
    
def erates_from_xml(result_xml:str) -> (pd.DataFrame, str, dict):
    """
    Tanslate <calcs_v3_erosion_data> xml string to pandas DataFrame.
    
    Returns
    -------
    df : pandas DataFrame with one row per sample and nuclide.
        Units are density: g/cm3, erosion rate: g/cm2/yr.
        Note that the m/Myr output of the browser version is not part of the
        xml string and can be calculated from the density.
    diagnostics : str
        diagnostics string.
    version : dict
        version info returned from the server.
    """
    
    """
    xml (v 2024-05-10)
    root: <calcs_v3_erosion_data>
    children :
        <erosionRateResult>   # one for each sample site
        <diagnostics>
        <version>
    
    <erosionRateResult> :
        item <sample_name>   # one for each sample site
        item <sample_density>   # one for each sample site
        item <nuclide_result>   # one or several if more than one analysis per sample site
            item <nuclide>
            item <E_gcm2_St>
            item <delE_gcm2_St> # int error
            item <delE_gcm2_St> # ext error
            item <E_gcm2_Lm>
            item <delE_gcm2_Lm>
            item <delE_gcm2_Lm>
            item <E_gcm2_LSDn>
            item <delE_gcm2_LSDn>
            item <delE_gcm2_LSDn>
    """
    
    version = {}
    # tags returned by the server are intentical for int and ext errors!
    cols = ['sample_name', 'nuclide', 'sample_density', 
            'E_gcm2_St',   'delE_int_gcm2_St',   'delE_ext_gcm2_St',
            'E_gcm2_Lm',   'delE_int_gcm2_Lm',   'delE_ext_gcm2_Lm',
            'E_gcm2_LSDn', 'delE_int_gcm2_LSDn', 'delE_ext_gcm2_LSDn']
    
    # empty dataframe formatted to match the output in case of an exception
    out = pd.DataFrame(columns=cols)
    out = out.rename(columns={
        'sample_name': 'sample', 'sample_density': 'density'})
    
    if len(result_xml)==0:
        ret_str = "Empty string returned from the server. This may happen if at least one sample has an extremely low nuclide concentration."
        return out, ret_str, version
    
    try:
        root = et.fromstring(result_xml)
    except Exception as err:
        print(err)
        return out, err, version

    tags = [i.tag for i in root] # ['erosionRateResult', 'erosionRateResult', ..., 'diagnostics', 'version']
    # note that 'erosionRateResult' may have more than one result
    
    iD = [i for i, tag in enumerate(tags) if tag=="diagnostics"]
    iV = [i for i, tag in enumerate(tags) if tag=="version"]
    assert len(iD)==1 # exactly 1 element <diagnostics>
    assert len(iV)<=1 # max 1 element <version>; may be 0 for invalid data
    
    # version dict:
    if len(iV)==1:
        versionElement = root[iV[0]]
        for item in versionElement:
            version[item.tag] = item.text
            
    # diagnostics:
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
        return out, diagnostics, version
        
    ## START xml conversion
    nuclide_dict = {'N10quartz': 'Be-10', 'N26quartz': 'Al-26'}    
    df = pd.DataFrame(columns=cols) # these are the tags actually returned from the server
    
    i = -1
    for element in root: # iterate over 'erosionRateResult's; each may have more than one analysis
        if element.tag=='erosionRateResult': # skip the 'dignostics' and 'version' elements
            sample_name = None
            sample_density = None
            
            for item in element:
                assert item.tag in {'sample_name', 'sample_density', 'nuclide_result'}
                
                if item.tag=='sample_name':
                    sample_name = item.text # name for one or several analyses
                
                if item.tag=='sample_density':
                    sample_density = item.text
                
                if item.tag=='nuclide_result':
                    i += 1 # new analysis i.e. next line in results
                    
                    df.loc[i, 'sample_name'] = sample_name
                    df.loc[i, 'sample_density'] = sample_density
                
                    # sort out internal and external errors with identical item tags
                    delInt = True # errors are reported first int then ext
                    for itm in item:
                        
                        if 'delE' in itm.tag:
                            if delInt==True:
                                tag = itm.tag[:4]+"_int"+itm.tag[4:]
                                delInt = False
                            else:
                                tag = itm.tag[:4]+"_ext"+itm.tag[4:]
                                delInt = True
                                
                            df.loc[i, tag] = itm.text # <delE_gcm2_XXX>
                        
                        elif itm.tag=='nuclide':        
                            df.loc[i, 'nuclide'] = nuclide_dict[itm.text] # translate to 'Be-10', 'Al-26'
                        
                        else:
                            df.loc[i, itm.tag] = itm.text # <E_gcm2_XX>
                            
                    
    df = df.rename(columns={
        'sample_name': 'sample', 'sample_density': 'density'})
    
    return df, diagnostics, version

    

def NofE_from_xml(result_xml:str) -> (pd.DataFrame, str, dict):
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