#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 11:00:51 2023

@author: Konstanze Stübner, kstueb@gmail.com

"""

import pytest
from pytest import mark

import sys, os
import pandas as pd
import numpy as np


from copy import deepcopy

import riversand
from riversand.calc import get_NofE_FullTable


class Test_get_textline():
    
    # - raise exception for incomplete / invalid topo data     
    # - incorrect standardizaion is silently overwritten 

    def test_get_textline__dict_dict_valid(self, rv_single):
        """ default values as dict """
        
        topo = {'lat':33, 'long':72, 'elevation':0}
        sample = {'N':1e5, 'delN':1e2}
        outstr = riversand.get_textline(sample, topo, shielding=1)
        assert outstr==("Test 33.00000 72.00000 0.000 std 0 2.65 1.00000 0.00000 2010 ; "+
                        "Test Be-10 quartz 100000 100 07KNSTD ;")
                
    def test_get_textline__dict_dict1(self, rv_single):
        """ shielding is in conflict with topo data """
        
        topo = {'lat':33, 'long':72, 'elevation':0}
        sample = {'N':1e5, 'delN':1e2}
        outstr = riversand.get_textline(sample, topo, shielding=1)
        
        with pytest.raises(ValueError) as e:
            riversand.get_textline(sample, topo, shielding='topo')
        assert e.value.args[0]=="Invalid shielding 'topo' (no shielding defined in 'topo')"
        
    def test_get_textline__dict_dict2(self, rv_single):
        """ shielding is in conflict with sample data """
        
        topo = {'lat':33, 'long':72, 'elevation':0}
        sample = {'N':1e5, 'delN':1e2}
        outstr = riversand.get_textline(sample, topo, shielding=1)
        
        with pytest.raises(ValueError) as e:
            riversand.get_textline(sample, topo, shielding='sample')
        assert e.value.args[0]=="Invalid shielding 'sample' (no shielding defined in 'sample')"
    
    def test_get_textline__dict_dict3(self, rv_single):
        """ input does not get modified """
        
        topo = {'lat':33, 'long':72, 'elevation':0}
        sample = {'N':1e5, 'delN':1e2}
        outstr = riversand.get_textline(sample, topo, shielding=1)
        assert topo=={'lat':33, 'long':72, 'elevation':0}
        assert sample=={'N':1e5, 'delN':1e2}

    def test_get_textline__dict_dict4(self, rv_single):
        """ invalid standardization is reset without warning """
        
        topo = {'lat':33, 'long':72, 'elevation':0}
        sample = {'N':1e5, 'delN':1e2, 'standardization': 'NIST'}
        outstr = riversand.get_textline(sample, topo, shielding=1)
        assert outstr==("Test 33.00000 72.00000 0.000 std 0 2.65 1.00000 0.00000 2010 ; "+
                        "Test Be-10 quartz 100000 100 07KNSTD ;")
        
    def test_get_textline__dict_dict5(self, rv_single):
        """ invalid sample data """
        
        topo = {'lat':33, 'long':72, 'elevation':0}
        sample = {'N':1e5, 'delN':1e2, 'density': 0}
        with pytest.raises(ValueError) as e:
            riversand.get_textline(sample, topo, shielding=1)
        assert e.value.args[0]=="Density is 0 or less"
    
    def test_get_textline__dict_dict6(self, rv_single):
        """ invalid topo data """
        
        topo = {'lat':-100, 'long':72, 'elevation':0}
        sample = {'N':1e5, 'delN':1e2, 'density': 1}
        with pytest.raises(ValueError) as e:
            riversand.get_textline(sample, topo, shielding=1)
        assert e.value.args[0]=="Latitude is out of bounds (-90..+90)"
    
    def test_get_textline__dict_dict7(self, rv_single):
        """ missing topo data """
        
        topo = {'long':72, 'elevation':0}
        sample = {'N':1e5, 'delN':1e2, 'density': 0}
        with pytest.raises(ValueError) as e:
            riversand.get_textline(sample, topo, shielding=1)
        assert e.value.args[0]=="Latitude is not defined"
    
    
    def test_get_textline__sample_as_Series(self, rv_single):
        """
        topo as dict, sample as Series
        - raise exception for invalid sample data
        - raise exception if shielding is in conflict with the topo/sample data 
        - incorrect standardizaion is silently overwritten 
        """
        # default values
        topo = {'lat':33, 'long':72, 'elevation':0}
        sample = rv_single.samples.iloc[0]
        assert 'standardization' not in sample.keys()
        assert 'thickness' not in sample.keys()
        
        assert isinstance(sample, pd.Series)
        outstr = riversand.get_textline(sample, topo, shielding=1)
        assert outstr==("DB02 33.00000 72.00000 0.000 std 0 2.7 1.00000 0.00000 2000 ; "+
                        "DB02 Be-10 quartz 23500 1400 07KNSTD ;")
        
        # input values do not get modified  
        assert 'standardization' not in sample.keys()
        assert 'thickness' not in sample.keys()
        
        with pytest.raises(ValueError) as e:
            riversand.get_textline(sample, topo, shielding='sample')
        assert e.value.args[0]=="Invalid shielding 'sample' (no shielding defined in 'sample')"
    
        # invalid standardization is reset without warning
        topo = {'lat':33, 'long':72, 'elevation':0}
        sample = deepcopy(rv_single.samples.iloc[0])
        sample['standardization']='NIST'
        outstr = riversand.get_textline(sample, topo, shielding=1)
        assert outstr==("DB02 33.00000 72.00000 0.000 std 0 2.7 1.00000 0.00000 2000 ; "+
                        "DB02 Be-10 quartz 23500 1400 07KNSTD ;")
        
        #invalid data:
        topo = {'lat':33, 'long':72, 'elevation':0}
        sample['density']=0
        with pytest.raises(ValueError) as e:
            riversand.get_textline(sample, topo, shielding=1)
        assert e.value.args[0]=="Density is 0 or less"
    
         
    def test_get_textline__topo_as_Series(self, rv_single):
        """
        sample as dict, topo as pd.Series
        - raise exception for incomplete / invalid topo data
        - raise exception if shielding is in conflict with the topo/sample data 
        """
        clips = rv_single.clip_all_rasters()
        topostats, summary = riversand.get_topostats(clips, bins=500, centroid='from_clipped')

        topo = topostats.iloc[0]
        assert isinstance(topo, pd.Series)
        #assert np.round(topo['wt'], 2)==0.01 # irrelevant for get_textline()
        assert isinstance(topo.index, pd.Index)
        
        sample = {'N':1e5, 'delN':1e2}
        
        # default values
        outstr = riversand.get_textline(sample, topo, shielding=1)
        assert outstr==("Test 45.59974 7.34807 877.024 std 0 2.65 1.00000 0.00000 2010 ; "+
                        "Test Be-10 quartz 100000 100 07KNSTD ;")
        
        t2 = topostats.drop(columns=['shielding'])
        topo = t2.iloc[0]
        with pytest.raises(ValueError) as e:
            riversand.get_textline(sample, topo, shielding='topo')
        assert e.value.args[0]=="Invalid shielding 'topo' (no shielding defined in 'topo')"
 
        t2 = topostats.drop(columns=['lat'])
        topo = t2.iloc[0]
        with pytest.raises(ValueError) as e:
            riversand.get_textline(sample, topo, shielding=1)
        assert e.value.args[0]=="Latitude is not defined"
        
    
    def test_get_textline__noQtzCorr(self, rv_single):
        """ clipping quartz-free has no impact """
    
        clips = rv_single.clip_all_rasters()    
        data = {'N': 1.2E+5, 'delN': 1.2E+4,
                'density' : "1.8", # set to 1.8
                'shielding' : 0.33, # ignored by shielding='topo'
                'erate' : 1.11, # IS NOT RESET to 0
                }    
        # quartz-free not clipped
        topostats, summary = riversand.get_topostats(clips, bins=500, centroid='from_clipped')
        summary.update({'elevation': summary['elevLo']})    
        textline = riversand.get_textline(data, summary, shielding='topo')
        
        assert textline==("Test 45.59974 7.34807 2300.000 std 0 1.8 0.92309 1.11000 2010 ; "+
                          "Test Be-10 quartz 120000 12000 07KNSTD ;")

    def test_get_textline__QtzCorr(self, rv_single):
        """ clipping quartz-free has no impact """
    
        clips = rv_single.clip_all_rasters()
        clips = riversand.eliminate_quartzfree(clips)
        data = {'N': 1.2E+5, 'delN': 1.2E+4,
                'density' : "1.8", # set to 1.8
                'shielding' : 0.33, # ignored by shielding='topo'
                'erate' : 1.11, # IS NOT RESET to 0
                }
        # clip quartz-free - note differences in textline
        topostats, summary = riversand.get_topostats(clips, bins=500, centroid='from_clipped')
        summary.update({'elevation': summary['elevLo']})    
        textline = riversand.get_textline(data, summary, shielding='topo')
        
        assert textline==("Test 45.61686 7.34994 2172.000 std 0 1.8 0.92685 1.11000 2010 ; "+
                          "Test Be-10 quartz 120000 12000 07KNSTD ;")
        
        # input does not get modified
        assert data['density']=='1.8'
        assert len(data.keys())==5
      
    def test_get_textline__DataFrame(self, rv_single):
        """ input data as pd.DataFrame """
        
        clips = rv_single.clip_all_rasters()
        topostats, summary = riversand.get_topostats(clips, bins=500, centroid='from_clipped')
        
        # single line df's
        topo = topostats.loc[topostats['bin']==500]
        sample = rv_single.samples.loc[rv_single.samples['name']=='DB02']
        assert isinstance(topo, pd.DataFrame)
        assert len(topo)==1
        assert isinstance(sample, pd.DataFrame)
        assert len(sample)==1
        
        textline = riversand.get_textline(sample, topo, shielding='topo')
        assert textline==("DB02 45.59974 7.34807 877.024 std 0 2.7 0.90576 0.00000 2000 ; "+
                          "DB02 Be-10 quartz 23500 1400 07KNSTD ;")
        
        # topo 2-line df / sample multi-line df
        topo2 = topostats.loc[topostats['bin']<2000]
        sample2 = rv_single.samples
        assert isinstance(topo2, pd.DataFrame)
        assert len(topo2)==3
        assert isinstance(sample2, pd.DataFrame)
        assert len(sample2)==10
        
        with pytest.raises(TypeError) as e:
            riversand.get_textline(sample2, topo, shielding='topo')
        assert e.value.args[0]=="get_textline() argument 'sample' must be length 1"
        
        with pytest.raises(TypeError) as e:
            riversand.get_textline(sample, topo2, shielding='topo')
        assert e.value.args[0]=="get_textline() argument 'topo' must be length 1"


def test_results_from_xml():
    """ split this into three functions for age, erosion rate and NofE and move into separate module """
    # riversand.results_from_xml()
    pass

#@mark.skip(reason="url request, slow")
class Test_get_erates_from_server():
    """
    get_erates_from_server() -> (pd.DataFrame, str, dict) # g/cm2/yr
    - takes textline or textblock of any valid calculator input
    - raises TypeError for 'textline' None or ''
    - returns empty df and server diagnostics if textline has errors
    - returns empty df and diagnostics="No response from {}" if server is down

    """
    def test_get_erates_from_server1(self):
        """ valid input, single textline """
        textline=("Test 45.62 7.35 1938 std 0 1.8 0.93 0. 2010 ; "+
                  "Test Be-10 quartz 120000 12000 07KNSTD ;")
        df, diagn, ver = riversand.get_erates_from_server(textline)
        assert len(df)==1
        assert set(df.keys())=={'E_LSDn', 'E_Lm', 'E_St',
                                'delE_ext_LSDn', 'delE_ext_Lm', 'delE_ext_St',
                                'delE_int_LSDn', 'delE_int_Lm', 'delE_int_St',
                                'density', 'name', 'nuclide'}
        
    def test_get_erates_from_server2(self):
        """ input is (valid) textblock """
        textblock=("Test 45.62 7.35 1938 std 0 1.8 0.93 0. 2010 ; "+
                   "Test Be-10 quartz 120000 12000 07KNSTD ;"+
                   "PH1 45.62 7.35 1938 std 0 1.8 0.93 0. 2010 ; "+
                   "PH1 Be-10 quartz 120000 12000 KNSTD ;" +
                   "PH1 Al-26 quartz 240000 12000 KNSTD ;")
        df, diagn, ver = riversand.get_erates_from_server(textblock)
        assert diagn=="No diagnostics"
        assert len(df)==3
        
    def test_get_erates_from_server_Invalid(self):
        """ diagnostics returned from server for invalid input """
        textline=("Test 45.62 7.35 1938 std 0 1.8 0.93 0. 2010 ; "+
                  "PH-1 Be-10 quartz 120000 12000 07KNSTD ;")
        df, diagn, ver = riversand.get_erates_from_server(textline)
        assert diagn=="validate_v3_input.m: Can't match line 2 to sample"
        assert df.empty
        
    @mark.xfail(reason="tests behaviour if server unresponsive / no internet connection")
    def test_get_erates_from_server_noInternet(self):
        textline=("Test 45.62 7.35 1938 std 0 1.8 0.93 0. 2010 ; "+
                  "Test Be-10 quartz 120000 12000 07KNSTD ;")
        
        df, diagn, ver = riversand.get_erates_from_server(textline)
        assert diagn=="No response from http://stoneage.hzdr.de/cgi/matweb"
        assert df.empty
        assert ver=={}
        assert set(df.keys())=={'E_LSDn', 'E_Lm', 'E_St',
                                'delE_ext_LSDn', 'delE_ext_Lm', 'delE_ext_St',
                                'delE_int_LSDn', 'delE_int_Lm', 'delE_int_St',
                                'density', 'name', 'nuclide'}

    
    
class Test_get_E():
    """
    get_E() -> dict # cm/yr
    - calls get_erates_from_server()
    - takes only a single sample!!
    - raises RuntimeError if more than one sampple
    - raises RuntimeError if diagnostics!="No diagnostics"
       (saturated, N very small, no server response)
    - returns dict {'St', 'Lm', 'LSDn'}
    
    """
    
    @mark.xfail(reason="tests behaviour if server unresponsive / no internet connection")
    def test_get_E_noInternet(self):
        textline=("Test 45.62 7.35 1938 std 0 1.8 0.93 0. 2010 ; "+
                  "Test Be-10 quartz 120000 12000 07KNSTD ;")
        with pytest.raises(RuntimeError) as e:
            riversand.get_E(textline)
        assert e.value.args[0]=="get_E() : No response from http://stoneage.hzdr.de/cgi/matweb"
        
    def test_get_E_valid(self):
        """ valid textline """
        textline=("Test 45.62 7.35 1938 std 0 1.8 0.93 0. 2010 ; "+
                  "Test Be-10 quartz 120000 12000 07KNSTD ;")
        out = riversand.get_E(textline)
        assert np.round(out['St'], 5)==0.01446
        assert np.round(out['Lm'], 5)==0.01473
        assert np.round(out['LSDn'], 5)==0.01497
        assert 'diagnostics' not in out.keys()
        
    def test_get_E_RuntimeError1(self):
        """ two-sample textline """
        textline=("Test 45.62 7.35 1938 std 0 1.8 0.93 0. 2010 ; "+
                  "Test Be-10 quartz 120000 12000 07KNSTD ;"+
                  "Test Al-26 quartz 120000 12000 KNSTD ;")
        with pytest.raises(RuntimeError) as e:
            riversand.get_E(textline)
        assert e.value.args[0]==("get_E() : invalid results, probably from invalid "+
                                 "argument 'textline'; must be a single sample and nuclide")
    
    def test_get_E_RuntimeError2(self):
        """ very low N """
        textline=("Test 45.62 7.35 1938 std 0 1.8 0.93 0. 2010 ; "+
                  "Test Be-10 quartz 12 12000 07KNSTD ;")
        with pytest.raises(RuntimeError) as e:
            riversand.get_E(textline)
        assert e.value.args[0]==("get_E() : nuclide concentration too low for calculation")

    def test_get_E_RuntimeError3(self):
        """ bogus textline """
        textline=("Test 45.62 7.35 1938 std 0 1.8 0.93 0. 2010 ; "+
                  "Test Be-10 q 120000 12000 07KNSTD ;")
        with pytest.raises(RuntimeError) as e:
            riversand.get_E(textline)
        assert e.value.args[0]==("get_E() : invalid results, probably from invalid "+
                                 "argument 'textline'; must be a single sample and nuclide")
        
    def test_get_E_RuntimeError4(self):
        """ bogus textline """
        textline=("Test 45.62 7.35 1938 std 0 1.8 0.93 0. 2010 ; "+
                  "Test Be-10 quartz 1e10 12000 07KNSTD ;")
        with pytest.raises(RuntimeError) as e:
            riversand.get_E(textline)
        assert e.value.args[0]==("get_E() : sample appears to be saturated")
        
 
class Test_get_NofE_from_server():
    """
    get_NofE_from_server() -> (pd.DataFrame, str, dict) # g/cm2/yr
    - takes textline or textblock of any valid calculator input
    - raises TypeError for 'textline' None or ''
    - returns empty df and server diagnostics if textline has errors
    - returns empty df and diagnostics="No response from {}" if server is down

    """
        
    def test_get_NofE_from_server_valid(self):
        text_block=("A1 33.30000 72.00000 1000.000 std 0 2.65 0.99000 0.00000 2010 ; "+
                    "A1 Be-10 quartz 100000 100 07KNSTD ;"+
                    "A2 33.30000 72.00000 1000.000 std 0 2.65 0.99000 0.00000 2010 ; "+
                    "A2 Be-10 quartz 100000 100 07KNSTD ;"+
                    "A3 33.30000 72.00000 500.000 std 0 2.65 0.99000 0.00000 2010 ; "+
                    "A3 Be-10 quartz 100000 100 07KNSTD ;")
        
        dfN, diagn, ver = riversand.get_NofE_from_server(text_block)
        assert diagn=="No diagnostics"
        assert len(dfN)==3
        
    def test_get_NofE_from_server_invalid(self):
        text_block=("A1 33.30000 72.00000 1000.000 std 0 2.65 0.99000 0.00000 2010 ; "+
                    "A1 Be-10 quartz 100000 100 07KNSTD ;"+
                    "A2 33.30000 72.00000 1000.000 std 0 2.65 0.99000 0.00000 2010 ; "+
                    "A2 Be-10 quartz 100000 100 07KNSTD ;"+
                    "A2 33.30000 72.00000 500.000 std 0 2.65 0.99000 0.00000 2010 ; "+
                    "A2 Be-10 quartz 100000 100 07KNSTD ;")
        
        dfN, diagn, ver = riversand.get_NofE_from_server(text_block)
        assert diagn=="validate_v3_input.m: Line 4 appears to match more than one sample"
        assert dfN.empty
        assert set(dfN.keys())=={'E_cmyr', 'NpredLSDn', 'NpredLm', 'NpredSt', 'name', 'nuclide'}

    @mark.xfail(reason="tests behaviour if server unresponsive / no internet connection")
    def test_get_NofE_from_server_noInternet(self):
        text_block=("Test 45.62 7.35 1938 std 0 1.8 0.93 0. 2010 ; "+
                    "Test Be-10 quartz 120000 12000 07KNSTD ;")
        
        dfN, diagn, ver = riversand.get_NofE_from_server(text_block)
        assert diagn=="No response from http://stoneage.hzdr.de/cgi/matweb"
        assert dfN.empty
        assert ver=={}
        assert set(dfN.keys())=={'E_cmyr', 'NpredLSDn', 'NpredLm', 'NpredSt', 'name', 'nuclide'}
    
    
#@mark.skip(reason="url request, slow")
class Test_get_NofE_FullTable():
    """
    Get predicted nuclide concentrations for given topographic statistics
    for a single sample and a suite of erates.
    
    - 'erates' numeric, list or erates = guess_erates()
    - 'topostats' pd.DataFrame, pd.Series (topostats.loc[n]) or dict (summary)
        raise ValueError if invalid / missing data
        'topostats' shold have a key 'wt' or is set wt=1 if missing
    - 'sample_data' pd.Series (samples.iloc[n]) or dict ({'N': xx, 'delN': xx})
        or single-row pd.DataFrame (samples.loc[samples['name']=='DB03',:])
        raise ValueError if invalid / missing data
        raise ValueError if invalid shielding
    """
    @pytest.fixture(autouse=True)
    def _clips_topostats_summary(self, rv_single):
        clips = rv_single.clip_all_rasters()
        topostats, summary = riversand.get_topostats(clips, bins=500, centroid='from_clipped')
        self._clips = clips
        self._topostats = topostats
        self._summary = summary

    def test_get_NofE_FullTable__erate_valid1(self, rv_single):
        """ 'erates' numeric, 'sample_data' pd.Series """
        erates = 0.5
        
        dfN, diagn, ver = get_NofE_FullTable(sample_data = rv_single.samples.iloc[0],
                                           topostats = self._topostats,
                                           shielding = 'topo',
                                           erates = erates)
        assert len(dfN)==(len(self._topostats))
        assert set(dfN.keys())=={'E_cmyr', 'LSDn', 'Lm', 'NpredLSDn',
                                  'NpredLm', 'NpredSt', 'St', 'name', 'nuclide'}
        assert diagn=="No diagnostics" 
        
    def test_get_NofE_FullTable__erate_valid2(self, rv_single):
        """ 'erates' [numeric] """
        erates = [0.5]
        
        dfN, diagn, ver = get_NofE_FullTable(sample_data = rv_single.samples.iloc[0],
                                           topostats = self._topostats,
                                           shielding = 'topo',
                                           erates = erates)
        assert len(dfN)==(len(self._topostats)*len(erates))
        assert diagn=="No diagnostics" 
    
    def test_get_NofE_FullTable__erate_valid3(self, rv_single):
        """ 'erates' list of numeric """
        erates = [0.5, 0.3, 0.1]
        
        dfN, diagn, ver = get_NofE_FullTable(sample_data = rv_single.samples.iloc[0],
                                           topostats = self._topostats,
                                           shielding = 'topo',
                                           erates = erates)
        assert len(dfN)==(len(self._topostats)*len(erates))
        assert diagn=="No diagnostics" 
        
    def test_get_NofE_FullTable__erate_valid4(self, rv_single):
        """ erates = riversand.guess_erates() """
        erates = riversand.guess_erates(0.1, 0.5, N=2)
        
        dfN, diagn, ver = get_NofE_FullTable(sample_data = rv_single.samples.iloc[0],
                                           topostats = self._topostats,
                                           shielding = 'topo',
                                           erates = erates)
        assert len(erates)==2
        assert len(dfN)==(len(self._topostats)*len(erates))
        assert diagn=="No diagnostics" 
        
    def test_get_NofE_FullTable__erate_valid5(self, rv_single):
        """ erates = 0 """
        erates = 0
        
        dfN, diagn, ver = get_NofE_FullTable(sample_data = rv_single.samples.iloc[0],
                                           topostats = self._topostats,
                                           shielding = 'topo',
                                           erates = erates)
        assert len(dfN)==(len(self._topostats))
        assert diagn=="No diagnostics" 
        
        
        
    def test_get_NofE_FullTable__topostats_valid1(self, rv_single):
        """ 'topostats' single row of topostats (pd.Series); preserve 'wt'"""
        
        topostats = deepcopy(self._topostats.iloc[2])
        topostats['wt'] = 1
        assert isinstance(topostats, pd.Series)
        assert np.round(topostats['elevation'], 1)==1771.9
        assert isinstance(topostats.index, pd.Index)
        dfN, diagn, ver = get_NofE_FullTable(sample_data = {'N':1e4, 'delN':1e2},
                                             topostats = topostats,
                                             shielding = .99,
                                             erates = 0.5)
        assert len(dfN)==1 # one erate, one topostats entry
        assert set(dfN.keys())=={'E_cmyr', 'LSDn', 'Lm', 'NpredLSDn',
                                  'NpredLm', 'NpredSt', 'St', 'name', 'nuclide'}
        
        E_low_elev = np.round(dfN.loc[0, 'LSDn'], 2)
        assert E_low_elev == 2363.3 # E_LSDn for third row in topostats
        
    def test_get_NofE_FullTable__topostats_valid2(self, rv_single):
        """ 'topostats' single row of topostats (pd.Series); preserve 'wt'"""
        
        topostats = deepcopy(self._topostats.iloc[-3])
        topostats['wt'] = 1
        assert np.round(topostats['elevation'], 1)==3188.2
        dfN, diagn, ver = get_NofE_FullTable(sample_data = {'N':1e4, 'delN':1e2},
                                            topostats = topostats,
                                            shielding = .99,
                                            erates = 0.5)
        E_hi_elev = np.round(dfN.loc[0, 'LSDn'], 2)
        assert E_hi_elev == 6199.1 # E_LSDn for third last row in topostats
        
    def test_get_NofE_FullTable__topostats_valid3(self, rv_single):
        """ 'topostats' single row of topostats (pd.Series); preserve 'wt'"""
        
        topostats = deepcopy(self._topostats.iloc[-3])
        assert topostats['wt']<1 # cf test_get_NofE_FullTable__topostats_valid2()
        dfN, diagn, ver = get_NofE_FullTable(sample_data = {'N':1e4, 'delN':1e2},
                                            topostats = topostats,
                                            shielding = .99,
                                            erates = 0.5)
        E_wt = np.round(dfN.loc[0, 'LSDn'], 2)
        assert E_wt == 830.5 # E_LSDn weighted by 'wt'
        assert np.round(topostats['wt'], 3)==0.134
        
        
        
    def test_get_NofE_FullTable__topostats_dict1(self, rv_single):
        """ 'topostats' dict; preserve 'wt' """
                
        # summary has no key 'elevation'
        summary = self._summary
        with pytest.raises(ValueError) as e:
            get_NofE_FullTable(sample_data = {'N':1e4, 'delN':1e2},
                               topostats = summary,
                               shielding = .99,
                               erates = 0.5)
        assert e.value.args[0]==("Cannot get valid input for the online calculator:\n"+
                                  "   Elevation is not defined")
        
    def test_get_NofE_FullTable__topostats_dict2(self, rv_single):
        """ 'topostats' dict; preserve 'wt' """
        
        summary = self._summary
        summary['elevation'] = summary['elevLo']
        dfN, diagn, ver = get_NofE_FullTable(sample_data = {'N':1e4, 'delN':1e2},
                                             topostats = summary,
                                             shielding = .99,
                                             erates = 0.5)
        assert 'wt' not in summary.keys() # 'wt' set to 1.0 by the function
        assert np.round(dfN.loc[0, 'LSDn'], 2)==3434.3

    def test_get_NofE_FullTable__topostats_dict3(self, rv_single):
        """ 'topostats' dict; preserve 'wt' """
            
        # custom 'wt' in summary 
        summary = self._summary
        summary['elevation'] = summary['elevLo']
        summary['wt'] = .5
        dfN, diagn, ver = get_NofE_FullTable(sample_data = {'N':1e4, 'delN':1e2},
                                             topostats = summary,
                                             shielding = .99,
                                             erates = 0.5)
        assert np.round(dfN.loc[0, 'LSDn'], 2)==1717.15
        # with 'wt'=0.5 the resulting erate is 50% lower
        #     compared to test_get_NofE_FullTable__topostats_dict3()
        
        
        
    def test_get_NofE_FullTable__standardization(self, rv_single):
        """ standardization is ignored / reset """
        dfN, diagn, ver = get_NofE_FullTable(sample_data = {'N':1e4, 'delN':1e2},
                                             topostats = self._topostats.iloc[3],
                                             shielding = .99,
                                             erates = 0.5)
        val1 = np.round(dfN.loc[0, 'LSDn'], 5)
        
        dfN, diagn, ver = get_NofE_FullTable(sample_data = {'N':1e4, 'delN':1e2,
                                                          'standardization': 'NIST'},
                                             topostats = self._topostats.iloc[3],
                                             shielding = .99,
                                             erates = 0.5)
        val2 = np.round(dfN.loc[0, 'LSDn'], 5)
        assert abs(1-val2/val1)<1e-5 # essentially identical
        
        
    def test_get_NofE_FullTable__sample_data(self, rv_single):
        """ 'sample_data' single-line pd.DataFrame vs. Series is all the same"""
        
        # sample_data1 : single-row DataFrame
        sample_data1 = rv_single.samples.loc[rv_single.samples['name']=='DB03',:]
        sample_data1.loc[1,'shielding']=1
        sample_data1.loc[1,'comment']=''
        dfN1, diagn, ver = get_NofE_FullTable(sample_data = sample_data1,
                                              topostats = {'lat': 33, 'long': 72, 'elevation': 0},
                                              shielding = 1,
                                              erates = 0.5)
        
        # sample_data2 : pd.Series
        sample_data2 = deepcopy(rv_single.samples.iloc[1])
        sample_data2['shielding']=1
        sample_data2['comment']=''
        dfN2, diagn, ver = get_NofE_FullTable(sample_data = sample_data2,
                                              topostats = {'lat': 33, 'long': 72, 'elevation': 0},
                                              shielding = 1,
                                              erates = 0.5)
        assert (sample_data1.iloc[0]==sample_data2).all() # input sample_data  is same
        assert (dfN1==dfN2).all().all() # result is same
        
    def test_get_NofE_FullTable__sample_data_Exception(self, rv_single):
        """ sample_data multi-line pd.DataFrame """
        
        sample_data = rv_single.samples.loc[rv_single.samples['name']=='DB04',:]
        assert len(sample_data)==2
        with pytest.raises(TypeError) as e:
            get_NofE_FullTable(sample_data = sample_data,
                               topostats = {'lat': 33, 'long': 72, 'elevation': 0},
                               shielding = 1,
                               erates = 0.5)
        assert e.value.args[0]==("get_NofE_FullTable() argument 'sample_data' "+
                                 "must be dict, pandas Series or single-row pandas Dataframe")    
        
                
    def test_get_NofE_FullTable__erate_Exception(self, rv_single):
        """ invalid value in 'erates' """
        
        with pytest.raises(ValueError) as e:
            get_NofE_FullTable(sample_data = rv_single.samples.iloc[0],
                               topostats = self._topostats,
                               shielding = 'topo',
                               erates = [-1, 0.1, 0.2])
        assert e.value.args[0]==("Cannot get valid input for the online calculator:\n"+
                                  "   Erosion rate is less than 0")
        
    def test_get_NofE_FullTable__sample_Exception(self, rv_single):
        """ invalid value in 'sample_data' """
        
        with pytest.raises(ValueError) as e:
            get_NofE_FullTable(sample_data = {'N':1e4, 'delN':1e2,
                                              'density': 0},
                               topostats = self._topostats,
                               shielding = 'topo',
                               erates = 0.5)
        assert e.value.args[0]==("Cannot get valid input for the online calculator:\n"+
                                  "   Density is 0 or less")
        
    def test_get_NofE_FullTable__shielding_Exception1(self, rv_single):
        """ invalid 'shielding' """
        
        with pytest.raises(ValueError) as e:
            get_NofE_FullTable(sample_data = {'N':1e4, 'delN':1e2},
                               topostats = {'lat': 33, 'long': 72, 'elevation': 0},
                               shielding = 'topo',
                               erates = 0.5)
        assert e.value.args[0]==("Cannot get valid input for the online calculator:\n"+
                                  "   Invalid shielding 'topo' (no shielding defined in 'topo')")
        
    def test_get_NofE_FullTable__shielding_Exception2(self, rv_single):
        """ invalid 'shielding' """
        
        with pytest.raises(ValueError) as e:
            get_NofE_FullTable(sample_data = {'N':1e4, 'delN':1e2},
                               topostats = {'lat': 33, 'long': 72, 'elevation': 0},
                               shielding = 'sample',
                               erates = 0.5)
        assert e.value.args[0]==("Cannot get valid input for the online calculator:\n"+
                                  "   Invalid shielding 'sample' (no shielding defined in 'sample')")


#@mark.skip(reason="url request, slow")
class Test_get_NofE():
    """
    Wrapper for get_NofE_FullTable:
    Predicted nuclide concentrations for a given catchment topography
    for a single sample and a suite of erosion rates.
    
    - 'topostats' must be pd.DataFrame
    - 'sample_data' pd.Series or dict or single-line pd.DataFrame
    - 'sample_data' and 'topostats' are validated (missing values set to default)
        raises ValueError "Invalid values in..." if bogus
    - raises RuntimeError if diagnostics is anything but "No diagnostics"
        (saturated is not an issue, only possibility is no url-response)
    - raises RunTimeError if get_NofE_FullTable() returns an empty DataFrame
        (not sure why this would ever happen)
    """
    @pytest.fixture(autouse=True)
    def _clips_topostats_summary(self, rv_single):
        clips = rv_single.clip_all_rasters()
        topostats, summary = riversand.get_topostats(clips, bins=500, centroid='from_clipped')
        self._clips = clips
        self._topostats = topostats
        self._summary = summary
    
    @mark.xfail(reason="tests behaviour if server unresponsive / no internet connection")
    def test_get_NofE_noInternet(self, rv_single):
        sample = {'N': 1e5, 'delN': 1e2}
        with pytest.raises(RuntimeError) as e:
            riversand.get_NofE(sample, self._topostats, shielding='topo', erates=0.5)
        assert e.value.args[0]=="get_NofE() : No response from http://stoneage.hzdr.de/cgi/matweb"
        
    def test_get_NofE_valid_sample1(self, rv_single):
        """ 'sample_data' dict """

        sample = {'N': 1e5, 'delN': 1e2}
        N = riversand.get_NofE(sample, self._topostats, shielding='topo', erates=0.5)
        assert len(N)==1
        
    def test_get_NofE_valid_sample2(self, rv_single):
        """ 'sample_data' pd.Series """

        sample = rv_single.samples.iloc[4] # has shielding data
        N = riversand.get_NofE(sample, self._topostats, shielding='sample', erates=0.5)
        assert len(N)==1
        
    def test_get_NofE_valid_sample3(self, rv_single):
        """ 'sample_data' single-row pd.DataFrame """

        sample = rv_single.samples.loc[rv_single.samples['name']=='DB05',:] # has shielding data
        N = riversand.get_NofE(sample, self._topostats, shielding='sample', erates=0.5)
        assert len(N)==1
    
    def test_get_NofE_valid_topostats(self, rv_single):
        """ 'topostats' pd.DataFrame with a single entry """
        topostats = self._topostats.loc[self._topostats['bin']==500]
        assert len(topostats)==1
        sample = rv_single.samples.iloc[4]
        N = riversand.get_NofE(sample, topostats, shielding='topo', erates=[0.2, 0.5])
        assert len(N)==2
        
    
    def test_get_NofE_shieldingException1(self, rv_single):
        """ invalid numeric 'shielding' """
        
        sample = {'N':1e5, 'delN': 1e3} 
        with pytest.raises(ValueError) as e:
            riversand.get_NofE(sample, self._topostats, shielding=1.2, erates=0.5)
        assert e.value.args[0]==("Cannot get valid input for the online calculator:\n"+
                                 "   Shielding is out of bounds (0..1)")

    def test_get_NofE_shieldingException2(self, rv_single):
        """ invalid 'shielding' """
        
        sample = {'N':1e5, 'delN': 1e3} 
        with pytest.raises(ValueError) as e:
            riversand.get_NofE(sample, self._topostats, shielding='sample', erates=0.5)
        assert e.value.args[0]=="shielding='sample' but no shielding factor defined in 'sample_data'"

    
    def test_get_NofE_shieldingException3(self, rv_single):
        """ invalid 'shielding' """
        
        sample = rv_single.samples.iloc[0] 
        with pytest.raises(ValueError) as e:
            riversand.get_NofE(sample, self._topostats, shielding='sample', erates=0.5)
        assert e.value.args[0]=="shielding='sample' but no shielding factor defined in 'sample_data'"
        
    def test_get_NofE_shieldingException4(self, rv_single):
        """ invalid 'shielding' """
        
        rv_single.shielding = None
        clips = rv_single.clip_all_rasters()
        topostats, summary = riversand.get_topostats(clips, bins=500, centroid='from_clipped')
        sample = rv_single.samples.iloc[0] 
        with pytest.raises(ValueError) as e:
            riversand.get_NofE(sample, topostats, shielding='topo', erates=0.5)
        assert e.value.args[0]=="shielding='topo' but no shielding factor defined in 'topostats'"
        
    def test_get_NofE_sampledataException1(self, rv_single):
        """ 'sample_data' multi-line pd.DataFrame """
        
        sample = rv_single.samples # passing the whole 'samples'
        with pytest.raises(TypeError) as e:
            riversand.get_NofE(sample, self._topostats, shielding='sample', erates=0.5)
        assert e.value.args[0]==("get_NofE() argument 'sample_data' must be dict, "+
                                 "pandas Series or single-row pandas Dataframe")
        
    def test_get_NofE_topostatsException1(self, rv_single):
        """ 'topostats' dict """
        
        # passing summary instead of topostats
        sample = rv_single.samples.iloc[0]
        with pytest.raises(TypeError) as e:
            riversand.get_NofE(sample, self._summary, shielding='sample', erates=0.5)
        assert e.value.args[0]==("get_NofE() argument 'topostats' must be "+
                                 "pandas DataFrame")

    def test_get_NofE_topostatsException2(self, rv_single):
        """ empty 'topostats' """
        
        sample = rv_single.samples.iloc[0]
        topostats = self._topostats.loc[self._topostats['bin']<500]
        assert len(topostats)==0
        with pytest.raises(ValueError) as e:
            riversand.get_NofE(sample, topostats, shielding='sample', erates=0.5)
        assert e.value.args[0]=="Invalid values in 'topostats' or empty table"
        
        

class Test_guess_erates():
    @pytest.fixture(autouse=True)
    def _E(self):
        self._Lo = {'St': 0.0624, 'Lm': 0.0640, 'LSDn': 0.0667}
        self._Hi = {'St': 0.1011, 'Lm': 0.1043, 'LSDn': 0.1136}

    def test_guess_erates_numerical1(self):
        """ two numerical values """
        out = riversand.guess_erates(self._Lo['Lm'], self._Hi['Lm'])
        assert len(out)==6
        assert np.round(out[0], 5)==0.06400
        assert np.round(out[-1], 5)==0.10430
        
    def test_guess_erates_numerical2(self):
        """ two numerical values, unnecessary kwarg 'scaling' """
        out = riversand.guess_erates(self._Hi['Lm'], self._Lo['Lm'], N=7, scaling='bogus')
        assert len(out)==7
        assert np.round(out[0], 5)==0.06400
        assert np.round(out[-1], 5)==0.10430
        
    def test_guess_erates_numerical3(self):
        """ one  numerical values, kwarg 'N' """
        out = riversand.guess_erates(0.08, N=5)
        assert len(out)==5
        assert np.round(out[0], 5)==0.04
        assert np.round(out[-1], 5)==0.16
        
    def test_guess_erates_numerical4(self):
        """ one  numerical value one dict, kwarg 'N' """
        with pytest.raises(ValueError) as e:
            riversand.guess_erates(self._Lo, self._Hi['Lm'], N=5)
        assert e.value.args[0]=="guess_erates() invalid arguments"
    
    
    def test_guess_erates_dict1(self):
        """ two dicts """
        out = riversand.guess_erates(self._Hi, self._Lo, scaling='Lm')
        assert len(out)==6
        assert np.round(out[0], 5)==0.06400
        assert np.round(out[-1], 5)==0.10430
        
    def test_guess_erates_dict2(self):
        """ one dict, no 'scaling' """
        with pytest.raises(ValueError) as e:
            riversand.guess_erates(self._Hi)
        assert e.value.args[0]=="guess_erates() invalid keyword argument scaling 'None'"
        
    def test_guess_erates_dict3(self):
        """ one dict, bogus 'scaling' """
        with pytest.raises(ValueError) as e:
            riversand.guess_erates(self._Hi, scaling='bogus')
        assert e.value.args[0]=="guess_erates() invalid keyword argument scaling 'bogus'"
        
    def test_guess_erates_dict4(self):
        """ one dict, valid 'scaling' and 'N' """
        out = riversand.guess_erates(self._Hi, scaling='Lm', N=5)
        assert len(out)==5
        assert np.round(out[0], 5)==0.05215
        assert np.round(out[-1], 5)==0.2086
        
    def test_guess_erates_dict5(self):
        """ one dict, valid 'scaling' and 'N' """
        out = riversand.guess_erates({'Lm': 0.08}, scaling='Lm', N=5)
        assert len(out)==5
        assert np.round(out[0], 5)==0.04
        assert np.round(out[-1], 5)==0.16
        
        

class Test_poly_E_results():
    """
    - sample, topo and shielding data are validated by get_NofE(),
        which raises ValueError, TypeError, RuntimeError
    - raise Exception if erates less than 4 values (override with strict=False)
    - list of errors; if "minE too high" or "maxE too low" then just single-element list
    
    """
    
    @pytest.fixture(autouse=True)
    def _polyE(self, rv_single):
        self._samples = rv_single.samples
        clips = rv_single.clip_all_rasters()
        topostats, summary = riversand.get_topostats(clips, bins=500, centroid='from_clipped')
        self._topostats = topostats
        self._summary = summary
        
    def test_polyE_results__valid(self, rv_single):
        """ standard input """
        erates = [0.073, 0.077, 0.082, 0.087, 0.093, 0.098]
        E, delE, NofE, err = riversand.poly_E_results(self._samples.iloc[0],
                                                      self._topostats,
                                                      shielding=1,
                                                      erates=erates,
                                                      scaling='LSDn')
        assert np.round(E, 4)==0.0864
        
    def test_polyE_results__TypeError1(self, rv_single):
        """ TypeError raised by get_NofE() """
        erates = [0.073, 0.077, 0.082, 0.087, 0.093, 0.098]
        
        with pytest.raises(TypeError) as e:
            riversand.poly_E_results(self._samples,
                                     self._topostats,
                                     shielding=1, erates=erates, scaling='LSDn')
        assert e.value.args[0]==("get_NofE() argument 'sample_data' must be "+
                                 "dict, pandas Series or single-row pandas Dataframe")
    
    def test_polyE_results__TypeError2(self, rv_single):
        """ TypeError raised by get_NofE() """
        erates = [0.073, 0.077, 0.082, 0.087, 0.093, 0.098]
        
        with pytest.raises(TypeError) as e:
            riversand.poly_E_results(self._samples.iloc[0],
                                     self._summary,
                                     shielding=1, erates=erates, scaling='LSDn')
        assert e.value.args[0]==("get_NofE() argument 'topostats' must be "+
                                 "pandas DataFrame")

    def test_polyE_results__ValueError1(self, rv_single):
        """ ValueError raised by get_NofE() """
        erates = [0.073, 0.077, 0.082, 0.087, 0.093, 0.098]
        
        with pytest.raises(ValueError) as e:
            riversand.poly_E_results({'N':23500, 'delN':1400, 'nuclide': 'Be'},
                                     self._topostats,
                                     shielding=1, erates=erates, scaling='LSDn')
        assert e.value.args[0]==("Invalid values in 'sample_data'")

    def test_polyE_results__ValueError2(self, rv_single):
        """ ValueError raised by get_NofE() """
        erates = [0.073, 0.077, 0.082, 0.087, 0.093, 0.098]
        
        with pytest.raises(ValueError) as e:
            riversand.poly_E_results({'N':23500, 'delN':1400},
                                     self._topostats,
                                     shielding='sample', erates=erates, scaling='LSDn')
        assert e.value.args[0]==("shielding='sample' but no shielding factor "+
                                 "defined in 'sample_data'")
    
    def test_polyE_results__RuntimeError_erates(self, rv_single):
        """ too few values in 'erates'"""
        erates = [0.073, 0.077, 0.098]
        
        with pytest.raises(RuntimeError) as e:
            riversand.poly_E_results(self._samples.iloc[0],
                                     self._topostats,
                                     shielding=1, erates=erates, scaling='LSDn')
        assert e.value.args[0]==("poly_E_results() : argument 'erates' should "+
                                 "have at least 4 values; use argument "+
                                 "'strict=False' to override")
    
    def test_polyE_results__RuntimeError_override(self, rv_single):
        """ too few values in 'erates' but 'strict=False'"""
        erates = [0.073, 0.077, 0.098]
        
        E, delE, NofE, err = riversand.poly_E_results(self._samples.iloc[0],
                                     self._topostats, shielding=1,
                                     erates=erates, scaling='LSDn',
                                     strict=False)
        assert np.round(E, 4)==0.0863
    
    
    def test_polyE_results__error1(self, rv_single):
        """ list of errors"""
        erates = [0.073, 0.077, 0.082, 0.087, 0.093, 0.098]
        
        E, delE, NofE, err = riversand.poly_E_results({'N': 8000, 'delN': 80},
                                     self._topostats, shielding=1,
                                     erates=erates, scaling='LSDn')
        assert err==["maxE too low"]
        assert np.isnan(E)
    
    
class Test_uncert_on_E():
    """ tuple delE-, delE+ """
        
    @pytest.fixture(autouse=True)
    def _clips_topostats_summary(self, rv_single):
        clips = rv_single.clip_all_rasters()
        topostats, summary = riversand.get_topostats(clips, bins=500, centroid='from_clipped')
        sample_data = rv_single.samples.iloc[0]
        erates = [0.06, 0.07, 0.08, 0.09, 0.1]
        E, delE, NofE, err = (riversand.
                      poly_E_results(sample_data, topostats, 
                                     shielding=1, erates=erates, scaling='Lm'))
        self._samples = rv_single.samples
        self._NofE = NofE
        self._E = E
        self._delE = delE
        
    def test_uncert_on_E_dict(self):
        """ 'sample_data' dict """
        sample_data = {'N': 23500, 'delN': 1400}
        delE = riversand.calc.uncert_on_E(self._E, self._NofE,
                                          sample_data)
        assert delE==self._delE
        
    def test_uncert_on_E_Series(self):
        """ 'sample_data' pd.Series """
        sample_data = self._samples.iloc[0]
        assert isinstance(sample_data, pd.Series)
        delE = riversand.calc.uncert_on_E(self._E, self._NofE,
                                          sample_data)
        assert delE==self._delE

    def test_uncert_on_E_DataFrame(self):
        """ 'sample_data' single-row pd.DataFrame """
        sample_data = self._samples.loc[self._samples['name']=='DB02',:]
        assert isinstance(sample_data, pd.DataFrame)
        delE = riversand.calc.uncert_on_E(self._E, self._NofE,
                                          sample_data)
        
        assert delE==self._delE
        
    def test_uncert_on_E_exception(self):
        """ 'sample_data' multi-row pd.DataFrame """
        sample_data = self._samples.loc[self._samples['name']=='DB04',:]
        assert len(sample_data)==2
        with pytest.raises(TypeError) as e:
            riversand.calc.uncert_on_E(self._E, self._NofE,
                                              sample_data)
            

def test_get_RMSE():
    """ return number, function of NofE only """
    pass

def test_get_version():
    """ dict version infor from server """
    # xfail: raises RuntimeError "get_version() : {diagnostics} if noInternet
    
    

    
    
    
    
