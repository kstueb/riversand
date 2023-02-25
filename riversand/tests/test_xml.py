#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 15:16:40 2023

@author: Konstanze Stübner, kstueb@gmail.com

"""

import pytest
from pytest import mark

import os
import sys

import numpy as np
import pandas as pd

import urllib.parse
import urllib.request
import xml.etree.ElementTree as et


import riversand
from riversand import xml
from riversand.params import url

@pytest.fixture()
def textlines():
    textline0 = ( # incomplete textline
        "T1 33.30000 72.10000 1000.000 std 0 2.65 0.95000 0.00000 2010 ; ")
        
    textline1 = ( # single textline
        "T1 33.30000 72.10000 1000.000 std 0 2.65 0.95000 0.00000 2010 ; "+
        "T1 Be-10 quartz 100000 1000 07KNSTD ; ")
    
    textline2 = ( # two locations, three samples
        "T1 33.30000 72.10000 1000.000 std 0 2.65 0.95000 0.00000 2010 ; "+
        "T1 Be-10 quartz 100000 1000 07KNSTD ; "+
        "T1 Al-26 quartz 1000 1000 KNSTD ; " +
        "T2 33.30000 72.10000 1000.000 std 0 2.65 0.95000 0.00000 2010 ; "+
        "T2 Be-10 quartz 100000 1000 07KNSTD ; ")
    
    textline2Fa = ( # one invalid entry
        "T1 33.30000 72.10000 1000.000 std 0 2.65 0.95000 0.00000 2010 ; "+
        "T1 Be-10 quartz 100000 1000 07KNSTD ; "+
        "T1 Al-26 feldspar 1000 1000 KNSTD ; " +
        "T2 33.30000 72.10000 1000.000 std 0 2.65 0.95000 0.00000 2010 ; "+
        "T2 Be-10 quartz 100000 1000 07KNSTD ; ")
    
    textline2Fb = ( # one invalid entry
        "T1 33.30000 72.10000 1000.000 std 0 2.65 0.95000 0.00000 2010 ; "+
        "T1 Be-10 quartz 100000 1000 07KNSTD ; "+
        "T1 Al-26 quartz 1000 0 KNSTD ; " +
        "T2 33.30000 72.10000 1000.000 std 0 2.65 0.95000 0.00000 2010 ; "+
        "T2 Be-10 quartz 100000 1000 07KNSTD ; ")

    textlineNofE = ( # suite of erates
        "T1 33.30000 72.10000 1000.000 std 0 2.65 0.95000 0.01000 2010 ; "+
        "T1 Be-10 quartz 100000 1000 07KNSTD ; "+
        "T2 33.30000 72.10000 1000.000 std 0 2.65 0.95000 0.02000 2010 ; "+
        "T2 Be-10 quartz 100000 1000 07KNSTD ; "+
        "T3 33.30000 72.10000 1000.000 std 0 2.65 0.95000 0.03000 2010 ; "+
        "T3 Be-10 quartz 100000 1000 07KNSTD ; "+
        "T4 33.30000 72.10000 1000.000 std 0 2.65 0.95000 0.04000 2010 ; "+
        "T4 Be-10 quartz 100000 1000 07KNSTD ; ")
    
    textlines = {
        'TL0' : textline0,
        'TL1' : textline1,
        'TL2' : textline2,
        'TL2Fa' : textline2Fa,
        'TL2Fb' : textline2Fb,
        'TLNofE' : textlineNofE
        }
    return textlines

@pytest.fixture()
def result_xml_NofE(textlines):
    result_xml = {}
    for key, textblock in textlines.items():
        request_vars = {"mlmfile" : "NofE_input_v3", "text_block" : textblock}
        form_data = urllib.parse.urlencode(request_vars).encode('ascii') # encode request
        result = urllib.request.urlopen(url, form_data) # send request
        result_xml[key] = result.read()
    return result_xml

@pytest.fixture()
def result_xml_E(textlines):
    result_xml = {}
    for key, textblock in textlines.items():
        request_vars = {"mlmfile" : "erosion_input_v3",
                        "reportType" : "XML",
                        "resultType" : "long",
                        "summary" : "no",
                        "plotFlag" : "no",
                        "text_block" : textblock }
        form_data = urllib.parse.urlencode(request_vars).encode('ascii') # encode request
        result = urllib.request.urlopen(url, form_data) # send request
        result_xml[key] = result.read()
    return result_xml

    
def test_read_erates0(result_xml_E):
    result_xml = result_xml_E['TL0']
    df, diagn, version = xml.read_erates(result_xml)
    assert len(df)==0
    assert diagn=="validate_v3_input.m: Only one line of data - nothing to calculate"
    assert set(df.keys())=={'E_LSDn', 'E_Lm', 'E_St',
                            'delE_ext_LSDn', 'delE_ext_Lm', 'delE_ext_St',
                            'delE_int_LSDn', 'delE_int_Lm', 'delE_int_St',
                            'density', 'name', 'nuclide'}
    assert version=={} # empty
    
def test_read_erates1(result_xml_E):
    result_xml = result_xml_E['TL1']
    df, diagn, version = xml.read_erates(result_xml)
    assert len(df)==1
    assert diagn=='No diagnostics'
    assert set(df.keys())=={'E_LSDn', 'E_Lm', 'E_St',
                            'delE_ext_LSDn', 'delE_ext_Lm', 'delE_ext_St',
                            'delE_int_LSDn', 'delE_int_Lm', 'delE_int_St',
                            'density', 'name', 'nuclide'}
    assert not version=={} # not empty
    
def test_read_erates2(result_xml_E):
    result_xml = result_xml_E['TL2']
    df, diagn, version = xml.read_erates(result_xml)
    assert len(df)==3
    assert diagn=='No diagnostics'
    assert set(df.keys())=={'E_LSDn', 'E_Lm', 'E_St',
                            'delE_ext_LSDn', 'delE_ext_Lm', 'delE_ext_St',
                            'delE_int_LSDn', 'delE_int_Lm', 'delE_int_St',
                            'density', 'name', 'nuclide'}
    assert not version=={} # not empty
    
def test_read_erates2Fa(result_xml_E):
    result_xml = result_xml_E['TL2Fa']
    df, diagn, version = xml.read_erates(result_xml)
    assert len(df)==0
    assert diagn=="validate_v3_input.m: Unknown mineral for Al-26 measurement - line 3"
    assert set(df.keys())=={'E_LSDn', 'E_Lm', 'E_St',
                            'delE_ext_LSDn', 'delE_ext_Lm', 'delE_ext_St',
                            'delE_int_LSDn', 'delE_int_Lm', 'delE_int_St',
                            'density', 'name', 'nuclide'}
    assert version=={} # empty

def test_read_erates2Fb(result_xml_E):
    result_xml = result_xml_E['TL2Fb']
    df, diagn, version = xml.read_erates(result_xml)
    assert len(df)==0
    assert diagn=="validate_v3_input.m: Al-26 uncertainty less than or equal to zero - line 3"
    assert set(df.keys())=={'E_LSDn', 'E_Lm', 'E_St',
                            'delE_ext_LSDn', 'delE_ext_Lm', 'delE_ext_St',
                            'delE_int_LSDn', 'delE_int_Lm', 'delE_int_St',
                            'density', 'name', 'nuclide'}
    assert version=={} # empty

def test_read_eratesNofE(result_xml_E):
    result_xml = result_xml_E['TLNofE']
    df, diagn, version = xml.read_erates(result_xml)
    assert len(df)==4
    assert diagn=='No diagnostics'
    assert set(df.keys())=={'E_LSDn', 'E_Lm', 'E_St',
                            'delE_ext_LSDn', 'delE_ext_Lm', 'delE_ext_St',
                            'delE_int_LSDn', 'delE_int_Lm', 'delE_int_St',
                            'density', 'name', 'nuclide'}
    assert not version=={} # not empty




def test_read_NofE0(result_xml_NofE):
    result_xml = result_xml_NofE['TL0']
    df, diagn, version = xml.read_NofE(result_xml)
    assert len(df)==0
    assert diagn=="validate_v3_input.m: Only one line of data - nothing to calculate"
    assert set(df.keys())=={'E_cmyr', 'NpredLSDn', 'NpredLm', 'NpredSt', 
                            'name', 'nuclide'}
    assert version=={} # empty
    
def test_read_NofE1(result_xml_NofE):
    result_xml = result_xml_NofE['TL1']
    df, diagn, version = xml.read_NofE(result_xml)
    assert len(df)==1
    assert diagn=='No diagnostics'
    assert set(df.keys())=={'E_cmyr', 'NpredLSDn', 'NpredLm', 'NpredSt', 
                            'name', 'nuclide'}
    assert not version=={} # not empty
    
def test_read_NofE2(result_xml_NofE):
    result_xml = result_xml_NofE['TL2']
    df, diagn, version = xml.read_NofE(result_xml)
    assert len(df)==3
    assert diagn=='No diagnostics'
    assert set(df.keys())=={'E_cmyr', 'NpredLSDn', 'NpredLm', 'NpredSt', 
                            'name', 'nuclide'}
    assert not version=={} # not empty
    
def test_read_NofE02Fa(result_xml_NofE):
    result_xml = result_xml_NofE['TL2Fa']
    df, diagn, version = xml.read_NofE(result_xml)
    assert len(df)==0
    assert diagn=="validate_v3_input.m: Unknown mineral for Al-26 measurement - line 3"
    assert set(df.keys())=={'E_cmyr', 'NpredLSDn', 'NpredLm', 'NpredSt', 
                            'name', 'nuclide'}
    assert version=={} # empty

def test_read_NofE2Fb(result_xml_NofE):
    result_xml = result_xml_NofE['TL2Fb']
    df, diagn, version = xml.read_NofE(result_xml)
    assert len(df)==0
    assert diagn=="validate_v3_input.m: Al-26 uncertainty less than or equal to zero - line 3"
    assert set(df.keys())=={'E_cmyr', 'NpredLSDn', 'NpredLm', 'NpredSt', 
                            'name', 'nuclide'}
    assert version=={} # empty

def test_read_NofENofE(result_xml_NofE):
    result_xml = result_xml_NofE['TLNofE']
    df, diagn, version = xml.read_NofE(result_xml)
    assert len(df)==4
    assert diagn=='No diagnostics'
    assert set(df.keys())=={'E_cmyr', 'NpredLSDn', 'NpredLm', 'NpredSt', 
                            'name', 'nuclide'}
    assert not version=={} # not empty
