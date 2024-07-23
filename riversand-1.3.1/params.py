#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

*******************************************************************************
params.py  :  default parameters

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

"""
import os

class Params:
    out_path = './plots'

    #url = "https://hess.ess.washington.edu/cgi-bin/matweb"
    url = "http://stoneage.hzdr.de/cgi/matweb"
    
    units = {'cm/yr': 1, 'cm/kyr': 1000, 'mm/yr': 10, 'mm/kyr': 10000, 'm/Myr': 10,
              'cm/a': 1, 'cm/ka': 1000, 'mm/a': 10, 'mm/ka': 10000, 'm/Ma': 10,
              }
    
    # # see http://hess.ess.washington.edu/math/docs/v3/v3_input_explained.html
    default_values = {'name'       : 'Test',
                      'press_flag' : 'std',    # 'std', 'ant'
                      'thickness'  : 0,        # in cm
                      'density'    : 2.65,     # in g/cm3
                      'shielding'  : 1,        # 0 to 1
                      'erate'      : 0,        # in cm/yr
                      'year'       : 2010,
                      'nuclide'    : 'Be-10',  # 'Be-10', 'Al-26'
                      'mineral'    : 'quartz', # 'quartz'
                      }
    
    _default_standards = {'Be-10': '07KNSTD',
                          'Al-26': 'KNSTD',
                         }
    



    Be_stds = {"07KNSTD":1.0000,"KNSTD":0.9042,"NIST_Certified": 1.0425,
               "LLNL31000":0.8761,"LLNL10000":0.9042,"LLNL3000": 0.8644,
               "LLNL1000":0.9313,"LLNL300":0.8562,"NIST_30000": 0.9313,
               "NIST_30200":0.9251,"NIST_30300":0.9221,"NIST_30600": 0.9130,
               "NIST_27900":1.0000,"S555":0.9124,"S2007":0.9124,
               "BEST433":0.9124,"BEST433N":1.0000,"S555N":1.0000,
               "S2007N":1.0000,"STD11":1.0000,"0":0,"NIST_30500":0.9124,
               "SMDBe12": 1.0000}
    Al_stds = {"KNSTD":1.0000,"ZAL94":0.9134,"SMAL11":1.021,"0":0,
               "ZAL94N":1.0000,"ASTER":1.021,"Z92-0222":1.0000,
               "SMDAl11":1.018}
    
    _scalingmethods = ['St', 'Lm', 'LSDn']
    
    _XML_scalingmethods = ['NpredSt', 'NpredLm', 'NpredLSDn'] # returned by the online calculator
    
    # keys for the display of sample data
    all_keys = ['name',
                'lat', 'latitude', 'long', 'longitude', 'elev', 'elevation',
                'press_flag', 'thickness', 'density', 'shielding', 'erate', 'year',
                'nuclide', 'mineral', 'N', 'delN', 'standardization']
            
    @classmethod
    def set_outpath(cls, value):
        try:
            os.makedirs(value)
        except FileExistsError:
            cls.out_path = value
        except PermissionError:
            print("No permission to make directory '{}'; saving plots to '{}'."
                  .format(value, cls.out_path))
        else:
            cls.out_path = value


def set_folder(value):
    Params.set_outpath(value)

def get_folder():
    print("Plots are saved in folder '{}'.".format(Params.out_path))

def print_units():
    print("Valid units of erosion rate are: {}.".format(', '.join(Params.units.keys())))
    
def print_defaults():
    print("\nDefault values:")
    for k,v in Params.default_values.items():
        print("  {0: <10} : {1: <10}".format(k, v))
        
    print("\nDefault standardizations:")
    for k,v in Params._default_standards.items():
        print("  {0: <10} : {1: <10}".format(k, v))
    print()
    
def print_standardizations():
    print("\nBe-10 standardizations:")
    for k,v in Params.Be_stds.items():
        print("  {0: <10} : {1: <10}".format(k, v))
    print("\nAl-26 standardizations:")
    for k,v in Params.Al_stds.items():
        print("  {0: <10} : {1: <10}".format(k, v))
    print()
