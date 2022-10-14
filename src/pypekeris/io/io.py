# -*- coding: utf-8 -*-
"""
Python module to read the .D binary data files

.. module:: py ocean acoustics data

:author:
    Gil Averbuch (gil.averbuch@whoi.edu)

:copyright:
    Gil Averbuch

:license:
    This code is distributed under the terms of the
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.en.html)
"""

import matplotlib.pyplot as plt
plt.rcParams['font.size'] = '16'
plt.rcParams['figure.dpi'] = 125
plt.rcParams['figure.facecolor'] = 'white'

import numpy as np
import pandas as pd


# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

class pekeris_():
    def __init__(self, c0,  c1,  c2, d, f, nq):
        '''
        This function initializes the pekeris class


        parameters
        ----------
            c0 = reference velocity. default 1500 m/s
            c1 = velocity layer 1. default 1500 m/s
            c2 = velocity layer 2. default 2000 m/s
            d  = depth of layer 1. default 1000 m
            f = frequency. default is 1 Hz
            nq = number of wavenumbers. default is 1e3

        Returns
        -------
        Header: header structure containing the parameters names and values (Pandas DataFrame)
        Waveforms: Obspy Stream object
        '''
        
        self.c0 = c0
        self.c1 = c1
        self.c2 = c2
        self.d  = d 
        self.freq = f
        self.nq = np
        
        return
    
    
