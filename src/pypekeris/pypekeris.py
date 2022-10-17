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

from .core.core import pekeris_


# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

def pekeris(zs=100, zr=100, r=5000, c0=1500,  c1=1500,  c2=2000, d=1000, dz=1, f=1, nq=1e3):
    '''
    This function initializes the pekeris problem


    parameters
    ----------
        zs = source depth. default 100 m
        zr = receiver depth, default 100m
        r  = receiver range, default 5000m
        c0 = reference velocity, default 1500 m/s
        c1 = velocity layer 1, default 1500 m/s
        c2 = velocity layer 2, default 2000 m/s
        d  = depth of layer 1, default 1000 m
        dz = depth interval, default 1 m
        f  = frequency, default is 1 Hz
        nq = number of wavenumbers, default is 1e3

    Returns
    -------
    Header: header structure containing the parameters names and values (Pandas DataFrame)
    Waveforms: Obspy Stream object

    '''


    P = pekeris_(zs, zr, r, c0,  c1,  c2, d, dz, f, nq)




    return P




# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------



