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

def pekeris(zs=10, zr=10, r=5000, dr=10, c0=1500,  c1=1500,  c2=2000, d=50, dz=0.5, f=50, nq=1e3):
    '''
    This function initializes the pekeris problem


    parameters
    ----------
        zs = source depth. default 10 m
        zr = receiver depth, default 10m
        r  = receiver range, default 5000m
        dr = range interval, default 10 m
        c0 = reference velocity, default 1500 m/s
        c1 = velocity layer 1, default 1500 m/s
        c2 = velocity layer 2, default 2000 m/s
        d  = depth of layer 1, default 50 m
        dz = depth interval, default 0.5 m
        f  = frequency, default is 50 Hz
        nq = number of wavenumbers, default is 1e3

    Returns
    -------
    Header: header structure containing the parameters names and values (Pandas DataFrame)
    Waveforms: Obspy Stream object

    '''


    P = pekeris_(zs, zr, r, dr, c0,  c1,  c2, d, dz, f, nq)




    return P




# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------



