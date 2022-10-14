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

from obspy import read_inventory, read, UTCDateTime, Stream, Trace
from .input.input import read_header, read_waveforms
from .output.output import save2mseed_


# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

def read_data(file_name, records_range, bit='24bit', sensitivity=170, gain=20):
    '''
    This function teads the data from a .D 16 or 24 bit binary file


    parameters
    ----------
        file_name: path to file
        records_range: range of record sections to extact
        bit: type of binry file. Can get: '24bit', and '16bit'. Default is for 24. Still need to add the pseudo 24 bit.
        sensitivity: default is 170. optional to set to different value or specify it per channel
        gain: default sensor gain is 20

    Returns
    -------
    Header: header structure containing the parameters names and values (Pandas DataFrame)
    Waveforms: Obspy Stream object

    '''


    
    Header = read_header(file_name)
    Waveforms = read_waveforms(file_name, Header, records_range, bit, sensitivity, gain)




    return Header, Waveforms

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------











# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

def save2mseed(waveforms, dir_name='./Results/'):
    '''
    This function teads the data from a .D 24 bit binary file


    parameters
    ----------
        waveforms: obspy strem 
        dir_name: directory to save the data to in a year/network/station structure. 
                  default directory is /Results/

    Returns
    -------
    Nothing


    '''

    save2mseed_(waveforms, dir_name)

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------



