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
from tqdm import tqdm
from obspy import read, UTCDateTime, Stream, Trace

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



def pekeris_broadband(t0=3, t_max=50, fmin=1e1, fmax=1e2, dt=1e-3, zs=5, zr=5, r=1e4, d=20, num_mode=3, c1=1500, c2=2000 ):
    '''
    This function is based on the Pekeris class, and conputes the broadband signal in a Pekeris waveguide.
    
    parameters
    ----------
        t0 = time delay of the signal, default is 3 sec
        t_max = maximum simulation time, default is 50 sec 
        fmin = minimum frequency, default is 1e1 Hz 
        fmax = maximum frequency, default is 1e2 Hz 
        dt = time intevals, default is 1e-3 sec 
        zs = source depth, default is 5 m 
        zr = receiver depth, default is 5 m 
        r = receiver distance, default is 1e4 m
        d = waveguide depth, default is 20 m
        num_mode = number of modes to compute, default is 3
        c1 = velocity layer 1, default 1500 m/s
        c2 = velocity layer 2, default 2000 m/s
    '''

        
    l=int(t_max/dt)
    time = np.arange(0,dt*l, dt)

    t0_samp = int(t0/dt)

    sig_ref = np.zeros(time.size, dtype=np.float32)
    sig_ref[t0_samp] = 1

    tr_s = Trace()
    tr_s.stats.network = 'Signal'
    tr_s.stats.station = 'src'
    tr_s.stats.channel = 'FDH' # 
    tr_s.stats.starttime = UTCDateTime.now()
    tr_s.stats.sampling_rate = 1/dt
    tr_s.stats.delta = dt
    tr_s.stats.npts = time.size
    tr_s.stats.calib = 1
    tr_s.stats.units = 'Norm.'
    tr_s.data = sig_ref

    tr_s.filter('bandpass', freqmin=fmin, freqmax=fmax, corners=4)

    sig_ref = tr_s.data/np.max(tr_s.data)

    l = len(sig_ref)
    frequency = np.fft.fftfreq(l, dt)
    df = frequency[1]
    freq_ref = np.fft.fft(sig_ref)

    pressure = np.zeros(frequency.size, dtype=np.complex64)

    fmin_idx = (np.abs(frequency - fmin*0.5)).argmin()
    fmax_idx = (np.abs(frequency - fmax*1.5)).argmin()
    

    for idx, f in enumerate(tqdm(frequency[fmin_idx:fmax_idx])):
        P = pekeris(f=f, nq=5e4, dr=1, zs=1, d=20, c1=c1, c2=c2)
        P.calc_parameters()
        P.calc_field(r_rec=r, z_rec=zr, num_mode=num_mode)
        pressure[idx] = P.Phi*freq_ref[idx]
    pressure[idx] = P.Phi*freq_ref[idx]
    
        
    pressure_time = np.real(np.fft.fft(freq_ref*pressure*2))
    T = r/c1
    shift_samp = int(T/dt)
    signal = np.roll(pressure_time, shift=shift_samp)
    
    tr_r = tr_s.copy()
    tr_r.stats.station = 'rcv'
    tr_r.data = signal

    
    return tr_s, tr_r


