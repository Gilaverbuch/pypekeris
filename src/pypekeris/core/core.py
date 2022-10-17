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
    def __init__(self,zs, zr, r, c0,  c1,  c2, d, dz, f, nq):
        '''
        This function initializes the pekeris class


        parameters
        ----------
        zs = source depth. default 100 m
        zr = receiver depth, default 100m
        r  = receiver range, default 5000m
        c0 = reference velocity, default 1500 m/s
        c1 = velocity layer 1, default 1500 m/s
        c2 = velocity layer 2, default 2000 m/s
        d  = depth of layer 1, default 1000 m
        f  = frequency, default is 1 Hz
        nq = number of wavenumbers, default is 1e3

        Returns
        -------
        Header: header structure containing the parameters names and values (Pandas DataFrame)
        Waveforms: Obspy Stream object
        '''

        # parameters = {
        #             'zs' : [zs],
        #             'zr' : [zr],
        #             'r'  : [r],
        #             'c0' : [c0],
        #             'c1' : [c1],
        #             'c2' : [c2],
        #             'd'  : [d],
        #             'freq' : [f],
        #             'nq' : [nq]
        #                    }

        # self.parameters = pd.DataFrame(parameters)

        self.zs = zs
        self.zr = zr
        self.r = r
        self.c0 = c0
        self.c1 = c1
        self.c2 = c2
        self.d  = d 
        self.dz = 1
        self.freq = f
        self.nq = nq




# -------------------------------------------------------------------------------------------------
    def _calc_parameters(self):
        '''
        This function calculates the parameters needed for the computations
        '''
        self.omega = 2*np.pi*self.freq #angular frequency
        self.k0 = self.omega/self.c0 # reference wavenumber
        self.k1 = self.c0/self.c1 # dimentionless wavenumber layer 1
        self.k2 = self.c0/self.c2 # dimentionless wavenumber layer 2
        self.dq = (self.k1-self.k2)/self.nq # horizontal dimentionless wavenumber resolution 
        self._q = np.arange(self.k2, self.k1, self.dq) # horizontal dimentionless wavenumber
        self._mu1 = np.sqrt(self.k1**2 - self._q**2) # dimentionless vertical wavenumber layer 1
        self._mu2_abs = np.sqrt(self._q**2 - self.k2**2) # dimentionless vertical wavenumber layer 2
        self._mu2 = 1j*self._mu2_abs
        self._mu_func = -self._mu2_abs/self._mu1 
        self._cot_func = 1/np.tan(self.k0*self._mu1*self.d)

        self._idx = np.argwhere(np.diff(np.sign(self._mu_func - self._cot_func))).flatten()
        self._idx = self._idx[self._cot_func[self._idx]<=10]

        self._z = np.arange(0, self.d, self.dz)

# -------------------------------------------------------------------------------------------------

    def _print_parameters(self):
        '''
        This function prints the single parameters used for the calculations
        '''

        for attr, value in self.__dict__.items():
            if not attr.startswith('_'):
                    print(attr, value)

    
# -------------------------------------------------------------------------------------------------

    def _plot_discrete_modes(self):
        '''
        This function plots the _cot_func, _mu_func, and their intersection points
        '''

        plt.figure(figsize = (10,5))
        plt.plot(self.k1, 0, 'k.', markersize=8)
        plt.plot(self.k2, 0, 'k.', markersize=8)
        plt.vlines(self.k1, -10, 10 , colors='k', linestyles='dashed')
        plt.vlines(self.k2, -10, 10 , colors='k', linestyles='dashed')
        plt.annotate('k1', [self.k1, 0])
        plt.annotate('k2', [self.k2, 0])
        plt.plot(self._q, self._mu_func, 'k')
        plt.plot(self._q, self._cot_func, 'k.')
        plt.plot(self._q[self._idx], self._cot_func[self._idx], 'ro')
        plt.ylim(-10,10)
        plt.title('Intersections of $\cot (k_0 \mu_1 d)$ and $-|\mu_2|/\mu_1$. Discrete eigenvalues')
        plt.xlabel('Non-dimensional horizontal wavenumber')
        plt.ylabel('Amplitude')
        plt.show()

# -------------------------------------------------------------------------------------------------


    def _plot_modes_shape(self):
        '''
        This function plots the shapes of the modes.
        '''

        plt.figure(figsize=(4,6))
        l = len(self._idx)
        for idx in range(0, l,1): 
            phi_1 = (2*self.k0*self._mu2_abs[self._idx[idx]])
            phi_2 = (1 + self.k0*self.d*self._mu2_abs[self._idx[idx]])
            phi_3 = np.sin(self.k0*self._mu1[self._idx[idx]]*self._z)
            phi_i = np.sqrt(phi_1/phi_2)*phi_3
            

            plt.plot(phi_i, self._z, label='mode ' + str(l-idx))

        plt.gca().invert_yaxis()
        plt.ylabel('Depth [m]')
        plt.xlabel('Amplitude')
        plt.title('Modes shape')
        plt.legend(bbox_to_anchor=(1, 1))
        plt.show()

# -------------------------------------------------------------------------------------------------

