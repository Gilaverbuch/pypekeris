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
    def __init__(self,zs, zr, r, dr, c0,  c1,  c2, d, dz, f, nq):
        '''
        This function initializes the pekeris class


        parameters
        ----------
        zs = source depth. default 100 m
        zr = receiver depth, default 100m
        r  = receiver range, default 5000m
        dr = range interval, default 10 m
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
        self.dr = dr
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
        self._r = np.arange(self.dr, self.r, self.dr)

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
        plt.title('Modes shape for f=%d' %self.freq)
        plt.legend(bbox_to_anchor=(1, 1))
        plt.show()

# -------------------------------------------------------------------------------------------------


    def _calc_field(self):
        '''
        This function plots the shapes of the modes.
        '''

        Phi = np.zeros((self._z.size, self._r.size), dtype=np.complex64)
        l = len(self._idx)
        for idx in range(0, l,1):
            for z in range(0, self._z.size):
                H = np.sqrt(2/(np.pi * self.k0 * self._q[self._idx[idx]] * self._r)) \
                    * np.exp(1j*(self.k0 * self._q[self._idx[idx]] * self._r - np.pi/4))     

                phi_1 = (2*self.k0*self._mu2_abs[self._idx[idx]])
                phi_2 = (1 + self.k0*self.d*self._mu2_abs[self._idx[idx]])
                phi_3i = np.sin(self.k0*self._mu1[self._idx[idx]]*self._z[z])
                phi_i = np.sqrt(phi_1/phi_2)*phi_3i
                phi_3s = np.sin(self.k0*self._mu1[self._idx[idx]]*self.zs)
                phi_s = np.sqrt(phi_1/phi_2)*phi_3s

                Phi[z,:]+= (1j/4)*(H * phi_i * phi_s)
            
        self.Phi = Phi   


# -------------------------------------------------------------------------------------------------

    def _plot_TL(self):
        '''
        This function plots the shapes of the modes.
        '''

        plt.figure(figsize=(10, 5))
        plt.pcolormesh(self._r, self._z, 20*np.log10(np.abs(self.Phi)))
        plt.colorbar()
        plt.gca().invert_yaxis()
        plt.ylabel('Depth [m]')
        plt.xlabel('Range [m]')
        plt.title('TL for f=%d'%self.freq)
        plt.show()

# -------------------------------------------------------------------------------------------------

    def _plot_disperssion(self, fmin=1, fmax=50, df=0.5):
        '''
        This function plots the modal disperssion relation.
        '''

        frequency = np.arange(fmin, fmax, df)

        c_phase = []
        k_r = []
        f_k_r = []
        for f_idx, f in enumerate(frequency):
            self.freq = f
            self._calc_parameters()
            
            for modes in range(0, self._idx.size):
                c_phase.append(2*np.pi*f/(self.k0*self._q[self._idx[modes]]))
                k_r.append(self._q[self._idx[modes]])
                f_k_r.append(f)
                
        k_r = np.asarray(k_r, dtype=np.float32)
        f_k_R = np.asarray(f_k_r, dtype=np.float32)
        c_phase = np.asarray(c_phase, dtype=np.float32)

        plt.figure(figsize=(5,5))
        plt.plot(c_phase, f_k_R, 'k.')
        plt.title('Modal disperssion relation')
        plt.ylabel('Frequency [Hz]')
        plt.xlabel('Phase velocity [m/s]')
        plt.show()


