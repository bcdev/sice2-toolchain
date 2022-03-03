# -*- coding: utf-8 -*-
"""
TODO
@author: bav@geus.dk, Olaf Danne (BC)
"""

import numpy as np
from sice2_constants import w, bai, xa, ya, f0, f1, f2, bet, gam, coef1, coef2, coef3, coef4


def funp(x, al, sph_calc, ak1):
    #     Spectral planar albedo
    # Inputs:
    # x                     input wavelength (should work with any)
    # ak1
    # al                    absorption length
    # sph_calc              sph_calc= 0 for planar =1 for spherical
    #
    # Constants:
    # xa(168),ya(168)       imaginary part (ya) of the refraction index at specified wavelength (xa)
    #
    # Outputs:
    # f1*funcs              ?
    #
    # bav 2020
    # using numpy interpolation

    y = np.interp(x, xa, ya)

    dega = 1000. * al * 4. * np.pi * y / x
    pow = np.sqrt(dega)
         
    if (pow >= 1.e-6): 
        rsd = np.exp(-pow)
    else: 
        rsd = 1.
         
    if (sph_calc == 0):     
        rs = rsd**ak1
    elif (sph_calc == 1):   
        rs = rsd

    if (x < 0.4):  
        x = 0.4
    funcs = f0 + f1 * np.exp(-x * bet) + f2 * np.exp(-x * gam)
     
    return rs * funcs


def qsimp(func, a, b):
    # integrate function between a and b using simpson's method. 
    # works as fast as scipy.integrate quad
    eps = 1.e-3
    jmax = 20
    ost = -1.e30
    os = -1.e30
    
    for j in range(jmax):
        
        if (j == 0):
            st = 0.5 * (b - a) * (func(a) + func(b))
        else:
            it = 2 ** (j - 1)
            tnm = it
            delta = (b - a) / tnm
            x = a + 0.5 * delta
            sum = 0.
            
            for jj in range(it):
                
                sum = sum + func(x)
                x = x + delta
            st = 0.5 * (st + (b - a) * sum / tnm)
        s = (4. * st - ost) / 3.
        if (j > 4):
            if (abs(s - os) < eps * abs(os)):
                return s
            if (s == 0) and (os == 0.):
                return s
        os = s
        ost = st
    print("Max iteration reached")
    
    return s
