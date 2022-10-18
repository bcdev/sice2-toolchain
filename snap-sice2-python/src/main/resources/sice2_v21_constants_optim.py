# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 16:58:31 2019
Update 07032019
This script contains the constants needed by the pySICEv1.1 library.
@author: bav@geus.dk
"""

import numpy as np
import xarray as xr


bandcoord = ('band', np.arange(21))


# solar spectrum constants
f0 = 32.38
f1 = -160140.33
f2 = 7959.53
bet = 1./(85.34*1.e-3)
gam = 1./(401.79*1.e-3)

# ICE REFRATIVE INDEX
xa = np.array([2.010E-001, 2.019E-001, 2.100E-001, 2.500E-001, 3.00E-001, 3.500E-001, 3.900E-001, 4.000E-001,
               4.100E-001, 4.200E-001, 4.300E-001, 4.400E-001, 4.500E-001, 4.600E-001, 4.700E-001, 4.800E-001, 4.900E-001, 5.000E-001,
               5.100E-001, 5.200E-001, 5.300E-001, 5.400E-001, 5.500E-001, 5.600E-001, 5.700E-001, 5.800E-001, 5.900E-001, 6.000E-001,
               6.100E-001, 6.200E-001, 6.300E-001, 6.400E-001, 6.500E-001, 6.600E-001, 6.700E-001, 6.800E-001, 6.900E-001, 7.000E-001,
               7.100E-001, 7.200E-001, 7.300E-001, 7.400E-001, 7.500E-001, 7.600E-001, 7.700E-001, 7.800E-001, 7.900E-001, 8.000E-001,
               8.100E-001, 8.200E-001, 8.300E-001, 8.400E-001, 8.500E-001, 8.600E-001, 8.700E-001, 8.800E-001, 8.900E-001, 9.000E-001,
               9.100E-001, 9.200E-001, 9.300E-001, 9.400E-001, 9.500E-001, 9.600E-001, 9.700E-001, 9.800E-001, 9.900E-001, 1.000E+000,
               1.010E+000, 1.020E+000, 1.030E+000, 1.040E+000, 1.050E+000, 1.060E+000, 1.070E+000, 1.080E+000, 1.090E+000, 1.100E+000,
               1.110E+000, 1.120E+000, 1.130E+000, 1.140E+000, 1.150E+000, 1.160E+000, 1.170E+000, 1.180E+000, 1.190E+000, 1.200E+000,
               1.210E+000, 1.220E+000, 1.230E+000, 1.240E+000, 1.250E+000, 1.260E+000, 1.270E+000, 1.280E+000, 1.290E+000, 1.300E+000,
               1.310E+000, 1.320E+000, 1.330E+000, 1.340E+000, 1.350E+000, 1.360E+000, 1.370E+000, 1.380E+000, 1.390E+000, 1.400E+000,
               1.410E+000, 1.420E+000, 1.430E+000, 1.440E+000, 1.449E+000, 1.460E+000, 1.471E+000, 1.481E+000, 1.493E+000, 1.504E+000,
               1.515E+000, 1.527E+000, 1.538E+000, 1.563E+000, 1.587E+000, 1.613E+000, 1.650E+000, 1.680E+000, 1.700E+000, 1.730E+000,
               1.760E+000, 1.800E+000, 1.830E+000, 1.840E+000, 1.850E+000, 1.855E+000, 1.860E+000, 1.870E+000, 1.890E+000, 1.905E+000,
               1.923E+000, 1.942E+000, 1.961E+000, 1.980E+000, 2.000E+000, 2.020E+000, 2.041E+000, 2.062E+000, 2.083E+000, 2.105E+000,
               2.130E+000, 2.150E+000, 2.170E+000, 2.190E+000, 2.220E+000, 2.240E+000, 2.245E+000, 2.250E+000, 2.260E+000, 2.270E+000,
               2.290E+000, 2.310E+000, 2.330E+000, 2.350E+000, 2.370E+000, 2.390E+000, 2.410E+000, 2.430E+000, 2.460E+000, 2.500E+000])

ya = np.array([3.249E-011, 2.0E-011, 2.0E-011, 2.0E-011, 2.0E-011, 2.0E-011, 2.0E-011, 2.365E-011, 2.669E-011, 3.135E-011,
               4.140E-011, 6.268E-011, 9.239E-011, 1.325E-010, 1.956E-010, 2.861E-010, 4.172E-010, 5.889E-010, 8.036E-010, 1.076E-009,
               1.409E-009, 1.813E-009, 2.289E-009, 2.839E-009, 3.461E-009, 4.159E-009, 4.930E-009, 5.730E-009, 6.890E-009, 8.580E-009,
               1.040E-008, 1.220E-008, 1.430E-008, 1.660E-008, 1.890E-008, 2.090E-008, 2.400E-008, 2.900E-008, 3.440E-008, 4.030E-008,
               4.300E-008, 4.920E-008, 5.870E-008, 7.080E-008, 8.580E-008, 1.020E-007, 1.180E-007, 1.340E-007, 1.400E-007, 1.430E-007,
               1.450E-007, 1.510E-007, 1.830E-007, 2.150E-007, 2.650E-007, 3.350E-007, 3.920E-007, 4.200E-007, 4.440E-007, 4.740E-007,
               5.110E-007, 5.530E-007, 6.020E-007, 7.550E-007, 9.260E-007, 1.120E-006, 1.330E-006, 1.620E-006, 2.000E-006, 2.250E-006,
               2.330E-006, 2.330E-006, 2.170E-006, 1.960E-006, 1.810E-006, 1.740E-006, 1.730E-006, 1.700E-006, 1.760E-006, 1.820E-006,
               2.040E-006, 2.250E-006, 2.290E-006, 3.040E-006, 3.840E-006, 4.770E-006, 5.760E-006, 6.710E-006, 8.660E-006, 1.020E-005,
               1.130E-005, 1.220E-005, 1.290E-005, 1.320E-005, 1.350E-005, 1.330E-005, 1.320E-005, 1.320E-005, 1.310E-005, 1.320E-005,
               1.320E-005, 1.340E-005, 1.390E-005, 1.420E-005, 1.480E-005, 1.580E-005, 1.740E-005, 1.980E-005, 3.442E-005, 5.959E-005,
               1.028E-004, 1.516E-004, 2.030E-004, 2.942E-004, 3.987E-004, 4.941E-004, 5.532E-004, 5.373E-004, 5.143E-004, 4.908E-004,
               4.594E-004, 3.858E-004, 3.105E-004, 2.659E-004, 2.361E-004, 2.046E-004, 1.875E-004, 1.650E-004, 1.522E-004, 1.411E-004,
               1.302E-004, 1.310E-004, 1.339E-004, 1.377E-004, 1.432E-004, 1.632E-004, 2.566E-004, 4.081E-004, 7.060E-004, 1.108E-003,
               1.442E-003, 1.614E-003, 1.640E-003, 1.566E-003, 1.458E-003, 1.267E-003, 1.023E-003, 7.586E-004, 5.255E-004, 4.025E-004,
               3.235E-004, 2.707E-004, 2.228E-004, 2.037E-004, 2.026E-004, 2.035E-004, 2.078E-004, 2.171E-004, 2.538E-004, 3.138E-004,
               3.858E-004, 4.591E-004, 5.187E-004, 5.605E-004, 5.956E-004, 6.259E-004, 6.820E-004, 7.530E-004])


# OLCI channels
wls = xr.DataArray([0.4000E+00, 0.4125E+00, 0.4425E+00, 0.4900E+00, 0.5100E+00, 0.5600E+00, 0.6200E+00,
                    0.6650E+00, 0.6737E+00,  0.6812E+00, 0.7088E+00, 0.7538E+00, 0.7613E+00, 0.7644E+00, 0.7675E+00, 0.7788E+00,
                    0.8650E+00, 0.8850E+00, 0.9000E+00, 0.9400E+00, 0.1020E+01], coords=[bandcoord])

# Imaginary part of ice refrative index at OLCI channels
bai = xr.DataArray([2.365E-11, 2.7E-11, 7.0E-11, 4.17E-10,
                    8.04E-10,  2.84E-09, 8.58E-09,  1.78E-08,  1.95E-08, 2.1E-08, 3.3E-08, 6.23E-08, 7.1E-08,  7.68E-08,  8.13E-08,
                    9.88E-08,  2.4E-07, 3.64E-07,  4.2E-07, 5.53e-07, 2.25E-06], coords=[bandcoord])

# %% Solar flux


def sol(x):
    # SOLAR SPECTRUM at GROUND level
    # Inputs:
    # x         wave length in micrometer
    # Outputs:
    # sol       solar spectrum in W m-2 micrometer-1 (?)
    #    if (x < 0.4):
    #            x=0.4
    sol1a = f0*x
    sol1b = - f1*np.exp(-bet*x)/bet
    sol1c = - f2*np.exp(-gam*x)/gam
    return sol1a+sol1b+sol1c


sol0 = (f0 + f1*np.exp(-bet * 0.4) + f2*np.exp(-gam * 0.4))*0.1

# solar flux calculation
# sol1      visible(0.3-0.7micron)
# somehow, a different sol1 needs to be used for clean snow and polluted snow
sol1_clean = sol(0.7) - sol(0.4) + sol0
sol1_pol = sol(0.7) - sol(0.3)
# sol2      near-infrared (0.7-2.4micron)
# same for clean and polluted
sol2 = sol(2.4) - sol(0.7)

# sol3      shortwave(0.3-2.4 micron)
# sol3 is also different for clean snow and polluted snow
sol3_clean = sol1_clean + sol2
sol3_pol = sol1_pol + sol2

# asol specific band
asol = sol(0.865) - sol(0.7)

# %% analystical integration of the solar flux


def analyt_func(z1, z2):
    # see BBA_calc_pol
    # compatible with array
    ak1 = (z2**2.-z1**2.)/2.
    ak2 = (z2/bet+1./bet/bet)*np.exp(-bet*z2) - (z1/bet+1./bet/bet)*np.exp(-bet*z1)
    ak3 = (z2/gam+1./gam**2)*np.exp(-gam*z2) - (z1/gam+1./gam**2)*np.exp(-gam*z1)

    am1 = (z2**3.-z1**3.)/3.
    am2 = (z2**2./bet+2.*z2/bet**2 + 2./bet**3) * np.exp(-bet*z2) \
          - (z1**2./bet+2.*z1/bet**2 + 2./bet**3) * np.exp(-bet*z1)
    am3 = (z2**2./gam+2.*z2/gam**2 + 2./gam**3.)*np.exp(-gam*z2) \
          - (z1**2./gam+2.*z1/gam**2 + 2./gam**3.)*np.exp(-gam*z1)

    return (f0*ak1 - f1*ak2 - f2*ak3), (f0*am1 - f1*am2 - f2*am3)


# %% solar constant
coef1, coef2 = analyt_func(0.3, 0.7)
coef3, coef4 = analyt_func(0.7, 0.865)