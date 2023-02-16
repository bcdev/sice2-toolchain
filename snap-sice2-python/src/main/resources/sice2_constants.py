# -*- coding: utf-8 -*-
"""
This script contains the constants needed by the SICE2 processors.

@author: Baptiste Vandecrux, GEUS (Geological Survey of Denmark and Greenland)
@author: Olaf Danne, BC (Brockmann Consult): SNAP integration.
"""

import os

import numpy as np
import xarray as xr

bandcoord = ("band", np.arange(21))

SNOW_TYPE_FLAG = 'SNOW_TYPE_FLAG'
POL_TYPE_FLAG = 'POL_TYPE'
SCDA_FLAG = 'SCDA_FLAG'

# OLCI channels
wls = xr.DataArray(
    [
        0.4000e00,
        0.4125e00,
        0.4425e00,
        0.4900e00,
        0.5100e00,
        0.5600e00,
        0.6200e00,
        0.6650e00,
        0.6737e00,
        0.6812e00,
        0.7088e00,
        0.7538e00,
        0.7613e00,
        0.7644e00,
        0.7675e00,
        0.7788e00,
        0.8650e00,
        0.8850e00,
        0.9000e00,
        0.9400e00,
        0.1020e01,
    ],
    coords=[bandcoord],
)

# Imaginary part of ice refrative index at OLCI channels
bai = xr.DataArray(
    [
        6.27e-10,
        5.78e-10,
        6.49e-10,
        1.08e-9,
        1.46e-9,
        3.35e-09,
        8.58e-09,
        1.78e-08,
        1.95e-08,
        2.1e-08,
        3.3e-08,
        6.23e-08,
        7.1e-08,
        7.68e-08,
        8.13e-08,
        9.88e-08,
        2.4e-07,
        3.64e-07,
        4.2e-07,
        5.53e-07,
        2.25e-06,
    ],
    coords=[bandcoord],
)

# BULK ice absorption coefficient at all OLCI channels  (1/nm)
alpha = 4 * np.pi * bai / wls

# ozone vertical optical density  at OLCI channels
cabsoz = xr.DataArray(
    [
        1.3782e-4,
        3.0488e-4,
        1.6457e-3,
        8.9359e-3,
        1.7505e-2,
        4.3471e-2,
        4.4871e-2,
        2.1016e-2,
        1.7162e-2,
        1.4663e-2,
        7.9830e-3,
        3.8797e-3,
        2.9338e-3,
        2.7992e-3,
        2.7297e-3,
        3.2560e-3,
        8.9569e-4,
        5.1888e-4,
        6.7158e-4,
        3.1278e-4,
        1.4088e-5,
    ],
    coords=[bandcoord],
)

# relative vertical optical density of ozone f at all OLCI channels
# (normalized to that at 620nm)
# f = cabsoz/cabsoz.isel(band=6)

# new constants defined by Alex
# aerosol properties(aot=aot500nm,ANNA=Angström exponent)
aot = 0.07
anna = 1.3

# patchy snow threshold for the channel R(400nm)
thv0 = 0.75
# threshold for the measured and simulated spectra differences
thv1 = 5.0
# threshold for the grain diameter
thv2 = 0.07
# threshold for the total ozone column retrieved
thv3 = 1000.0

# %% Solar flux
# solar spectrum constants
f0 = 32.38
f1 = -160140.33
f2 = 7959.53
bet = 11.72  # 1.0 / 0.08534
gam = 2.49  # 1.0 / 0.40179


def sol(x):
    # antiderivative of the SOLAR SPECTRUM at GROUND level
    # Inputs:
    # x         wave length in micrometer
    # Outputs:
    # sol       solar spectrum in W m-2 micrometer-1 (?)
    # if (x < 0.4):
    #         x=0.4
    sol1a = f0 * x
    sol1b = -f1 * np.exp(-bet * x) / bet
    sol1c = -f2 * np.exp(-gam * x) / gam
    return sol1a + sol1b + sol1c


# solar flux calculation
# sol1      visible(0.3-0.7micron)
# Update 2022: same for clean and polluted
sol_vis = sol(0.7) - sol(0.33)
# sol2      near-infrared (0.7-2.4micron)
# Update 2022: same for clean and polluted
sol_nir = sol(2.4) - sol(0.7)
# sol3      shortwave(0.3-2.4 micron)
sol_sw = sol_vis + sol_nir

# asol specific band
asol = sol(0.865) - sol(0.7)

# %% ICE REFRATIVE INDEX
xa = np.array(
    [
        0.32,
        0.325,
        0.33,
        0.335,
        0.34,
        0.345,
        0.35,
        0.355,
        0.36,
        0.365,
        0.37,
        0.375,
        0.38,
        0.385,
        0.39,
        0.395,
        0.4,
        0.405,
        0.41,
        0.415,
        0.42,
        0.425,
        0.43,
        0.435,
        0.44,
        0.445,
        0.45,
        0.455,
        0.46,
        0.465,
        0.47,
        0.475,
        0.48,
        0.485,
        0.49,
        0.495,
        0.5,
        0.505,
        0.51,
        0.515,
        0.52,
        0.525,
        0.53,
        0.535,
        0.54,
        0.545,
        0.55,
        0.555,
        0.56,
        0.565,
        0.57,
        0.575,
        0.58,
        0.585,
        0.59,
        0.595,
        0.6,
        0.61,
        0.62,
        0.63,
        0.64,
        0.65,
        0.66,
        0.67,
        0.68,
        0.69,
        0.7,
        0.71,
        0.72,
        0.73,
        0.74,
        0.75,
        0.76,
        0.77,
        0.78,
        0.79,
        0.8,
        0.81,
        0.82,
        0.83,
        0.84,
        0.85,
        0.86,
        0.87,
        0.88,
        0.89,
        0.9,
        0.91,
        0.92,
        0.93,
        0.94,
        0.95,
        0.96,
        0.97,
        0.98,
        0.99,
        1.0,
        1.01,
        1.02,
        1.03,
        1.04,
        1.05,
        1.06,
        1.07,
        1.08,
        1.09,
        1.1,
        1.11,
        1.12,
        1.13,
        1.14,
        1.15,
        1.16,
        1.17,
        1.18,
        1.19,
        1.2,
        1.21,
        1.22,
        1.23,
        1.24,
        1.25,
        1.26,
        1.27,
        1.28,
        1.29,
        1.3,
        1.31,
        1.32,
        1.33,
        1.34,
        1.35,
        1.36,
        1.37,
        1.38,
        1.39,
        1.4,
        1.41,
        1.42,
        1.43,
        1.44,
        1.449,
        1.46,
        1.471,
        1.481,
        1.493,
        1.504,
        1.515,
        1.527,
        1.538,
        1.563,
        1.587,
        1.613,
        1.65,
        1.68,
        1.7,
        1.73,
        1.76,
        1.8,
        1.83,
        1.84,
        1.85,
        1.855,
        1.86,
        1.87,
        1.89,
        1.905,
        1.923,
        1.942,
        1.961,
        1.98,
        2.0,
        2.02,
        2.041,
        2.062,
        2.083,
        2.105,
        2.13,
        2.15,
        2.17,
        2.19,
        2.22,
        2.24,
        2.245,
        2.25,
        2.26,
        2.27,
        2.29,
        2.31,
        2.33,
        2.35,
        2.37,
        2.39,
        2.41,
        2.43,
        2.46,
        2.5,
    ]
)

ya = np.array(
    [
        8.38875029e-10,
        7.95734990e-10,
        7.90943827e-10,
        7.71359229e-10,
        7.47344657e-10,
        7.14322838e-10,
        6.99103321e-10,
        6.78074139e-10,
        6.63775953e-10,
        6.55295670e-10,
        6.47804333e-10,
        6.39936528e-10,
        6.42039090e-10,
        6.27801214e-10,
        6.33959963e-10,
        6.16145700e-10,
        6.27511173e-10,
        6.06103626e-10,
        5.84866807e-10,
        5.71551780e-10,
        5.76446512e-10,
        5.86493711e-10,
        5.97242486e-10,
        6.19603170e-10,
        6.37426183e-10,
        6.61162140e-10,
        6.89556417e-10,
        7.18627374e-10,
        7.51969388e-10,
        7.89963283e-10,
        8.40212398e-10,
        8.92133312e-10,
        9.46872015e-10,
        1.00892822e-09,
        1.07911676e-09,
        1.15933631e-09,
        1.24566244e-09,
        1.34463623e-09,
        1.45679374e-09,
        1.56520969e-09,
        1.69596619e-09,
        1.83341271e-09,
        1.98568407e-09,
        2.15727694e-09,
        2.34704124e-09,
        2.56573477e-09,
        2.80077577e-09,
        3.06346813e-09,
        3.34610379e-09,
        3.64478242e-09,
        3.95975067e-09,
        4.30195534e-09,
        4.65194061e-09,
        5.04250741e-09,
        5.47611863e-09,
        5.96820618e-09,
        5.73000000e-09,
        6.89000000e-09,
        8.58000000e-09,
        1.04000000e-08,
        1.22000000e-08,
        1.43000000e-08,
        1.66000000e-08,
        1.89000000e-08,
        2.09000000e-08,
        2.40000000e-08,
        2.90000000e-08,
        3.44000000e-08,
        4.03000000e-08,
        4.30000000e-08,
        4.92000000e-08,
        5.87000000e-08,
        7.08000000e-08,
        8.58000000e-08,
        1.02000000e-07,
        1.18000000e-07,
        1.34000000e-07,
        1.40000000e-07,
        1.43000000e-07,
        1.45000000e-07,
        1.51000000e-07,
        1.83000000e-07,
        2.15000000e-07,
        2.65000000e-07,
        3.35000000e-07,
        3.92000000e-07,
        4.20000000e-07,
        4.44000000e-07,
        4.74000000e-07,
        5.11000000e-07,
        5.53000000e-07,
        6.02000000e-07,
        7.55000000e-07,
        9.26000000e-07,
        1.12000000e-06,
        1.33000000e-06,
        1.62000000e-06,
        2.00000000e-06,
        2.25000000e-06,
        2.33000000e-06,
        2.33000000e-06,
        2.17000000e-06,
        1.96000000e-06,
        1.81000000e-06,
        1.74000000e-06,
        1.73000000e-06,
        1.70000000e-06,
        1.76000000e-06,
        1.82000000e-06,
        2.04000000e-06,
        2.25000000e-06,
        2.29000000e-06,
        3.04000000e-06,
        3.84000000e-06,
        4.77000000e-06,
        5.76000000e-06,
        6.71000000e-06,
        8.66000000e-06,
        1.02000000e-05,
        1.13000000e-05,
        1.22000000e-05,
        1.29000000e-05,
        1.32000000e-05,
        1.35000000e-05,
        1.33000000e-05,
        1.32000000e-05,
        1.32000000e-05,
        1.31000000e-05,
        1.32000000e-05,
        1.32000000e-05,
        1.34000000e-05,
        1.39000000e-05,
        1.42000000e-05,
        1.48000000e-05,
        1.58000000e-05,
        1.74000000e-05,
        1.98000000e-05,
        3.44200000e-05,
        5.95900000e-05,
        1.02800000e-04,
        1.51600000e-04,
        2.03000000e-04,
        2.94200000e-04,
        3.98700000e-04,
        4.94100000e-04,
        5.53200000e-04,
        5.37300000e-04,
        5.14300000e-04,
        4.90800000e-04,
        4.59400000e-04,
        3.85800000e-04,
        3.10500000e-04,
        2.65900000e-04,
        2.36100000e-04,
        2.04600000e-04,
        1.87500000e-04,
        1.65000000e-04,
        1.52200000e-04,
        1.41100000e-04,
        1.30200000e-04,
        1.31000000e-04,
        1.33900000e-04,
        1.37700000e-04,
        1.43200000e-04,
        1.63200000e-04,
        2.56600000e-04,
        4.08100000e-04,
        7.06000000e-04,
        1.10800000e-03,
        1.44200000e-03,
        1.61400000e-03,
        1.64000000e-03,
        1.56600000e-03,
        1.45800000e-03,
        1.26700000e-03,
        1.02300000e-03,
        7.58600000e-04,
        5.25500000e-04,
        4.02500000e-04,
        3.23500000e-04,
        2.70700000e-04,
        2.22800000e-04,
        2.03700000e-04,
        2.02600000e-04,
        2.03500000e-04,
        2.07800000e-04,
        2.17100000e-04,
        2.53800000e-04,
        3.13800000e-04,
        3.85800000e-04,
        4.59100000e-04,
        5.18700000e-04,
        5.60500000e-04,
        5.95600000e-04,
        6.25900000e-04,
        6.82000000e-04,
        7.53000000e-04,
    ]
)

SNAPPY_DIR = os.path.expanduser('~') + os.sep + '.snap' + os.sep + 'snap-python'

OLCI_NUM_SPECTRAL_BANDS = 21

AOT = 0.1

# OLCI channels
w = np.array((
    0.4000E+00, 0.4125E+00, 0.4425E+00, 0.4900E+00, 0.5100E+00, 0.5600E+00, 0.6200E+00, 0.6650E+00, 0.6737E+00,
    0.6812E+00, 0.7088E+00, 0.7538E+00, 0.7613E+00, 0.7644E+00, 0.7675E+00, 0.7788E+00, 0.8650E+00, 0.8850E+00,
    0.9000E+00, 0.9400E+00, 0.1020E+01)).astype('float32')

# Imaginary part of ice refrative index at OLCI channels
bai = np.array((2.365E-11, 2.7E-11, 7.0E-11, 4.17E-10, 8.04E-10, 2.84E-09, 8.58E-09, 1.78E-08, 1.95E-08, 2.1E-08,
                3.3E-08, 6.23E-08, 7.1E-08, 7.68E-08, 8.13E-08, 9.88E-08, 2.4E-07, 3.64E-07, 4.2E-07, 5.53e-07,
                2.25E-06)).astype('float32')

tg_vod = np.array((
    1.378170469E-004,
    3.048780958E-004,
    1.645714060E-003,
    8.935947110E-003,
    1.750535146E-002,
    4.347104369E-002,
    4.487130794E-002,
    2.101591797E-002,
    1.716230955E-002,
    1.466298300E-002,
    7.983028470E-003,
    3.879744653E-003,
    2.923775641E-003,
    2.792211429E-003,
    2.729651478E-003,
    3.255969698E-003,
    8.956858078E-004,
    5.188799343E-004,
    6.715773241E-004,
    3.127781417E-004,
    1.408798425E-005
)).astype('float64')

tg_water_vod = np.array((
    3.506319806E-005,
    6.967526062E-006,
    3.505636418E-004,
    1.381931600E-004,
    1.324224122E-003,
    5.099066205E-005,
    2.363042944E-005,
    6.440468944E-004,
    2.377942129E-005,
    3.303735591E-004,
    1.158079797E-002,
    1.154852173E-004,
    3.316322348E-006,
    2.465775843E-007,
    1.137753056E-007,
    3.191069529E-004,
    5.601462094E-004,
    3.026113840E-003,
    3.265557172E-001,
    2.804599956E+000,
    2.858762309E-003
)).astype('float64')


# %% Solar flux
def sol(x):
    # SOLAR SPECTRUM at GROUND level
    # Inputs:
    # x         wave length in micrometer
    # Outputs:
    # sol       solar spectrum in W m-2 micrometer-1 (?)
    #    if (x < 0.4):
    #            x=0.4
    sol1a = f0 * x
    sol1b = - f1 * np.exp(-bet * x) / bet
    sol1c = - f2 * np.exp(-gam * x) / gam
    return sol1a + sol1b + sol1c


sol0 = (f0 + f1 * np.exp(-bet * 0.4) + f2 * np.exp(-gam * 0.4)) * 0.1

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


# %% analystical function and coefficients used in the polluted snow BBA calculation
def analyt_func(z1, z2):
    # see BBA_calc_pol
    # compatible with array
    ak1 = (z2 ** 2. - z1 ** 2.) / 2.
    ak2 = (z2 / bet + 1. / bet / bet) * np.exp(-bet * z2) - (z1 / bet + 1. / bet / bet) * np.exp(-bet * z1)
    ak3 = (z2 / gam + 1. / gam ** 2) * np.exp(-gam * z2) - (z1 / gam + 1. / gam ** 2) * np.exp(-gam * z1)

    am1 = (z2 ** 3. - z1 ** 3.) / 3.
    am2 = (z2 ** 2. / bet + 2. * z2 / bet ** 2 + 2. / bet ** 3) * np.exp(-bet * z2) \
        - (z1 ** 2. / bet + 2. * z1 / bet ** 2 + 2. / bet ** 3) * np.exp(-bet * z1)
    am3 = (z2 ** 2. / gam + 2. * z2 / gam ** 2 + 2. / gam ** 3.) * np.exp(-gam * z2) \
        - (z1 ** 2. / gam + 2. * z1 / gam ** 2 + 2. / gam ** 3.) * np.exp(-gam * z1)

    return (f0 * ak1 - f1 * ak2 - f2 * ak3), (f0 * am1 - f1 * am2 - f2 * am3)


coef1, coef2 = analyt_func(0.3, 0.7)
coef3, coef4 = analyt_func(0.7, 0.865)

mandatory_tif_inputs_for_retrieval = ['r_TOA_01.tif', 'r_TOA_04.tif', 'r_TOA_17.tif', 'r_TOA_21.tif',
                                      'SZA.tif', 'SAA.tif', 'OZA.tif', 'OAA.tif',
                                      'O3.tif', 'height.tif']

pol_type_flags_map = {
    str(0.0): 1,
    str(1.0): 2,
    str(2.0): 4,
    str(3.0): 8
}

snow_type_flags_map = {
    str(1.0): 1,
    str(2.0): 2,
    str(3.0): 4,
    str(100.0): 8,
    str(102.0): 16,
    str(103.0): 32,
    str(104.0): 64,
    str(105.0): 128
}

OLCI_L2_IPF_BITMASK_FLAG_CODING_DICT = {
    'WQSF.INVALID': 1,
    'WQSF.WATER': 2,
    'WQSF.LAND': 4,
    'WQSF.CLOUD': 8,
    'WQSF.CLOUD_AMBIGUOUS': 8388608,
    'WQSF.CLOUD_MARGIN': 16777216,
    'WQSF.SNOW_ICE': 16,
    'WQSF.INLAND_WATER': 32,
    'WQSF.TIDAL': 64,
    'WQSF.COSMETIC': 128,
    'WQSF.SUSPECT': 256,
    'WQSF.HISOLZEN': 512,
    'WQSF.SATURATED': 1024,
    'WQSF.MEGLINT': 2048,
    'WQSF.HIGHGLINT': 4096,
    'WQSF.WHITECAPS': 8192,
    'WQSF.ADJAC': 16384,
    'WQSF.WV_FAIL': 32768,
    'WQSF.PAR_FAIL': 65536,
    'WQSF.AC_FAIL': 131072,
    'WQSF.OC4ME_FAIL': 262144,
    'WQSF.OCNN_FAIL': 524288,
    'WQSF.KDM_FAIL': 2097152,
    'WQSF.BPAC_ON': 33554432,
    'WQSF.WHITE_SCATT': 67108864,
    'WQSF.LOWRW': 134217728,
    'WQSF.HIGHRW': 268435456,
    'WQSF.ANNOT_ANGSTROM': 4294967296,
    'WQSF.ANNOT_AERO_B': 8589934592,
    'WQSF.ANNOT_ABSO_D': 17179869184,
    'WQSF.ANNOT_ACLIM': 34359738368,
    'WQSF.ANNOT_ABSO_A': 68719476736,
    'WQSF.ANNOT_MIXR1': 137438953472,
    'WQSF.ANNOT_DROUT': 274877906944,
    'WQSF.ANNOT_TAU06': 549755813888,
    'WQSF.RWNEG_O1': 1099511627776,
    'WQSF.RWNEG_O2': 2199023255552,
    'WQSF.RWNEG_O3': 4398046511104,
    'WQSF.RWNEG_O4': 8796093022208,
    'WQSF.RWNEG_O5': 17592186044416,
    'WQSF.RWNEG_O6': 35184372088832,
    'WQSF.RWNEG_O7': 70368744177664,
    'WQSF.RWNEG_O8': 140737488355328,
    'WQSF.RWNEG_O9': 281474976710656,
    'WQSF.RWNEG_O10': 562949953421312,
    'WQSF.RWNEG_O11': 1125899906842624,
    'WQSF.RWNEG_O12': 2251799813685248,
    'WQSF.RWNEG_O16': 4503599627370496,
    'WQSF.RWNEG_O17': 9007199254740992,
    'WQSF.RWNEG_O18': 18014398509481984,
    'WQSF.RWNEG_O21': 36028797018963968
}

IDEPIX_PROCESSOR_ID = 'IDEPIX'

IDEPIX_FLAG_BAND_NAME = 'pixel_classif_flags'

IDEPIX_BITMASK_FLAG_CODING_DICT = {
    'pixel_classif_flags.IDEPIX_INVALID': 1,
    'pixel_classif_flags.IDEPIX_CLOUD': 2,
    'pixel_classif_flags.IDEPIX_CLOUD_AMBIGUOUS': 4,
    'pixel_classif_flags.IDEPIX_CLOUD_SURE': 8,
    'pixel_classif_flags.IDEPIX_CLOUD_BUFFER': 16,
    'pixel_classif_flags.IDEPIX_CLOUD_SHADOW': 32,
    'pixel_classif_flags.IDEPIX_SNOW_ICE': 64,
    'pixel_classif_flags.IDEPIX_BRIGHT': 128,
    'pixel_classif_flags.IDEPIX_WHITE': 256,
    'pixel_classif_flags.IDEPIX_COASTLINE': 512,
    'pixel_classif_flags.IDEPIX_LAND': 1024,
    'pixel_classif_flags.IDEPIX_MOUNTAIN_SHADOW': 2048
}

DEFAULT_IDEPIX_VALID_PIXEL_EXPR = 'not pixel_classif_flags.IDEPIX_INVALID ' \
                                  'and not pixel_classif_flags.IDEPIX_CLOUD ' \
                                  'and not pixel_classif_flags.IDEPIX_CLOUD_AMBIGUOUS ' \
                                  'and not pixel_classif_flags.IDEPIX_CLOUD_SURE ' \
                                  'and not pixel_classif_flags.IDEPIX_CLOUD_BUFFER ' \
                                  'and not pixel_classif_flags.IDEPIX_CLOUD_SHADOW ' \
                                  'and not pixel_classif_flags.IDEPIX_MOUNTAIN_SHADOW '

SCDA_PROCESSOR_ID = 'SCDA'

SCDA_FLAG_BAND_NAME = 'scda_cloud_mask'

SCDA_BITMASK_FLAG_CODING_DICT = {
    'scda_cloud_mask.SCDA_INVALID': 1,
    'scda_cloud_mask.SCDA_CLEAR': 2,
    'scda_cloud_mask.SCDA_CLOUDY': 4,
}

DEFAULT_SCDA_VALID_PIXEL_EXPR = 'not scda_cloud_mask.SCDA_INVALID ' \
                                'and not scda_cloud_mask.SCDA_CLOUDY '
