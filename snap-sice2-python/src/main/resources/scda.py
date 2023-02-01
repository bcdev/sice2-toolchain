# -*- coding: utf-8 -*-
"""

@author: Adrien Wehrlé, GEUS (Geological Survey of Denmark and Greenland)

Implementation of the Simple Cloud Detection Algorithm (SCDA) v2.0
using SLSTR bands, described in Fig. 5 of METSÄMÄKI et al, 2015.

METSÄMÄKI, Sari, PULLIAINEN, Jouni, SALMINEN, Miia, et al. Introduction 
to GlobSnow Snow Extent products with considerations for accuracy assessment. 
Remote Sensing of Environment, 2015, vol. 156, p. 96-108.

The original syntax has been preserved to easily link back to the description
of the algorithm.

@author: Olaf Danne, BC (Brockmann Consult): Technical adaptations for usage from SNAP operators (01/2023).

"""

import numpy as np
from numpy import asarray as ar


def radiometric_calibration(r16_data):
    """
    Sentinel-3 Product Notice – SLSTR:
    "Based on the analysis performed to-date, a recommendation has been put forward to users to
    adjust the S5 and S6 reflectances by factors of 1.12 and 1.20 respectively in the nadir view and
    1.15 and 1.26 in the oblique view. Uncertainty estimates on these differences are still to be
    evaluated and comparisons with other techniques have yet to be included."

    INPUTS:
        r16_data: numpy array: Dataset for Top of Atmosphere (TOA) reflectance channel S5.
             Central wavelengths at 1.6um. [rasterio.io.DatasetReader]

    OUTPUTS:
        r16_data_rc: numpy array: Adjusted Top of Atmosphere (TOA) reflectance for channel S5.
    """

    factor = 1.12
    return r16_data * factor


def scda_v20(r550, r16, bt37, bt11, bt12):
    """

    INPUTS: numpy arrays for:
        R550, R16: Top of Atmosphere (TOA) reflectances for channels S1 and S5.
                   Central wavelengths at 550nm and 1.6um. [arrays]
        BT37, BT11, BT12: Gridded pixel Brightness Temperatures (BT) for channels
                          S7, S8 and S9 (1km TIR grid, nadir view). Central
                          wavelengths at 3.7, 11 and 12 um. [arrays]

    OUTPUTS: numpy arrays for:
        {inpath}/scda: Simple Cloud Detection Algorithm (SCDA) results (clouds=1, clear=0)
        {inpath}/ndsi: Normalized Difference Snow Index (NDSI)

    """

    # radiometric calibration of R16:
    r16_rc = radiometric_calibration(r16)

    # determining the NDSI, needed for the cloud detection
    ndsi = (r550 - r16_rc) / (r550 + r16_rc)

    # initializing thresholds
    base = np.empty((r550.shape[0], r550.shape[1]))
    thr = base.copy()
    thr[:] = np.nan
    th_rmax = base.copy()
    th_rmax[:] = -5.5
    s = base.copy()
    s[:] = 1.1

    # masking nan values
    mask_invalid = np.isnan(r550)

    # tests 1 to 5, only based on inputs
    t1 = ar(r550 > 0.30) * ar(ndsi / r550 < 0.8) * ar(bt12 <= 290)
    t2 = ar(bt11 - bt37 < -13) * ar(r550 > 0.15) * ar(ndsi >= -0.30) * ar(r16_rc > 0.10) * ar(bt12 <= 293)
    t3 = ar(bt11 - bt37 < -30)
    t4 = ar(r550 < 0.75) * ar(bt12 > 265)
    t5 = ar(r550 > 0.75)

    cloud_detection = t1
    cloud_detection[cloud_detection == False] = t2[cloud_detection == False]
    cloud_detection[cloud_detection == False] = t3[cloud_detection == False]

    thr1 = 0.5 * bt12 - 133

    th_rmax[t4 == False] = -8
    thr = np.minimum(thr1, th_rmax)
    s[t5 == False] = 1.5

    # test 6, based on fluctuating thresholds
    t6_1 = ar(bt11 - bt37 < thr)
    t6_2 = ar(ndsi / r550 < s)
    t6_3 = ar((ndsi >= -0.02) & (ndsi <= 0.75))
    t6_4 = ar(bt12 <= 270)
    t6_5 = ar(r550 > 0.18)

    t6 = t6_1 * t6_2 * t6_3 * t6_4 * t6_5

    cloud_detection[cloud_detection == False] = t6[cloud_detection == False]

    # set cloudy to 1, clear to 0, and invalid to 255:
    cloud_detection = np.where(cloud_detection == True, 1.0, 0.0)
    cloud_detection[mask_invalid] = 255.0

    return cloud_detection, ndsi
