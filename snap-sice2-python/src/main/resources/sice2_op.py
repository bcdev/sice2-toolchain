import datetime

import sice2_algo
import platform
import tempfile
import sys
import os
import time
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
import rasterio as rio

import snappy
from snappy import ProductIO
from snappy import FlagCoding

# If a Java type is needed which is not imported by snappy by default it can be retrieved manually.
# First import jpy
from snappy import jpy

# and then import the type
import sice2_constants
import sice2_utils

Float = jpy.get_type('java.lang.Float')
Color = jpy.get_type('java.awt.Color')

NDSI_HIGH_THRESHOLD = 0.8
NDSI_LOW_THRESHOLD = -0.5


class Sice2Op:
    """
    The Sice2 GPF operator
    TODO: set up repository!!

    Authors: O.Danne, 2022
    """

    def __init__(self):
        pass

    def initialize(self, context):
        """
        GPF initialize method

        :param operator
        :return:
        """
        t_start = datetime.datetime.now()
        start_time = time.process_time()

        resource_root = os.path.dirname(__file__)
        f = open(tempfile.gettempdir() + '/sice2_py.log', 'w')

        sys.path.append(resource_root)

        f.write('Python module location: ' + __file__ + '\n')
        f.write('Python module location parent: ' + resource_root + '\n')

        print('platform.system(): ' + platform.system() + '\n')
        print('sys.version_info(): ' + str(sys.version_info) + '\n')
        print('sys.version_info > (3, 0): ' + str(sys.version_info > (3, 0)))

        # Tif directory parameter defined in ndsi_op-info.xml
        self.tif_input_directory = context.getParameter('tifInputDir')
        self.tif_output_directory = context.getParameter('tifOutputDir')

        # Via the context object the source product which shall be processed can be retrieved
        # source_product = context.getSourceProduct('source')
        # source_product = snappy.Product('dummy', 'dummy', 1, 1)
        # width = source_product.getSceneRasterWidth()
        # height = source_product.getSceneRasterHeight()
        # TODO: setup list of source products (all tif products in tif test dir)
        # print('initialize: source product location is', source_product.getFileLocation())

        # print('initialize: tif_directory =', self.tif_directory)
        # tif_files = os.listdir(self.tif_directory)

        np.seterr(invalid='ignore')
        start_time = time.process_time()

        # open first rtoa tif to get metadata:
        rtoa_01_tif_product = rio.open(self.tif_input_directory + os.sep + 'r_TOA_01.tif')
        tif_metadata = rtoa_01_tif_product.meta
        self.product_width = tif_metadata['width']
        self.product_height = tif_metadata['height']

        rtoa_data = np.full([sice2_constants.OLCI_NUM_SPECTRAL_BANDS, self.product_height, self.product_width], np.nan)
        for i in range(sice2_constants.OLCI_NUM_SPECTRAL_BANDS):
            try:
                rtoa_tif_product = rio.open(
                    (self.tif_input_directory + os.sep + 'r_TOA_' + str(i + 1).zfill(2) + '.tif'))
                rtoa_data[i, :, :] = rtoa_tif_product.read(1).astype('float32')
            except:
                rtoa_data[i, :, :] = np.nan

        sza_data = rio.open(self.tif_input_directory + os.sep + 'SZA.tif').read(1).astype('float32')
        saa_data = rio.open(self.tif_input_directory + os.sep + 'SAA.tif').read(1).astype('float32')
        oza_data = rio.open(self.tif_input_directory + os.sep + 'OZA.tif').read(1).astype('float32')
        oaa_data = rio.open(self.tif_input_directory + os.sep + 'OAA.tif').read(1).astype('float32')
        ozone_data = rio.open(self.tif_input_directory + os.sep + 'O3.tif').read(1).astype('float32')
        wv_data = rio.open(self.tif_input_directory + os.sep + 'WV.tif').read(1).astype('float32')
        height_data = rio.open(self.tif_input_directory + os.sep + 'height.tif').read(1).astype('float32')

        # todo: clarify if we really want this:
        sza_data[np.isnan(rtoa_data[0, :, :])] = np.nan
        saa_data[np.isnan(rtoa_data[0, :, :])] = np.nan
        oza_data[np.isnan(rtoa_data[0, :, :])] = np.nan
        oaa_data[np.isnan(rtoa_data[0, :, :])] = np.nan

        water_vod = genfromtxt(resource_root + os.sep + 'auxdata' + os.sep + 'tg_water_vod.dat', delimiter='   ')
        self.voda = water_vod[range(sice2_constants.OLCI_NUM_SPECTRAL_BANDS), 1]

        ozone_vod = genfromtxt(resource_root + os.sep + 'auxdata' + os.sep + 'tg_vod.dat', delimiter='   ')
        self.tozon = ozone_vod[range(sice2_constants.OLCI_NUM_SPECTRAL_BANDS), 1]

        # As it is always a good idea to separate responsibilities the algorithmic methods are put
        # into an other class
        self.sice2algo = sice2_algo.Sice2Algo(NDSI_LOW_THRESHOLD, NDSI_HIGH_THRESHOLD)

        # =========== ozone scattering  ====================================
        # variable renaming from breadboard:
        # BXXX --> o3_sice
        # toa_cor_o3 --> toa_cor_o3
        o3_sice, rtoa_cor_o3 = self.sice2algo.ozone_scattering(ozone_data, self.tozon, sza_data, oza_data, rtoa_data)

        # Filtering pixels unsuitable for retrieval
        isnow = np.full(sza_data.shape, np.nan, dtype=float)
        isnow[sza_data > 75] = 100
        isnow[rtoa_cor_o3[20, :, :] < 0.1] = 102

        for i_channel in range(sice2_constants.OLCI_NUM_SPECTRAL_BANDS):
            rtoa_cor_o3[i_channel, ~np.isnan(isnow)] = np.nan

        sza_data[~np.isnan(isnow)] = np.nan
        saa_data[~np.isnan(isnow)] = np.nan
        oza_data[~np.isnan(isnow)] = np.nan
        oaa_data[~np.isnan(isnow)] = np.nan
        height_data[~np.isnan(isnow)] = np.nan

        # print('sza_data[500,500]: ' + str(sza_data[500,500]))

        # =========== view geometry and atmosphere propeties  ==============
        raa, am1, am2, ak1, ak2, amf, co = self.sice2algo.view_geometry(oaa_data, saa_data, sza_data, oza_data)
        tau, p, g, gaer, taumol, tauaer = self.sice2algo.aerosol_properties(sice2_constants.AOT, height_data, co)

        # =========== snow properties  ====================================

        D, area, al, r0, bal = self.sice2algo.snow_properties(rtoa_cor_o3, ak1, ak2)
        # filtering small D
        D_thresh = 0.1
        isnow[D < D_thresh] = 104

        for i in range(sice2_constants.OLCI_NUM_SPECTRAL_BANDS):
            rtoa_cor_o3[i, D < D_thresh] = np.nan

        area[D < D_thresh] = np.nan
        al[D < D_thresh] = np.nan
        r0[D < D_thresh] = np.nan
        bal[D < D_thresh] = np.nan
        am1[D < D_thresh] = np.nan
        am2[D < D_thresh] = np.nan

        # =========== clean snow  ====================================

        # for that we calculate the theoretical reflectance at band 1 of a surface with:
        # r0 = 1, a (albedo) = 1, ak1 = 1, ak2 = 1
        # t1 and t2 are the backscattering fraction
        t1, t2, ratm, r, astra, rms = self.sice2algo.prepare_coef(tau, g, p, am1, am2, amf, gaer,
                                                      taumol, tauaer)
        rs_1 = self.sice2algo.alb2rtoa(1, t1[0, :, :], t2[0, :, :], np.ones_like(r0), np.ones_like(ak1),
                           np.ones_like(ak2), ratm[0, :, :], r[0, :, :])

        # we then compare it to the observed toa[0] value
        ind_clean = rtoa_cor_o3[0, :, :] >= rs_1
        isnow[ind_clean] = 0

        # STEP 4a: clean snow retrieval
        # the spherical albedo derivation: alb_sph

        alb_sph = np.exp(-np.sqrt(1000. * 4. * np.pi
                                  * sice2_utils.mult_channel(sice2_constants.bai / sice2_constants.w,
                                                             np.tile(al, (sice2_constants.OLCI_NUM_SPECTRAL_BANDS, 1, 1)))))
        alb_sph[alb_sph > 0.999] = 1

        # ========== very dirty snow  ====================================

        ind_pol = rtoa_cor_o3[0, :, :] < rs_1

        isnow[ind_pol] = 1

        ind_very_dark = np.logical_and(rtoa_cor_o3[20] < 0.4, ind_pol)
        isnow[ind_very_dark] = 6

        am11 = np.sqrt(1. - am1[ind_very_dark] ** 2.)
        am12 = np.sqrt(1. - am2[ind_very_dark] ** 2.)

        tz = np.arccos(-am1[ind_very_dark] * am2[ind_very_dark] + am11 * am12
                       * np.cos(raa[ind_very_dark] * 3.14159 / 180.)) * 180. / np.pi

        pz = 11.1 * np.exp(-0.087 * tz) + 1.1 * np.exp(-0.014 * tz)

        rclean = 1.247 + 1.186 * (am1[ind_very_dark] + am2[ind_very_dark]) \
                 + 5.157 * am1[ind_very_dark] * am2[ind_very_dark] + pz

        rclean = rclean / 4. / (am1[ind_very_dark] + am2[ind_very_dark])
        r0[ind_very_dark] = rclean

        # =========== polluted snow  ====================================

        ind_pol = np.logical_or(ind_very_dark, ind_pol)

        if np.any(ind_pol):
            subs_pol = np.argwhere(ind_pol)

            # approximation of the transcendental equation allowing closed-from solution
            # alb_sph[:, ind_pol] = (toa_cor_o3[:, ind_pol] - r[:, ind_pol]) \
            # /(t1[:,ind_pol]*t2[:,ind_pol]*r0[ind_pol] + ratm[:,ind_pol]*(toa_cor_o3[:,ind_pol] - r[:,ind_pol]))

            # solving iteratively the transcendental equation
            alb_sph[:, ind_pol] = 1

            def solver_wrapper(toa_cor_o3, tau, t1, t2, r0, ak1, ak2, ratm, r):
                def func_solv(albedo):
                    return toa_cor_o3 - self.sice2algo.alb2rtoa(albedo, t1, t2, r0, ak1, ak2, ratm, r)
                # it is assumed that albedo is in the range 0.1-1.0

                return self.sice2algo.zbrent(func_solv, 0.1, 1, 100, 1.e-6)

            solver_wrapper_v = np.vectorize(solver_wrapper)
            # loop over all bands except band 19, 20
            for i_channel in np.append(np.arange(18), [20]):

                alb_sph[i_channel, ind_pol] = solver_wrapper_v(
                    rtoa_cor_o3[i_channel, ind_pol], tau[i_channel, ind_pol],
                    t1[i_channel, ind_pol], t2[i_channel, ind_pol],
                    r0[ind_pol], ak1[ind_pol], ak2[ind_pol], ratm[i_channel, ind_pol],
                    r[i_channel, ind_pol])

                ind_bad = alb_sph[i_channel, :, :] == -999
                alb_sph[i_channel, ind_bad] = np.nan
                isnow[ind_bad] = -i_channel

            # INTERNal CHECK FOR CLEAN PIXELS
            # Are reprocessed as clean
            ind_clear_pol1 = np.logical_and(ind_pol, alb_sph[0, :, :] > 0.98)
            ind_clear_pol2 = np.logical_and(ind_pol, alb_sph[1, :, :] > 0.98)
            ind_clear_pol = np.logical_or(ind_clear_pol1, ind_clear_pol2)
            isnow[ind_clear_pol] = 7

            for i_channel in range(sice2_constants.OLCI_NUM_SPECTRAL_BANDS):
                alb_sph[i_channel, ind_clear_pol] = np.exp(-np.sqrt(4. * 1000.
                                                                    * al[ind_clear_pol]
                                                                    * np.pi * sice2_constants.bai[i_channel]
                                                                    / sice2_constants.w[i_channel]))

                # re-defining polluted pixels
            ind_pol = np.logical_and(ind_pol, isnow != 7)

            # retrieving snow impurities
            ntype, bf, conc = self.sice2algo.snow_impurities(alb_sph, bal)

            # alex   09.06.2019
            # reprocessing of albedo to remove gaseous absorption using linear polynomial
            # approximation in the range 753-778nm.
            # Meaning: alb_sph[12],alb_sph[13] and alb_sph[14] are replaced by a linear
            # interpolation between alb_sph[11] and alb_sph[15]
            afirn = (alb_sph[15, ind_pol] - alb_sph[11, ind_pol]) / (sice2_constants.w[15] - sice2_constants.w[11])
            bfirn = alb_sph[15, ind_pol] - afirn * sice2_constants.w[15]
            alb_sph[12, ind_pol] = bfirn + afirn * sice2_constants.w[12]
            alb_sph[13, ind_pol] = bfirn + afirn * sice2_constants.w[13]
            alb_sph[14, ind_pol] = bfirn + afirn * sice2_constants.w[14]

            # BAV 09-02-2020: 0.5 to 0.35
            # pixels that are clean enough in channels 18 19 20 and 21 are not affected
            # by pollution, the analytical equation can then be used
            ind_ok = np.logical_and(ind_pol, rtoa_cor_o3[20, :, :] > 0.35)

            for i_channel in range(17, sice2_constants.OLCI_NUM_SPECTRAL_BANDS):
                alb_sph[i_channel, ind_ok] = np.exp(-np.sqrt(4. * 1000. * al[ind_ok]
                                                             * np.pi * sice2_constants.bai[i_channel]
                                                             / sice2_constants.w[i_channel]))
            # Alex, SEPTEMBER 26, 2019
            # to avoid the influence of gaseous absorption (water vapor) we linearly
            # interpolate in the range 885-1020nm for bare ice cases only (low toa[20])
            # Meaning: alb_sph[18] and alb_sph[19] are replaced by a linear interpolation
            # between alb_sph[17] and alb_sph[20]
            delx = sice2_constants.w[20] - sice2_constants.w[17]
            bcoef = (alb_sph[20, ind_pol] - alb_sph[17, ind_pol]) / delx
            acoef = alb_sph[20, ind_pol] - bcoef * sice2_constants.w[20]
            alb_sph[18, ind_pol] = acoef + bcoef * sice2_constants.w[18]
            alb_sph[19, ind_pol] = acoef + bcoef * sice2_constants.w[19]

        # ========= derivation of plane albedo and reflectance ===========

        rp = np.power(alb_sph, ak1)
        refl = r0 * np.power(alb_sph, (ak1 * ak2 / r0))

        ind_all_clean = np.logical_or(ind_clean, isnow == 7)

        # CalCULATION OF BBA of clean snow

        # approximation
        # planar albedo
        # rp1 and rp2 not derived anymore
        rp3 = np.full(sza_data.shape, np.nan, dtype=float)
        rp3[ind_all_clean] = self.sice2algo.plane_albedo_sw_approx(D[ind_all_clean],
                                                       am1[ind_all_clean])
        # spherical albedo
        # rs1 and rs2 not derived anymore
        rs3 = np.full(sza_data.shape, np.nan, dtype=float)
        rs3[ind_all_clean] = self.sice2algo.spher_albedo_sw_approx(D[ind_all_clean])

        # calculation of the BBA for the polluted snow
        rp1 = np.full(sza_data.shape, np.nan, dtype=float)
        rs1 = np.full(sza_data.shape, np.nan, dtype=float)
        rp2 = np.full(sza_data.shape, np.nan, dtype=float)
        rs2 = np.full(sza_data.shape, np.nan, dtype=float)
        rp1[ind_pol], rp2[ind_pol], rp3[ind_pol] = self.sice2algo.BBA_calc_pol(
            rp[:, ind_pol], sice2_constants.asol, sice2_constants.sol1_pol, sice2_constants.sol2, sice2_constants.sol3_pol)
        rs1[ind_pol], rs2[ind_pol], rs3[ind_pol] = self.sice2algo.BBA_calc_pol(
            alb_sph[:, ind_pol], sice2_constants.asol, sice2_constants.sol1_pol, sice2_constants.sol2, sice2_constants.sol3_pol)

        # =========== Output  ====================================
        self.write_output_tif('O3_SICE', o3_sice, self.tif_output_directory, tif_metadata)
        self.write_output_tif('grain_diameter', D, self.tif_output_directory, tif_metadata)
        self.write_output_tif('snow_specific_surface_area', area, self.tif_output_directory, tif_metadata)
        self.write_output_tif('al', al, self.tif_output_directory, tif_metadata)
        self.write_output_tif('r0', r0, self.tif_output_directory, tif_metadata)
        self.write_output_tif('diagnostic_retrieval', isnow, self.tif_output_directory, tif_metadata)
        self.write_output_tif('conc', conc, self.tif_output_directory, tif_metadata)
        self.write_output_tif('albedo_bb_planar_sw', rp3, self.tif_output_directory, tif_metadata)
        self.write_output_tif('albedo_bb_spherical_sw', rs3, self.tif_output_directory, tif_metadata)

        for i in np.append(np.arange(11), np.arange(15, sice2_constants.OLCI_NUM_SPECTRAL_BANDS)):
            self.write_output_tif('albedo_spectral_spherical_'+ str(i + 1).zfill(2), alb_sph[i, :, :],
                                  self.tif_output_directory, tif_metadata)
            self.write_output_tif('albedo_spectral_planar_'+ str(i + 1).zfill(2), rp[i, :, :],
                                  self.tif_output_directory, tif_metadata)
            self.write_output_tif('rBRR_'+ str(i + 1).zfill(2), refl[i, :, :],
                                  self.tif_output_directory, tif_metadata)

        # Create the target product
        # todo: we do not want a target product. Check how 'autoWriteDisabled' works with snappy
        snow_product = snappy.Product('py_SNOW', 'py_SNOW', 1, 1)  # dummy
        #
        # snappy.ProductUtils.copyGeoCoding(source_product, snow_product)
        # snappy.ProductUtils.copyMetadata(source_product, snow_product)
        # For copying the time information no helper method exists yet, but will come in SNAP 5.0
        # snow_product.setStartTime(source_product.getStartTime())
        # snow_product.setEndTime(source_product.getEndTime())

        # test: compare rio with snappy ProductIO --> rio is ~10x faster...
        t1 = datetime.datetime.now()
        self.sza_band = self._get_source_band_from_geotif_product(self.tif_input_directory, 'SZA.tif', 'band_1')
        self.product_width = self.sza_band.getRasterWidth()
        self.product_height = self.sza_band.getRasterHeight()
        sza_data = np.zeros(self.product_width * self.product_height, np.float32)
        self.sza_band.readPixels(0, 0, self.product_width, self.product_height, sza_data)
        t2 = datetime.datetime.now()
        print('GPT t2 - t1 (microseconds): ' + str((t2 - t1).microseconds))

        t1 = datetime.datetime.now()
        saa_test = rio.open(self.tif_input_directory + os.sep + 'SAA.tif').read(1).astype('float32')
        t2 = datetime.datetime.now()
        print('RIO t2 - t1 (microseconds): ' + str((t2 - t1).microseconds))

        # Provide the created target product to the framework so the computeTileStack method can be called
        # properly and the data can be written to disk.
        context.setTargetProduct(snow_product)
        f.write('end initialize.')
        print('end initialize.')
        t_end = datetime.datetime.now()
        print('SNAPPY SICE2 processing time (seconds): ' + str((t_end - t_start).seconds))

        end_time = time.process_time()
        print('SNAPPY SICE2 processing time (CPU seconds): ' + str((end_time - start_time)))

        f.close()

    def computeTileStack(self, context, target_tiles, target_rectangle):
        pass

    def dispose(self, operator):
        """
        The GPF dispose method. Nothing to do here.
        :param operator:
        :return:
        """
        pass

    def _get_source_band_from_geotif_product(self, source_dir, tif_file_name, band_name):
        """
        Gets band from input geotif product by name
        :param source_dir - the source product directory
        :param tif_file_name
        :param band_name

        :return: band
        """

        source_file_path = os.path.join(self.tif_input_directory, tif_file_name)
        if not os.path.exists(source_file_path):
            raise RuntimeError('TIF source product ' + source_file_path + ' does not exist.')

        geotif_product = ProductIO.readProduct(source_file_path)
        band = geotif_product.getBand(band_name)
        if not band:
            raise RuntimeError('Product has no band with name: ' + band_name)

        return band

    def write_output_tif(self, var_name, var_data, output_folder, metadata):
        """
        Writes output tif file for given variable
        :param var_name: variable/file name
        :param var_data: variable data (numpy array)
        :param output_folder: putput folder
        :param metadata
        :return:
        """

        with rio.open(output_folder + os.sep + var_name + '.tif', 'w+', **metadata) as dst:
            dst.write(var_data.astype('float32'), 1)
