import datetime
from math import ceil
import datetime
import os
import platform
import sys
import tempfile
import time

import esa_snappy
import numpy as np
# If a Java type is needed which is not imported by snappy by default it can be retrieved manually.
# First import jpy
from esa_snappy import jpy

Float = jpy.get_type('java.lang.Float')

import scda

SLSTR_S1_SOLAR_FLUX = 1837.39
SLSTR_S5_SOLAR_FLUX = 248.33


class Sice2ScdaOp:
    """
    The Sice2 GPF operator

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
        f = open(tempfile.gettempdir() + '/sice2_scda_op.log', 'w')

        sys.path.append(resource_root)

        f.write('Python module location: ' + __file__ + '\n')
        f.write('Python module location parent: ' + resource_root + '\n')

        print('platform.system(): ' + platform.system() + '\n')
        print('sys.version_info(): ' + str(sys.version_info) + '\n')
        print('sys.version_info > (3, 0): ' + str(sys.version_info > (3, 0)))

        self.source_product = context.getSourceProduct('l1bProduct')
        print('initialize: source product location is', self.source_product.getFileLocation())

        self.width = self.source_product.getSceneRasterWidth()
        self.height = self.source_product.getSceneRasterHeight()

        print('Start sice2_scda_op...')

        #####
        self.r550_radiance_band = self._get_band(self.source_product, "S1_radiance_an")
        self.r1600_radiance_band = self._get_band(self.source_product, "S5_radiance_an")
        self.bt37_band = self._get_band(self.source_product, "S7_BT_in")
        self.bt11_band = self._get_band(self.source_product, "S8_BT_in")
        self.bt12_band = self._get_band(self.source_product, "S9_BT_in")
        self.sza_tpg = self._get_band(self.source_product, "solar_zenith_tn")

        ##############################

        # Create the target product
        scda_product = esa_snappy.Product('py_SICE21_scda', 'py_SICE21_scda', self.width, self.height)
        esa_snappy.ProductUtils.copyGeoCoding(self.source_product, scda_product)
        esa_snappy.ProductUtils.copyMetadata(self.source_product, scda_product)
        scda_product.setStartTime(self.source_product.getStartTime())
        scda_product.setEndTime(self.source_product.getEndTime())

        context.setTargetProduct(scda_product)
        self.add_target_bands(scda_product)

        f.write('end initialize.')
        print('end initialize.')
        t_end = datetime.datetime.now()
        print('SNAPPY SICE2 processing time (seconds): ' + str((t_end - t_start).seconds))

        end_time = time.process_time()
        print('SNAPPY SICE2 processing time (CPU seconds): ' + str((end_time - start_time)))

        f.close()

    def add_target_bands(self, scda_product):
        """

        :param scda_product:
        :return:
        """
        self.scda_band = scda_product.addBand('scda_cloud_mask', esa_snappy.ProductData.TYPE_UINT8)
        self.scda_band.setDescription('SCDA binary cloud mask')
        self.scda_band.setNoDataValue(255)
        self.scda_band.setNoDataValueUsed(True)

        self.ndsi_band = scda_product.addBand('ndsi', esa_snappy.ProductData.TYPE_FLOAT32)
        self.ndsi_band.setDescription('NDSI: (r550 - r1600) / (r550 + r1600)')
        self.ndsi_band.setNoDataValue(Float.NaN)
        self.ndsi_band.setNoDataValueUsed(True)

    def computeTileStack(self, context, target_tiles, target_rectangle):

        r550_rad_tile = context.getSourceTile(self.r550_radiance_band, target_rectangle)
        r1600_rad_tile = context.getSourceTile(self.r1600_radiance_band, target_rectangle)
        bt37_tile = context.getSourceTile(self.bt37_band, target_rectangle)
        bt11_tile = context.getSourceTile(self.bt11_band, target_rectangle)
        bt12_tile = context.getSourceTile(self.bt12_band, target_rectangle)
        sza_tile = context.getSourceTile(self.sza_tpg, target_rectangle)

        rw = target_rectangle.width
        rh = target_rectangle.height
        r550_rad_data = np.array(r550_rad_tile.getSamplesFloat(), dtype=np.float32).reshape((rh, rw))
        r1600_rad_data = np.array(r1600_rad_tile.getSamplesFloat(), dtype=np.float32).reshape((rh, rw))
        bt37_data = np.array(bt37_tile.getSamplesFloat(), dtype=np.float32).reshape((rh, rw))
        bt11_data = np.array(bt11_tile.getSamplesFloat(), dtype=np.float32).reshape((rh, rw))
        bt12_data = np.array(bt12_tile.getSamplesFloat(), dtype=np.float32).reshape((rh, rw))
        sza_data = np.array(sza_tile.getSamplesFloat(), dtype=np.float32).reshape((rh, rw))

        r550_refl_data = (r550_rad_data * np.pi) / (SLSTR_S1_SOLAR_FLUX * np.cos(np.deg2rad(sza_data)))
        r1600_refl_data = (r1600_rad_data * np.pi) / (SLSTR_S5_SOLAR_FLUX * np.cos(np.deg2rad(sza_data)))

        scda_data, ndsi_data = scda.scda_v20(r550_refl_data, r1600_refl_data, bt37_data, bt11_data, bt12_data)

        # The target tiles which shall be filled with data are provided as parameter to this method
        # Set the results to the target tiles
        target_tiles.get(self.scda_band).setSamples(scda_data.flatten())
        target_tiles.get(self.ndsi_band).setSamples(ndsi_data.flatten())

    def dispose(self, operator):
        """
        The GPF dispose method. Nothing to do here.
        :param operator:
        :return:
        """
        pass

    @staticmethod
    def _get_band(input_product, band_name):
        """
        Gets band from input product by name
        :param input_product
        :param band_name
        :return:
        """
        band = input_product.getBand(band_name)
        if not band:
            band = input_product.getTiePointGrid(band_name)
            if not band:
                raise RuntimeError('Product has no band or tpg with name', band_name)
        return band
