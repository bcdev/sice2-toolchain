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
from esa_snappy import ProductIO

Float = jpy.get_type('java.lang.Float')

import scda


class Sice2ScdaTifdirsOp:
    """
    The Sice2 GPF operator

    Authors: O.Danne, 2022
    """

    def __init__(self):
        pass

    def initialize(self, context):
        """
        GPF initialize method

        :param context
        :return:
        """
        t_start = datetime.datetime.now()
        start_time = time.process_time()

        resource_root = os.path.dirname(__file__)
        f = open(tempfile.gettempdir() + '/sice2_scda_tifdirs_op.log', 'w')

        sys.path.append(resource_root)

        f.write('Python module location: ' + __file__ + '\n')
        f.write('Python module location parent: ' + resource_root + '\n')

        print('platform.system(): ' + platform.system() + '\n')
        print('sys.version_info(): ' + str(sys.version_info) + '\n')
        print('sys.version_info > (3, 0): ' + str(sys.version_info > (3, 0)))

        # Tif directory parameter defined in ndsi_op-info.xml
        tif_input_directory = context.getParameter('tifInputDir')

        print('Start sice2_scda_tifdirs_op...')
        print('Input folder:', tif_input_directory)

        self.called_compute_tile_stack = 0

        r550_product = ProductIO.readProduct(tif_input_directory + os.sep + "r_TOA_S1.tif")
        r1600_product = ProductIO.readProduct(tif_input_directory + os.sep + "r_TOA_S5.tif")
        bt37_product = ProductIO.readProduct(tif_input_directory + os.sep + "BT_S7.tif")
        bt11_product = ProductIO.readProduct(tif_input_directory + os.sep + "BT_S8.tif")
        bt12_product = ProductIO.readProduct(tif_input_directory + os.sep + "BT_S9.tif")

        self.width = r550_product.getSceneRasterWidth()
        self.height = r550_product.getSceneRasterHeight()

        self.r550_band = self._get_band(r550_product, "S1_reflectance_an")
        self.r1600_band = self._get_band(r1600_product, "S5_reflectance_an")
        self.bt37_band = self._get_band(bt37_product, "S7_BT_in")
        self.bt11_band = self._get_band(bt11_product, "S8_BT_in")
        self.bt12_band = self._get_band(bt12_product, "S9_BT_in")

        ##############################

        # Create the target product
        scda_product = esa_snappy.Product('py_SICE21_scda', 'py_SICE21_scda', self.width, self.height)
        esa_snappy.ProductUtils.copyGeoCoding(r550_product, scda_product)
        esa_snappy.ProductUtils.copyMetadata(r550_product, scda_product)
        scda_product.setStartTime(r550_product.getStartTime())
        scda_product.setEndTime(r550_product.getEndTime())

        context.setTargetProduct(scda_product)
        self.add_target_bands(scda_product)

        f.write('end initialize.')
        print('end initialize.')
        t_end = datetime.datetime.now()
        print('SNAPPY SCDA processing time (seconds): ' + str((t_end - t_start).seconds))

        end_time = time.process_time()
        print('SNAPPY SCDA processing time (CPU seconds): ' + str((end_time - start_time)))

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

        # write NDSI in standard op only
        # self.ndsi_band = snow_product.addBand('ndsi', esa_snappy.ProductData.TYPE_FLOAT32)
        # self.ndsi_band.setDescription('NDSI')
        # self.ndsi_band.setNoDataValue(Float.NaN)
        # self.ndsi_band.setNoDataValueUsed(True)

    def computeTileStack(self, context, target_tiles, target_rectangle):

        # debugging:
        # num_pixels = target_rectangle.width * target_rectangle.height
        # print('Call computeTileStack: num_pixels=' + str(num_pixels))
        # print('target_rectangle.x =' + str(target_rectangle.x))
        # print('target_rectangle.y =' + str(target_rectangle.y))
        # print('target_rectangle.width =' + str(target_rectangle.width))
        # print('target_rectangle.height =' + str(target_rectangle.height))
        # self.called_compute_tile_stack = self.called_compute_tile_stack + 1
        # print('Tile #' + str(self.called_compute_tile_stack))

        r550_tile = context.getSourceTile(self.r550_band, target_rectangle)
        r1600_tile = context.getSourceTile(self.r1600_band, target_rectangle)
        bt37_tile = context.getSourceTile(self.bt37_band, target_rectangle)
        bt11_tile = context.getSourceTile(self.bt11_band, target_rectangle)
        bt12_tile = context.getSourceTile(self.bt12_band, target_rectangle)

        rw = target_rectangle.width
        rh = target_rectangle.height
        r550_data = np.array(r550_tile.getSamplesFloat(), dtype=np.float32).reshape((rh, rw))
        r1600_data = np.array(r1600_tile.getSamplesFloat(), dtype=np.float32).reshape((rh, rw))
        bt37_data = np.array(bt37_tile.getSamplesFloat(), dtype=np.float32).reshape((rh, rw))
        bt11_data = np.array(bt11_tile.getSamplesFloat(), dtype=np.float32).reshape((rh, rw))
        bt12_data = np.array(bt12_tile.getSamplesFloat(), dtype=np.float32).reshape((rh, rw))

        scda_data, ndsi_data = scda.scda_v20(r550_data, r1600_data, bt37_data, bt11_data, bt12_data)

        # The target tiles which shall be filled with data are provided as parameter to this method
        # Set the results to the target tiles
        target_tiles.get(self.scda_band).setSamples(scda_data.flatten())

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
                if input_product.getNumBands() == 1:
                    # special case: check tif files with 1 band named 'band_1'
                    band = input_product.getBand("band_1")
                if not band:
                    raise RuntimeError('Product has no band or tpg with name', band_name)
        return band


