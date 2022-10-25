import datetime

import sice2_v21_algo
# import sice2_v21_io
import sice2_v21_io2

import platform
import tempfile
import sys
import os
import time
# import numpy as np
# from numpy import genfromtxt
# import rasterio as rio

import xarray as xr

import esa_snappy
from esa_snappy import ProductIO

# If a Java type is needed which is not imported by snappy by default it can be retrieved manually.
# First import jpy
from esa_snappy import jpy

# and then import the type
import sice2_constants
import sice2_utils

Float = jpy.get_type('java.lang.Float')
Color = jpy.get_type('java.awt.Color')

NDSI_HIGH_THRESHOLD = 0.8
NDSI_LOW_THRESHOLD = -0.5


class Sice2V21TifdirsOp:
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

        print('Start sice.py')
        print('Input folder:', self.tif_input_directory)
        print('Output folder:', self.tif_output_directory)

        olci_reader = sice2_v21_io2.Sice2V21TifIo(self.tif_input_directory)
        olci_reader.open_tif()

        olci_scene = olci_reader.olci_scene

        start_time = time.process_time()

        if len(olci_scene.xy) < 1000000:
            # snow = process(olci_scene)
            snow = sice2_v21_algo.process(olci_scene)
        else:
            # snow = process_by_chunk(olci_scene, chunk_size=500000)
            snow = sice2_v21_algo.process_by_chunk(olci_scene, chunk_size=500000)

        duration = time.process_time() - start_time
        print("Time elapsed: ", duration)

        olci_reader.write_output(snow, self.tif_output_directory)

        # Create the target product
        # todo: we do not want a target product for the moment. Check how 'autoWriteDisabled' works with snappy
        snow_product = esa_snappy.Product('py_SNOW', 'py_SNOW', 1, 1)  # dummy
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


