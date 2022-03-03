import sice2_algo
import platform
import tempfile
import sys
import os
import numpy
import snappy
from snappy import FlagCoding

# If a Java type is needed which is not imported by snappy by default it can be retrieved manually.
# First import jpy
from snappy import jpy

# and then import the type
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
        resource_root = os.path.dirname(__file__)
        f = open(tempfile.gettempdir() + '/sice2_py.log', 'w')

        sys.path.append(resource_root)

        f.write('Python module location: ' + __file__ + '\n')
        f.write('Python module location parent: ' + resource_root + '\n')

        print('platform.system(): ' + platform.system() + '\n')
        print('sys.version_info(): ' + str(sys.version_info) + '\n')
        print ('sys.version_info > (3, 0): ' + str(sys.version_info > (3, 0)))

        # Tif directory parameter defined in ndsi_op-info.xml
        self.tif_directory = context.getParameter('tifDir')
        # TODO: generate test data: run script https://github.com/GEUS-SICE/SICE/blob/master/S3.xml
        # on suitable OLCI product from S3Snow

        # Via the context object the source product which shall be processed can be retrieved
        # source_product = context.getSourceProduct('source')
        source_product = snappy.Product('dummy', 'dummy', 1, 1)
        width = source_product.getSceneRasterWidth()
        height = source_product.getSceneRasterHeight()
        # TODO: setup list of source products (all tif products in tif test dir)
        # print('initialize: source product location is', source_product.getFileLocation())

        print('initialize: tif_directory =', self.tif_directory)

        # As it is always a good idea to separate responsibilities the algorithmic methods are put
        # into an other class
        self.sice2algo = sice2_algo.Sice2Algo(NDSI_LOW_THRESHOLD, NDSI_HIGH_THRESHOLD)

        # Create the target product
        snow_product = snappy.Product('py_SNOW', 'py_SNOW', width, height)
        # snappy.ProductUtils.copyGeoCoding(source_product, snow_product)
        # snappy.ProductUtils.copyMetadata(source_product, snow_product)
        # For copying the time information no helper method exists yet, but will come in SNAP 5.0
        # snow_product.setStartTime(source_product.getStartTime())
        # snow_product.setEndTime(source_product.getEndTime())

        # Adding new bands to the target product is straight forward.
        self.ndsi_band = snow_product.addBand('ndsi', snappy.ProductData.TYPE_FLOAT32)
        self.ndsi_band.setDescription('The Normalized Difference Snow Index')
        self.ndsi_band.setNoDataValue(Float.NaN)
        self.ndsi_band.setNoDataValueUsed(True)
        self.ndsi_flags_band = snow_product.addBand('ndsi_flags', snappy.ProductData.TYPE_UINT8)
        self.ndsi_flags_band.setDescription('The flag information')

        # Create a flagCoding for the flag band. This helps to display the information properly within SNAP.
        ndsiFlagCoding = FlagCoding('ndsi_flags')
        # The NDSI_LOW flag shall be at bit position 0 and has therefor the value 1, NDSI_HIGH has the
        # value 2 (bit 1) and so one
        low_flag = ndsiFlagCoding.addFlag("NDSI_LOW", 1, "NDSI below " + str(NDSI_LOW_THRESHOLD))
        high_flag = ndsiFlagCoding.addFlag("NDSI_HIGH", 2, "NDSI above " + str(NDSI_HIGH_THRESHOLD))
        neg_flag = ndsiFlagCoding.addFlag("NDSI_NEG", 4, "NDSI negative")
        pos_flag = ndsiFlagCoding.addFlag("NDSI_POS", 8, "NDSI positive")
        snow_product.getFlagCodingGroup().add(ndsiFlagCoding)
        self.ndsi_flags_band.setSampleCoding(ndsiFlagCoding)

        # Also for each flag a layer should be created
        snow_product.addMask('mask_' + low_flag.getName(), 'ndsi_flags.' + low_flag.getName(),
                             low_flag.getDescription(), Color.YELLOW, 0.3)
        snow_product.addMask('mask_' + high_flag.getName(), 'ndsi_flags.' + high_flag.getName(),
                             high_flag.getDescription(), Color.GREEN, 0.3)
        snow_product.addMask('mask_' + neg_flag.getName(), 'ndsi_flags.' + neg_flag.getName(),
                             neg_flag.getDescription(), Color(255, 0, 0), 0.3)
        snow_product.addMask('mask_' + pos_flag.getName(), 'ndsi_flags.' + pos_flag.getName(),
                             pos_flag.getDescription(), Color.BLUE, 0.3)

        # Provide the created target product to the framework so the computeTileStack method can be called
        # properly and the data can be written to disk.
        context.setTargetProduct(snow_product)
        f.write('end initialize.')
        f.close()

    def computeTileStack(self, context, target_tiles, target_rectangle):
        # The operator is asked by the framework to provide the data for a rectangle when the data is needed.
        # The required source data for the computation can be retrieved by getSourceTile(...) via the context object.
        # lower_tile = context.getSourceTile(self.lower_band, target_rectangle)
        # upper_tile = context.getSourceTile(self.upper_band, target_rectangle)

        # The actual data can be retrieved from the tiles by getSampleFloats(), getSamplesDouble() or getSamplesInt()
        # lower_samples = lower_tile.getSamplesFloat()
        # upper_samples = upper_tile.getSamplesFloat()
        # Values at specific pixel locations can be retrieved for example by lower_tile.getSampleFloat(x, y)

        # Convert the data into numpy data. It is easier and faster to work with as if you use plain python arithmetic
        # lower_data = numpy.array(lower_samples, dtype=numpy.float32) * self.lower_factor
        # upper_data = numpy.array(upper_samples, dtype=numpy.float32) * self.upper_factor
        lower_data = numpy.array([0.8])  # test!
        upper_data = numpy.array([0.7])

        # Doing the actual computation
        ndsi = self.sice2algo.compute_ndsi(lower_data, upper_data)
        ndsi_flags = self.sice2algo.compute_flags(ndsi)

        # The target tile which shall be filled with data are provided as parameter to this method
        ndsi_tile = target_tiles.get(self.ndsi_band)
        ndsi_flags_tile = target_tiles.get(self.ndsi_flags_band)

        # Set the result to the target tiles
        ndsi_tile.setSamples(ndsi)
        ndsi_flags_tile.setSamples(ndsi_flags)

    def dispose(self, operator):
        """
        The GPF dispose method. Nothing to do here.
        :param operator:
        :return:
        """
        pass

    def _get_band(self, input_product, band_name):
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

