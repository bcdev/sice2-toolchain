import sice2_test_algo
import sice2_constants
import platform
import tempfile
import sys
import os
import numpy
# import snappy
sys.path.append(sice2_constants.SNAPPY_DIR)
import snappy_esa
# from snappy import FlagCoding
from snappy_esa import FlagCoding

NDVI_HIGH_THRESHOLD = 0.8
NDVI_LOW_THRESHOLD = -0.5

# If a Java type is needed which is not imported by snappy by default it can be retrieved manually.
# First import jpy
# from snappy import jpy
from snappy_esa import jpy

# and then import the type
Float = jpy.get_type('java.lang.Float')
Color = jpy.get_type('java.awt.Color')

import xarray as xr


class Sice2TestOp:
    """
    TEST operator - does an NDVI only (from SNAP examples)
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
        sys.path.append(sice2_constants.SNAPPY_DIR)

        f.write('Python module location: ' + __file__ + '\n')
        f.write('Python module location parent: ' + resource_root + '\n')

        print('platform.system(): ' + platform.system() + '\n')
        print('sys.version_info(): ' + str(sys.version_info) + '\n')
        print ('sys.version_info > (3, 0): ' + str(sys.version_info > (3, 0)))

        # Via the context object the source product which shall be processed can be retrieved
        source_product = context.getSourceProduct('source')
        print('initialize: source product location is', source_product.getFileLocation())

        width = source_product.getSceneRasterWidth()
        height = source_product.getSceneRasterHeight()

        # Retrieve a parameters defined in ndvi_op-info.xml
        lower_band_name = context.getParameter('lowerName')
        self.lower_factor = context.getParameter('lowerFactor')
        upper_band_name = context.getParameter('upperName')
        self.upper_factor = context.getParameter('upperFactor')

        self.lower_band = self._get_band(source_product, lower_band_name)
        self.upper_band = self._get_band(source_product, upper_band_name)

        print('initialize: lower_band =', self.lower_band, ', upper_band =', self.upper_band)
        print('initialize: lower_factor =', self.lower_factor, ', upper_factor =', self.upper_factor)

        # As it is always a good idea to separate responsibilities the algorithmic methods are put
        # into an other class
        self.sice2algo = sice2_test_algo.Sice2TestAlgo(NDVI_LOW_THRESHOLD, NDVI_HIGH_THRESHOLD)

        # Create the target product
        # ndvi_product = snappy.Product('py_NDVI', 'py_NDVI', width, height)
        ndvi_product = snappy_esa.Product('py_NDVI', 'py_NDVI', width, height)
        # ProductUtils provides several useful helper methods to build the target product.
        # In most cases it is sufficient to copy the information from the source to the target.
        # That's why mainly copy methods exist like copyBand(...), copyGeoCoding(...), copyMetadata(...)
        # snappy.ProductUtils.copyGeoCoding(source_product, ndvi_product)
        snappy_esa.ProductUtils.copyGeoCoding(source_product, ndvi_product)
        # snappy.ProductUtils.copyMetadata(source_product, ndvi_product)
        snappy_esa.ProductUtils.copyMetadata(source_product, ndvi_product)
        # For copying the time information no helper method exists yet, but will come in SNAP 5.0
        ndvi_product.setStartTime(source_product.getStartTime())
        ndvi_product.setEndTime(source_product.getEndTime())

        # Adding new bands to the target product is straight forward.
        # self.ndvi_band = ndvi_product.addBand('ndvi', snappy.ProductData.TYPE_FLOAT32)
        self.ndvi_band = ndvi_product.addBand('ndvi', snappy_esa.ProductData.TYPE_FLOAT32)
        self.ndvi_band.setDescription('The Normalized Difference Vegetation Index')
        self.ndvi_band.setNoDataValue(Float.NaN)
        self.ndvi_band.setNoDataValueUsed(True)
        # self.ndvi_flags_band = ndvi_product.addBand('ndvi_flags', snappy.ProductData.TYPE_UINT8)
        self.ndvi_flags_band = ndvi_product.addBand('ndvi_flags', snappy_esa.ProductData.TYPE_UINT8)
        self.ndvi_flags_band.setDescription('The flag information')

        # Create a flagCoding for the flag band. This helps to display the information properly within SNAP.
        ndviFlagCoding = FlagCoding('ndvi_flags')
        # The NDVI_LOW flag shall be at bit position 0 and has therefor the value 1, NDVI_HIGH has the
        # value 2 (bit 1) and so one
        low_flag = ndviFlagCoding.addFlag("NDVI_LOW", 1, "NDVI below " + str(NDVI_LOW_THRESHOLD))
        high_flag = ndviFlagCoding.addFlag("NDVI_HIGH", 2, "NDVI above " + str(NDVI_HIGH_THRESHOLD))
        neg_flag = ndviFlagCoding.addFlag("NDVI_NEG", 4, "NDVI negative")
        pos_flag = ndviFlagCoding.addFlag("NDVI_POS", 8, "NDVI positive")
        ndvi_product.getFlagCodingGroup().add(ndviFlagCoding)
        self.ndvi_flags_band.setSampleCoding(ndviFlagCoding)

        # Also for each flag a layer should be created
        ndvi_product.addMask('mask_' + low_flag.getName(), 'ndvi_flags.' + low_flag.getName(),
                             low_flag.getDescription(), Color.YELLOW, 0.3)
        ndvi_product.addMask('mask_' + high_flag.getName(), 'ndvi_flags.' + high_flag.getName(),
                             high_flag.getDescription(), Color.GREEN, 0.3)
        ndvi_product.addMask('mask_' + neg_flag.getName(), 'ndvi_flags.' + neg_flag.getName(),
                             neg_flag.getDescription(), Color(255, 0, 0), 0.3)
        ndvi_product.addMask('mask_' + pos_flag.getName(), 'ndvi_flags.' + pos_flag.getName(),
                             pos_flag.getDescription(), Color.BLUE, 0.3)

        # Provide the created target product to the framework so the computeTileStack method can be called
        # properly and the data can be written to disk.
        context.setTargetProduct(ndvi_product)
        f.write('end initialize.')
        f.close()

    def computeTileStack(self, context, target_tiles, target_rectangle):
        # The operator is asked by the framework to provide the data for a rectangle when the data is needed.
        # The required source data for the computation can be retrieved by getSourceTile(...) via the context object.
        lower_tile = context.getSourceTile(self.lower_band, target_rectangle)
        upper_tile = context.getSourceTile(self.upper_band, target_rectangle)

        # The actual data can be retrieved from the tiles by getSampleFloats(), getSamplesDouble() or getSamplesInt()
        lower_samples = lower_tile.getSamplesFloat()
        upper_samples = upper_tile.getSamplesFloat()
        # Values at specific pixel locations can be retrieved for example by lower_tile.getSampleFloat(x, y)

        # Convert the data into numpy data. It is easier and faster to work with as if you use plain python arithmetic
        lower_data = numpy.array(lower_samples, dtype=numpy.float32) * self.lower_factor
        upper_data = numpy.array(upper_samples, dtype=numpy.float32) * self.upper_factor

        # Doing the actual computation
        ndvi = self.sice2algo.compute_ndvi(lower_data, upper_data)
        ndvi_flags = self.sice2algo.compute_flags(ndvi)

        # The target tile which shall be filled with data are provided as parameter to this method
        ndvi_tile = target_tiles.get(self.ndvi_band)
        ndvi_flags_tile = target_tiles.get(self.ndvi_flags_band)

        # Set the result to the target tiles
        ndvi_tile.setSamples(ndvi)
        ndvi_flags_tile.setSamples(ndvi_flags)

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

