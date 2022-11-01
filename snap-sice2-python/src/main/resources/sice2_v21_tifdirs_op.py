import datetime
from math import ceil

import sice2_constants
import sice2_v21_algo

import platform
import tempfile
import sys
import os
import time
import numpy as np

import xarray as xr

import esa_snappy
from esa_snappy import ProductIO
from esa_snappy import RsMathUtils

# If a Java type is needed which is not imported by snappy by default it can be retrieved manually.
# First import jpy
from esa_snappy import jpy

# and then import the type
import sice2_constants
import sice2_utils
import sice2_v21_utils

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

        self.called_compute_tile_stack = 0

        tif_source_product_paths = sice2_v21_utils.Sice2V21Utils.get_tif_source_product_paths(self.tif_input_directory)

        sza_product = ProductIO.readProduct(self.tif_input_directory + os.sep + "SZA.tif")
        saa_product = ProductIO.readProduct(self.tif_input_directory + os.sep + "SAA.tif")
        vza_product = ProductIO.readProduct(self.tif_input_directory + os.sep + "OZA.tif")
        vaa_product = ProductIO.readProduct(self.tif_input_directory + os.sep + "OAA.tif")
        total_ozone_product = ProductIO.readProduct(self.tif_input_directory + os.sep + "O3.tif")
        altitude_product = ProductIO.readProduct(self.tif_input_directory + os.sep + "height.tif")

        self.width = sza_product.getSceneRasterWidth()
        self.height = sza_product.getSceneRasterHeight()

        #####
        self.sza_band = self._get_band(sza_product, "band_1")
        self.saa_band = self._get_band(saa_product, "band_1")
        self.vza_band = self._get_band(vza_product, "band_1")
        self.vaa_band = self._get_band(vaa_product, "band_1")
        self.total_ozone_band = self._get_band(total_ozone_product, "band_1")
        self.altitude_band = self._get_band(altitude_product, "band_1")

        self.rtoa_bands = []
        self.rtoa_file_names = [f'r_TOA_{i:02}' for i in range(1, sice2_constants.OLCI_NUM_SPECTRAL_BANDS + 1)]
        for i in range(1, len(self.rtoa_file_names) + 1):
            rtoa_path = self.tif_input_directory + os.sep + self.rtoa_file_names[i-1] + '.tif'
            if rtoa_path in tif_source_product_paths:
                rtoa_product = ProductIO.readProduct(rtoa_path)
                self.rtoa_bands.append(self._get_band(rtoa_product, "band_1"))
            else:
                self.rtoa_bands.append(None)

        # Create the target product
        snow_product = esa_snappy.Product('py_SICE21_snow', 'py_SICE21_snow', self.width, self.height)
        # snow_product.setPreferredTileSize(self.width, self.height)
        self.pref_tile_width = 1000
        self.pref_tile_height = 1000
        snow_product.setPreferredTileSize(self.pref_tile_width, self.pref_tile_height)
        self.num_tiles_to_process = ceil(self.width/self.pref_tile_width) * ceil(self.height/self.pref_tile_height)
        # self.tiles_processed = []
        # self.tiles_processed_list_cleared = False
        esa_snappy.ProductUtils.copyGeoCoding(sza_product, snow_product)
        esa_snappy.ProductUtils.copyMetadata(sza_product, snow_product)
        snow_product.setStartTime(sza_product.getStartTime())
        snow_product.setEndTime(sza_product.getEndTime())

        context.setTargetProduct(snow_product)

        # grain_diameter_data = snow['diameter'].values
        # snow_specific_area_data = snow['area'].values
        # al_data = snow['al'].values
        # r0_data = snow['r0'].values
        # isnow_data = snow['isnow'].values
        # conc_data = snow['conc'].values
        # albedo_bb_planar_sw_data = snow['rp3'].values
        # albedo_bb_spherical_sw_data = snow['rs3'].values
        # factor_data = snow['factor'].values

        self.grain_diameter_band = snow_product.addBand('grain_diameter', esa_snappy.ProductData.TYPE_FLOAT32)
        self.grain_diameter_band.setDescription('grain_diameter')
        self.grain_diameter_band.setNoDataValue(Float.NaN)
        self.grain_diameter_band.setNoDataValueUsed(True)

        self.snow_specific_area_band = snow_product.addBand('snow_specific_area', esa_snappy.ProductData.TYPE_FLOAT32)
        self.snow_specific_area_band.setDescription('snow_specific_area')
        self.snow_specific_area_band.setNoDataValue(Float.NaN)
        self.snow_specific_area_band.setNoDataValueUsed(True)

        self.al_band = snow_product.addBand('al', esa_snappy.ProductData.TYPE_FLOAT32)
        self.al_band.setDescription('al')
        self.al_band.setNoDataValue(Float.NaN)
        self.al_band.setNoDataValueUsed(True)

        self.r0_band = snow_product.addBand('r0', esa_snappy.ProductData.TYPE_FLOAT32)
        self.r0_band.setDescription('r0')
        self.r0_band.setNoDataValue(Float.NaN)
        self.r0_band.setNoDataValueUsed(True)

        self.isnow_band = snow_product.addBand('isnow', esa_snappy.ProductData.TYPE_FLOAT32)
        self.isnow_band.setDescription('isnow')
        self.isnow_band.setNoDataValue(Float.NaN)
        self.isnow_band.setNoDataValueUsed(True)

        self.albedo_bb_planar_sw_band = snow_product.addBand('albedo_bb_planar_sw', esa_snappy.ProductData.TYPE_FLOAT32)
        self.albedo_bb_planar_sw_band.setDescription('Planar SW albedo')
        self.albedo_bb_planar_sw_band.setNoDataValue(Float.NaN)
        self.albedo_bb_planar_sw_band.setNoDataValueUsed(True)

        self.albedo_bb_spherical_sw_band = snow_product.addBand('albedo_bb_spherical_sw', esa_snappy.ProductData.TYPE_FLOAT32)
        self.albedo_bb_spherical_sw_band.setDescription('Spherical SW albedo')
        self.albedo_bb_spherical_sw_band.setNoDataValue(Float.NaN)
        self.albedo_bb_spherical_sw_band.setNoDataValueUsed(True)

        self.factor_band = snow_product.addBand('factor', esa_snappy.ProductData.TYPE_FLOAT32)
        self.factor_band.setDescription('factor')
        self.factor_band.setNoDataValue(Float.NaN)
        self.factor_band.setNoDataValueUsed(True)

        f.write('end initialize.')
        print('end initialize.')
        t_end = datetime.datetime.now()
        print('SNAPPY SICE2 processing time (seconds): ' + str((t_end - t_start).seconds))

        end_time = time.process_time()
        print('SNAPPY SICE2 processing time (CPU seconds): ' + str((end_time - start_time)))

        f.close()

    def computeTileStack(self, context, target_tiles, target_rectangle):

        num_pixels = target_rectangle.width * target_rectangle.height
        print('Call computeTileStack: num_pixels=' + str(num_pixels))
        print('target_rectangle.x =' + str(target_rectangle.x))
        print('target_rectangle.y =' + str(target_rectangle.y))
        print('target_rectangle.width =' + str(target_rectangle.width))
        print('target_rectangle.height =' + str(target_rectangle.height))
        self.called_compute_tile_stack = self.called_compute_tile_stack + 1
        print('Tile ' + str(self.called_compute_tile_stack) + ' of ' + str(self.num_tiles_to_process))

        sza_tile = context.getSourceTile(self.sza_band, target_rectangle)
        saa_tile = context.getSourceTile(self.saa_band, target_rectangle)
        vza_tile = context.getSourceTile(self.vza_band, target_rectangle)
        vaa_tile = context.getSourceTile(self.vaa_band, target_rectangle)
        total_ozone_tile = context.getSourceTile(self.total_ozone_band, target_rectangle)
        altitude_tile = context.getSourceTile(self.altitude_band, target_rectangle)

        sza_data = np.array(sza_tile.getSamplesFloat(), dtype=np.float32)
        saa_data = np.array(saa_tile.getSamplesFloat(), dtype=np.float32)
        vza_data = np.array(vza_tile.getSamplesFloat(), dtype=np.float32)
        vaa_data = np.array(vaa_tile.getSamplesFloat(), dtype=np.float32)
        total_ozone_data = np.array(total_ozone_tile.getSamplesFloat(), dtype=np.float32)
        altitude_data = np.array(altitude_tile.getSamplesFloat(), dtype=np.float32)

        variables = {
            'solar_zenith_angle': ['sza', 'deg', sza_data],
            'solar_azimuth_angle': ['saa', 'deg', saa_data],
            'satellite_zenith_angle': ['vza','deg', vza_data],
            'satellite_azimuth_angle': ['vaa','deg', vaa_data],
            'total_ozone': ['ozone','DU', total_ozone_data],
            'altitude': ['elevation', 'm', altitude_data]
        }

        olci_scene = xr.Dataset()

        for variable in variables:
            var_name = variables[variable][0]
            var_unit = variables[variable][1]
            var_data = variables[variable][2]
            olci_scene[var_name] = self._get_var(var_data, target_rectangle.width, target_rectangle.height, var_data, var_unit)

        bands = [f'Oa{i:02}' for i in range(1, sice2_constants.OLCI_NUM_SPECTRAL_BANDS + 1)]

        # dat = xr.full_like(self.toa[0], fill_value=np.nan)
        toa = []
        for i in range(len(self.rtoa_bands)):
            if self.rtoa_bands[i] is not None:
                rtoa_tile = context.getSourceTile(self.rtoa_bands[i], target_rectangle)
                rtoa_data = np.array(rtoa_tile.getSamplesFloat(), dtype=np.float32)
            else:
                rtoa_data = np.empty(target_rectangle.width * target_rectangle.height)
                rtoa_data[:] = np.nan
            toa.append(self._get_var(rtoa_data, target_rectangle.width, target_rectangle.height, bands[i], 'dl'))

        olci_scene['toa'] = xr.concat(toa, dim='band')

        # snow = sice2_v21_algo.process(olci_scene)
        if num_pixels == 1000000:
            chunk_size = 250000
        else:
            chunk_size = num_pixels
        print('Call process_by_chunk: chunksize=' + str(num_pixels))
        snow = sice2_v21_algo.process_by_chunk(olci_scene, chunk_size=chunk_size)

        grain_diameter_data = snow['diameter'].values
        snow_specific_area_data = snow['area'].values
        al_data = snow['al'].values
        r0_data = snow['r0'].values
        isnow_data = snow['isnow'].values
        # conc_data = snow['conc'].values
        albedo_bb_planar_sw_data = snow['rp3'].values
        albedo_bb_spherical_sw_data = snow['rs3'].values
        factor_data = snow['factor'].values

        # The target tile which shall be filled with data are provided as parameter to this method
        grain_diameter_tile = target_tiles.get(self.grain_diameter_band)
        snow_specific_area_tile = target_tiles.get(self.snow_specific_area_band)
        al_tile = target_tiles.get(self.al_band)
        r0_tile = target_tiles.get(self.r0_band)
        isnow_tile = target_tiles.get(self.isnow_band)
        albedo_bb_planar_sw_tile = target_tiles.get(self.albedo_bb_planar_sw_band)
        albedo_bb_spherical_sw_tile = target_tiles.get(self.albedo_bb_spherical_sw_band)
        factor_tile = target_tiles.get(self.factor_band)

        # Set the result to the target tiles
        grain_diameter_tile.setSamples(grain_diameter_data)
        snow_specific_area_tile.setSamples(snow_specific_area_data)
        al_tile.setSamples(al_data)
        r0_tile.setSamples(r0_data)
        isnow_tile.setSamples(isnow_data)
        albedo_bb_planar_sw_tile.setSamples(albedo_bb_planar_sw_data)
        albedo_bb_spherical_sw_tile.setSamples(albedo_bb_spherical_sw_data)
        factor_tile.setSamples(factor_data)


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


    def _get_var(self, data, width, height, long_name, unit):
        data_da = xr.DataArray(
            data=data.reshape(width, height),
            dims=["x", "y"],
            attrs={'long_name': long_name,
                   'unit': unit}
        )
        return data_da.compute().stack(xy=("x", "y"))