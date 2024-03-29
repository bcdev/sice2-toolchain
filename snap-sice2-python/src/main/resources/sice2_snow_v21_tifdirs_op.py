import datetime
import os
import platform
import sys
import tempfile
import time
from math import ceil

import esa_snappy
import numpy as np
import xarray as xr
from esa_snappy import ProductIO
# If a Java type is needed which is not imported by snappy by default it can be retrieved manually.
# First import jpy
from esa_snappy import jpy

# and then import the type
import sice2_constants
import sice2_snow_v21_algo
import sice2_utils

Float = jpy.get_type('java.lang.Float')
Color = jpy.get_type('java.awt.Color')


class Sice2SnowV21TifdirsOp:
    """
    SICE2 operator for the retrieval of snow properties. (Uses OLCI single TIFs, GEUS toolchain mode.)

    @author: Olaf Danne, BC (Brockmann Consult)
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
        f = open(tempfile.gettempdir() + '/sice2_py.log', 'w')

        sys.path.append(resource_root)

        f.write('Python module location: ' + __file__ + '\n')
        f.write('Python module location parent: ' + resource_root + '\n')

        print('platform.system(): ' + platform.system() + '\n')
        print('sys.version_info(): ' + str(sys.version_info) + '\n')
        print('sys.version_info > (3, 0): ' + str(sys.version_info > (3, 0)))

        # Tif directory parameter defined in ndsi_op-info.xml
        tif_input_directory = context.getParameter('tifInputDir')

        print('Start sice.py')
        print('Input folder:', tif_input_directory)

        self.called_compute_tile_stack = 0

        tif_source_product_paths = sice2_utils.get_tif_source_product_paths(tif_input_directory)

        sza_product = ProductIO.readProduct(tif_input_directory + os.sep + "SZA.tif")
        saa_product = ProductIO.readProduct(tif_input_directory + os.sep + "SAA.tif")
        vza_product = ProductIO.readProduct(tif_input_directory + os.sep + "OZA.tif")
        vaa_product = ProductIO.readProduct(tif_input_directory + os.sep + "OAA.tif")
        total_ozone_product = ProductIO.readProduct(tif_input_directory + os.sep + "O3.tif")
        altitude_product = ProductIO.readProduct(tif_input_directory + os.sep + "height.tif")

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
            rtoa_path = tif_input_directory + os.sep + self.rtoa_file_names[i - 1] + '.tif'
            if rtoa_path in tif_source_product_paths:
                rtoa_product = ProductIO.readProduct(rtoa_path)
                self.rtoa_bands.append(self._get_band(rtoa_product, "band_1"))
            else:
                self.rtoa_bands.append(None)

        # Create the target product
        snow_product = esa_snappy.Product('py_SICE21_snow', 'py_SICE21_snow', self.width, self.height)
        # snow_product.setPreferredTileSize(self.width, self.height)
        pref_tile_width = 1000
        pref_tile_height = 1000
        snow_product.setPreferredTileSize(pref_tile_width, pref_tile_height)
        self.num_tiles_to_process = ceil(self.width / pref_tile_width) * ceil(self.height / pref_tile_height)
        esa_snappy.ProductUtils.copyGeoCoding(sza_product, snow_product)
        esa_snappy.ProductUtils.copyMetadata(sza_product, snow_product)
        snow_product.setStartTime(sza_product.getStartTime())
        snow_product.setEndTime(sza_product.getEndTime())
        snow_product.setAutoGrouping("albedo_spectral_spherical:albedo_spectral_planar:rBRR")

        context.setTargetProduct(snow_product)

        self.add_target_bands(snow_product)

        f.write('end initialize.')
        print('end initialize.')
        t_end = datetime.datetime.now()
        print('SNAPPY SICE2 processing time (seconds): ' + str((t_end - t_start).seconds))

        end_time = time.process_time()
        print('SNAPPY SICE2 processing time (CPU seconds): ' + str((end_time - start_time)))

        f.close()

    def add_target_bands(self, snow_product):
        """
         Adds bands to snow target product.

         :param snow_product: snow target product
         :return: void
         """
        self.grain_diameter_band = snow_product.addBand('grain_diameter', esa_snappy.ProductData.TYPE_FLOAT32)
        self.grain_diameter_band.setDescription('Snow grain diameter')
        self.grain_diameter_band.setNoDataValue(Float.NaN)
        self.grain_diameter_band.setNoDataValueUsed(True)

        self.snow_specific_area_band = snow_product.addBand('snow_specific_area', esa_snappy.ProductData.TYPE_FLOAT32)
        self.snow_specific_area_band.setDescription('Snow specific surface area')
        self.snow_specific_area_band.setNoDataValue(Float.NaN)
        self.snow_specific_area_band.setNoDataValueUsed(True)

        self.al_band = snow_product.addBand('al', esa_snappy.ProductData.TYPE_FLOAT32)
        self.al_band.setDescription('Effective absorption length')
        self.al_band.setNoDataValue(Float.NaN)
        self.al_band.setNoDataValueUsed(True)

        self.r0_band = snow_product.addBand('r0', esa_snappy.ProductData.TYPE_FLOAT32)
        self.r0_band.setDescription('Reflectance of a semi-infinite non-absorbing snow layer')
        self.r0_band.setNoDataValue(Float.NaN)
        self.r0_band.setNoDataValueUsed(True)

        self.isnow_band = snow_product.addBand('isnow', esa_snappy.ProductData.TYPE_FLOAT32)
        self.isnow_band.setDescription('Snow retrieval flag')
        self.isnow_band.setNoDataValue(Float.NaN)
        self.isnow_band.setNoDataValueUsed(True)

        self.albedo_bb_planar_sw_band = snow_product.addBand('albedo_bb_planar_sw', esa_snappy.ProductData.TYPE_FLOAT32)
        self.albedo_bb_planar_sw_band.setDescription('Shortwave broadband planar albedo')
        self.albedo_bb_planar_sw_band.setNoDataValue(Float.NaN)
        self.albedo_bb_planar_sw_band.setNoDataValueUsed(True)

        self.albedo_bb_spherical_sw_band = snow_product.addBand('albedo_bb_spherical_sw',
                                                                esa_snappy.ProductData.TYPE_FLOAT32)
        self.albedo_bb_spherical_sw_band.setDescription('Shortwave broadband spherical albedo')
        self.albedo_bb_spherical_sw_band.setNoDataValue(Float.NaN)
        self.albedo_bb_spherical_sw_band.setNoDataValueUsed(True)

        self.factor_band = snow_product.addBand('factor', esa_snappy.ProductData.TYPE_FLOAT32)
        self.factor_band.setDescription('Snow covered fraction within mixed pixels')
        self.factor_band.setNoDataValue(Float.NaN)
        self.factor_band.setNoDataValueUsed(True)

        self.o3_sice_band = snow_product.addBand('O3_SICE', esa_snappy.ProductData.TYPE_FLOAT32)
        self.o3_sice_band.setDescription('Total ozone product (OLCI) corrected for ozone scattering')
        self.o3_sice_band.setNoDataValue(Float.NaN)
        self.o3_sice_band.setNoDataValueUsed(True)

        self.cv1_band = snow_product.addBand('cv1', esa_snappy.ProductData.TYPE_FLOAT32)
        self.cv1_band.setDescription('Quality check 1')
        self.cv1_band.setNoDataValue(Float.NaN)
        self.cv1_band.setNoDataValueUsed(True)

        self.cv2_band = snow_product.addBand('cv2', esa_snappy.ProductData.TYPE_FLOAT32)
        self.cv2_band.setDescription('Quality check 2')
        self.cv2_band.setNoDataValue(Float.NaN)
        self.cv2_band.setNoDataValueUsed(True)

        self.pol_type_band = snow_product.addBand('pol_type', esa_snappy.ProductData.TYPE_FLOAT32)
        self.pol_type_band.setDescription('Type of pollutant: 1(soot), 2( dust), 3 and 4 (other or mixture)')
        self.pol_type_band.setNoDataValue(Float.NaN)
        self.pol_type_band.setNoDataValueUsed(True)

        self.impurity_load_band = snow_product.addBand('impurity_load', esa_snappy.ProductData.TYPE_FLOAT32)
        self.impurity_load_band.setDescription('Pollutant load')
        self.impurity_load_band.setNoDataValue(Float.NaN)
        self.impurity_load_band.setNoDataValueUsed(True)

        self.albedo_spectral_spherical_bands = []
        self.albedo_spectral_planar_bands = []
        self.r_brr_bands = []
        for i in np.append(np.arange(11), np.arange(15, sice2_constants.OLCI_NUM_SPECTRAL_BANDS)):
            albedo_spectral_spherical_band = snow_product.addBand('albedo_spectral_spherical_' + str(i + 1).zfill(2),
                                                                  esa_snappy.ProductData.TYPE_FLOAT32)
            self.albedo_spectral_spherical_bands.append(albedo_spectral_spherical_band)
            albedo_spectral_planar_band = snow_product.addBand('albedo_spectral_planar_' + str(i + 1).zfill(2),
                                                               esa_snappy.ProductData.TYPE_FLOAT32)
            self.albedo_spectral_planar_bands.append(albedo_spectral_planar_band)
            r_brr_band = snow_product.addBand('rBRR_' + str(i + 1).zfill(2),
                                              esa_snappy.ProductData.TYPE_FLOAT32)
            self.r_brr_bands.append(r_brr_band)

    def computeTileStack(self, context, target_tiles, target_rectangle):
        """
        The GPF computeTileStack implementation.

        :param context: operator context
        :param target_tiles: target tiles
        :param target_rectangle: target rectangle
        :return: void
        """

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
            'satellite_zenith_angle': ['vza', 'deg', vza_data],
            'satellite_azimuth_angle': ['vaa', 'deg', vaa_data],
            'total_ozone': ['ozone', 'DU', total_ozone_data],
            'altitude': ['elevation', 'm', altitude_data]
        }

        olci_scene = xr.Dataset()

        for variable in variables:
            var_name = variables[variable][0]
            var_unit = variables[variable][1]
            var_data = variables[variable][2]
            olci_scene[var_name] = self._get_var(var_data, target_rectangle.width, target_rectangle.height, var_data,
                                                 var_unit)

        bands = [f'Oa{i:02}' for i in range(1, sice2_constants.OLCI_NUM_SPECTRAL_BANDS + 1)]

        toa = []
        for i in range(sice2_constants.OLCI_NUM_SPECTRAL_BANDS):
            if self.rtoa_bands[i] is not None:
                rtoa_tile = context.getSourceTile(self.rtoa_bands[i], target_rectangle)
                rtoa_data = np.array(rtoa_tile.getSamplesFloat(), dtype=np.float32)
            else:
                rtoa_data = np.empty(target_rectangle.width * target_rectangle.height)
                rtoa_data[:] = np.nan
            toa.append(self._get_var(rtoa_data, target_rectangle.width, target_rectangle.height, bands[i], 'dl'))

        olci_scene['toa'] = xr.concat(toa, dim='band')

        # ### SNOW RETRIEVAL:
        if num_pixels == 1000000:
            chunk_size = 250000
        else:
            chunk_size = num_pixels

        print('Call process_by_chunk: chunksize=' + str(num_pixels))
        snow = sice2_snow_v21_algo.process_by_chunk(olci_scene, chunk_size=chunk_size)
        ###

        # Extract output from 'snow' xarray.Dataset:
        grain_diameter_data = snow['diameter'].values
        snow_specific_area_data = snow['area'].values
        albedo_bb_planar_sw_data = snow['rp3'].values
        albedo_bb_spherical_sw_data = snow['rs3'].values
        isnow_data = snow['isnow'].values
        r0_data = snow['r0'].values
        al_data = snow['al'].values
        factor_data = snow['factor'].values
        o3_sice_data = snow['tocos'].values
        cv1_data = snow['cv1'].values
        cv2_data = snow['cv2'].values
        pol_type_data = snow['ntype'].values
        impurity_load_data = snow['aload_ppm'].values

        albedo_spectral_spherical_data = []
        albedo_spectral_planar_data = []
        r_brr_data = []
        for i in np.append(np.arange(11), np.arange(15, sice2_constants.OLCI_NUM_SPECTRAL_BANDS)):
            albedo_spectral_spherical_data.append(snow.alb_sph.sel(band=i).values)
            albedo_spectral_planar_data.append(snow.rp.sel(band=i).values)
            r_brr_data.append(snow.refl.sel(band=i).values)

        # The target tiles which shall be filled with data are provided as parameter to this method
        # Set the results to the target tiles
        target_tiles.get(self.grain_diameter_band).setSamples(grain_diameter_data)
        target_tiles.get(self.snow_specific_area_band).setSamples(snow_specific_area_data)
        target_tiles.get(self.al_band).setSamples(al_data)
        target_tiles.get(self.r0_band).setSamples(r0_data)
        target_tiles.get(self.isnow_band).setSamples(isnow_data)
        target_tiles.get(self.albedo_bb_planar_sw_band).setSamples(albedo_bb_planar_sw_data)
        target_tiles.get(self.albedo_bb_spherical_sw_band).setSamples(albedo_bb_spherical_sw_data)
        target_tiles.get(self.factor_band).setSamples(factor_data)
        target_tiles.get(self.o3_sice_band).setSamples(o3_sice_data)
        target_tiles.get(self.cv1_band).setSamples(cv1_data)
        target_tiles.get(self.cv2_band).setSamples(cv2_data)
        target_tiles.get(self.pol_type_band).setSamples(pol_type_data)
        target_tiles.get(self.impurity_load_band).setSamples(impurity_load_data)

        for i in range(sice2_constants.OLCI_NUM_SPECTRAL_BANDS - 4):
            target_tiles.get(self.albedo_spectral_spherical_bands[i]).setSamples(albedo_spectral_spherical_data[i])
            target_tiles.get(self.albedo_spectral_planar_bands[i]).setSamples(albedo_spectral_planar_data[i])
            target_tiles.get(self.r_brr_bands[i]).setSamples(r_brr_data[i])

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

    @staticmethod
    def _get_var(data, width, height, long_name, unit):
        data_da = xr.DataArray(
            data=data.reshape(width, height),
            dims=["x", "y"],
            attrs={'long_name': long_name,
                   'unit': unit}
        )
        return data_da.compute().stack(xy=("x", "y"))
