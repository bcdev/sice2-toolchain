import datetime
from math import ceil

import platform
import tempfile
import sys
import os
import time
import numpy as np

import xarray as xr

# Append esa_snappy installation dir to path:
sys.path.append(os.path.expanduser('~') + os.sep + '.snap' + os.sep + 'snap-python')

import esa_snappy
from esa_snappy import ProductData
from esa_snappy import ProductIO
from esa_snappy import FlagCoding

# If a Java type is needed which is not imported by snappy by default it can be retrieved manually.
# First import jpy
from esa_snappy import jpy

# and then import the type
Float = jpy.get_type('java.lang.Float')
Color = jpy.get_type('java.awt.Color')

BitSetter = jpy.get_type('org.esa.snap.core.util.BitSetter')
BandMathsType = jpy.get_type('org.esa.snap.core.datamodel.Mask$BandMathsType')
HashMap = jpy.get_type('java.util.HashMap')

from esa_snappy import Product
from esa_snappy import ProductUtils
from esa_snappy import GPF

import sice2_constants
import sice2_v21_algo
import sice2_v21_io
import sice2_v21_utils


class Sice2V21Op:
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

        # Via the context object the source product which shall be processed can be retrieved
        self.source_product = context.getSourceProduct('l1bProduct')
        self.cloud_product = context.getSourceProduct('cloudProduct')

        self.width = self.source_product.getSceneRasterWidth()
        self.height = self.source_product.getSceneRasterHeight()

        print('initialize: source product location is', self.source_product.getFileLocation())

        self.cloud_mask_band = None
        if self.cloud_product is not None:
            print('initialize: cloud product location is', self.cloud_product.getFileLocation())

            # if SCDA cloud product:
            if self.cloud_product.containsBand(sice2_constants.SCDA_FLAG_BAND_NAME):
                olci_rad_subset_product = Product('Rad subset', 'Rad subset', self.width, self.height)
                ProductUtils.copyGeoCoding(self.source_product, olci_rad_subset_product)
                ProductUtils.copyBand('Oa10_radiance', self.source_product, olci_rad_subset_product, True)

                # collocate the two subsets, OLCI as master:
                input_products = HashMap()
                input_products.put('master', olci_rad_subset_product)
                input_products.put('slave', self.cloud_product)
                parameters = HashMap()
                parameters.put('masterComponentPattern', '${ORIGINAL_NAME}')
                parameters.put('slaveComponentPattern', '${ORIGINAL_NAME}')
                operator_name = 'Collocate'
                scda_cloud_product_on_olci_raster = GPF.createProduct(operator_name, parameters, input_products)
                self.cloud_mask_band = self._get_band(scda_cloud_product_on_olci_raster, "scda_cloud_mask")
                self.cloud_product = scda_cloud_product_on_olci_raster
            elif self.cloud_product.containsBand(sice2_constants.IDEPIX_FLAG_BAND_NAME):
                # Idepix
                self.cloud_mask_band = self._get_band(self.cloud_product, sice2_constants.IDEPIX_FLAG_BAND_NAME)
                pass
            else:
                raise Exception('Selected cloud product is not valid.')

        self.write_spectral_albedos = context.getParameter('writeSpectralAlbedos')

        self.called_compute_tile_stack = 0

        #####
        self.input_products_all_band_names = []
        for i in range(len(self.source_product.getBands())):
            self.input_products_all_band_names.append(self.source_product.getBandAt(i).getName())
        for i in range(len(self.source_product.getTiePointGrids())):
            self.input_products_all_band_names.append(self.source_product.getTiePointGridAt(i).getName())
        if self.cloud_product is not None:
            for i in range(len(self.cloud_product.getBands())):
                self.input_products_all_band_names.append(self.cloud_product.getBandAt(i).getName())
            for i in range(len(self.cloud_product.getTiePointGrids())):
                self.input_products_all_band_names.append(self.cloud_product.getTiePointGridAt(i).getName())

        self.sza_band = self._get_band(self.source_product, "SZA")
        self.saa_band = self._get_band(self.source_product, "SAA")
        self.vza_band = self._get_band(self.source_product, "OZA")
        self.vaa_band = self._get_band(self.source_product, "OAA")
        self.total_ozone_band = self._get_band(self.source_product, "total_ozone")
        self.altitude_band = self._get_band(self.source_product, "altitude")

        self.latitude_band = self._get_band(self.source_product, "latitude")
        self.longitude_band = self._get_band(self.source_product, "longitude")

        self.radiance_bands = []
        self.radiance_band_names = [f'Oa{i:02}_radiance' for i in range(1, sice2_constants.OLCI_NUM_SPECTRAL_BANDS + 1)]
        for i in range(1, len(self.radiance_band_names) + 1):
            self.radiance_bands.append(self._get_band(self.source_product, self.radiance_band_names[i - 1]))

        self.solar_flux_bands = []
        self.solar_flux_band_names = [f'solar_flux_band_{i:01}' for i in
                                      range(1, sice2_constants.OLCI_NUM_SPECTRAL_BANDS + 1)]
        for i in range(1, len(self.solar_flux_band_names) + 1):
            self.solar_flux_bands.append(self._get_band(self.source_product, self.solar_flux_band_names[i - 1]))

        # Create the target product
        snow_product = esa_snappy.Product('py_SICE21_snow', 'py_SICE21_snow', self.width, self.height)
        self.pref_tile_width = 1000
        self.pref_tile_height = 1000
        snow_product.setPreferredTileSize(self.pref_tile_width, self.pref_tile_height)
        self.num_tiles_to_process = ceil(self.width / self.pref_tile_width) * ceil(self.height / self.pref_tile_height)
        esa_snappy.ProductUtils.copyGeoCoding(self.source_product, snow_product)
        esa_snappy.ProductUtils.copyMetadata(self.source_product, snow_product)
        snow_product.setStartTime(self.source_product.getStartTime())
        snow_product.setEndTime(self.source_product.getEndTime())
        snow_product.setAutoGrouping("albedo_spectral_spherical:albedo_spectral_planar:rBRR")

        if self.cloud_product is not None:
            if self.cloud_product.containsBand(sice2_constants.SCDA_FLAG_BAND_NAME):
                ProductUtils.copyBand(sice2_constants.SCDA_FLAG_BAND_NAME, self.cloud_product, snow_product, True)
                sice2_v21_io.create_scda_bitmask(snow_product)
            elif self.cloud_product.containsBand(sice2_constants.IDEPIX_FLAG_BAND_NAME):
                ProductUtils.copyBand(sice2_constants.IDEPIX_FLAG_BAND_NAME, self.cloud_product, snow_product, True)
                ProductUtils.copyFlagCoding(self.cloud_mask_band.getFlagCoding(), snow_product)
                sice2_v21_io.create_idepix_bitmask(snow_product)

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

        :param snow_product:
        :return:
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

        self.isnow_band = snow_product.addBand('isnow', esa_snappy.ProductData.TYPE_INT16)
        self.isnow_band.setDescription('Snow retrieval flag')
        snow_retrieval_flag_coding = sice2_v21_io.create_snow_retrieval_flag_coding(sice2_constants.SNOW_TYPE_FLAG)
        self.isnow_band.setSampleCoding(snow_retrieval_flag_coding)
        snow_product.getFlagCodingGroup().add(snow_retrieval_flag_coding)
        sice2_v21_io.create_snow_retrieval_bitmask(snow_product)

        self.pol_type_band = snow_product.addBand('pol_type', esa_snappy.ProductData.TYPE_INT16)
        self.pol_type_band.setDescription('Type of pollutant')
        pol_type_index_coding, image_info = sice2_v21_io.create_pol_type_index_coding(sice2_constants.POL_TYPE_FLAG)
        self.pol_type_band.setImageInfo(image_info)
        snow_product.getIndexCodingGroup().add(pol_type_index_coding)
        self.pol_type_band.setSampleCoding(pol_type_index_coding)

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

        self.impurity_load_band = snow_product.addBand('impurity_load', esa_snappy.ProductData.TYPE_FLOAT32)
        self.impurity_load_band.setDescription('Pollutant load')
        self.impurity_load_band.setNoDataValue(Float.NaN)
        self.impurity_load_band.setNoDataValueUsed(True)

        if self.write_spectral_albedos:
            self.albedo_spectral_spherical_bands = []
            self.albedo_spectral_planar_bands = []
            self.r_brr_bands = []
            for i in np.append(np.arange(11), np.arange(15, sice2_constants.OLCI_NUM_SPECTRAL_BANDS)):
                albedo_spectral_spherical_band = snow_product.addBand(
                    'albedo_spectral_spherical_' + str(i + 1).zfill(2),
                    esa_snappy.ProductData.TYPE_FLOAT32)
                albedo_spectral_spherical_band.setDescription('Spectral spherical albedo')
                albedo_spectral_spherical_band.setNoDataValue(Float.NaN)
                albedo_spectral_spherical_band.setNoDataValueUsed(True)
                self.albedo_spectral_spherical_bands.append(albedo_spectral_spherical_band)
                albedo_spectral_planar_band = snow_product.addBand('albedo_spectral_planar_' + str(i + 1).zfill(2),
                                                                   esa_snappy.ProductData.TYPE_FLOAT32)
                albedo_spectral_planar_band.setDescription('Spectral planar albedo')
                albedo_spectral_planar_band.setNoDataValue(Float.NaN)
                albedo_spectral_planar_band.setNoDataValueUsed(True)
                self.albedo_spectral_planar_bands.append(albedo_spectral_planar_band)
                r_brr_band = snow_product.addBand('rBRR_' + str(i + 1).zfill(2),
                                                  esa_snappy.ProductData.TYPE_FLOAT32)
                r_brr_band.setDescription('Bottom-of-Rayleigh reflectance')
                r_brr_band.setNoDataValue(Float.NaN)
                r_brr_band.setNoDataValueUsed(True)
                self.r_brr_bands.append(r_brr_band)

    def computeTileStack_test(self, context, target_tiles, target_rectangle):
        pass

    def computeTileStack(self, context, target_tiles, target_rectangle):

        num_pixels = target_rectangle.width * target_rectangle.height
        print('Call computeTileStack: num_pixels=' + str(num_pixels))
        print('target_rectangle.x =' + str(target_rectangle.x))
        print('target_rectangle.y =' + str(target_rectangle.y))
        print('target_rectangle.width =' + str(target_rectangle.width))
        print('target_rectangle.height =' + str(target_rectangle.height))
        self.called_compute_tile_stack = self.called_compute_tile_stack + 1
        print('Tile ' + str(self.called_compute_tile_stack) + ' of ' + str(self.num_tiles_to_process))

        # filter invalid (i.e. cloudy) pixels:
        valid_expression_filter_array = self.setup_valid_expression_filter(context, target_rectangle)

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
            rad_tile = context.getSourceTile(self.radiance_bands[i], target_rectangle)
            rad_data = np.array(rad_tile.getSamplesFloat(), dtype=np.float32)
            solar_flux_tile = context.getSourceTile(self.solar_flux_bands[i], target_rectangle)
            solar_flux_data = np.array(solar_flux_tile.getSamplesFloat(), dtype=np.float32)

            _rad = self._get_var(rad_data, target_rectangle.width, target_rectangle.height, bands[i], 'dl')
            _toa = (_rad * np.pi) / (solar_flux_data * np.cos(np.deg2rad(sza_data)))
            _toa[np.where(~valid_expression_filter_array)] = np.nan
            toa.append(np.clip(_toa, 0, 1))
        olci_scene['toa'] = xr.concat(toa, dim='band')

        ##### SNOW RETRIEVAL:
        chunk_size = int(min(num_pixels, 250000))
        print('Call process_by_chunk: chunksize=' + str(chunk_size))
        snow = sice2_v21_algo.process_by_chunk(olci_scene, chunk_size=chunk_size)
        #####

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

        if self.write_spectral_albedos:
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
        target_tiles.get(self.albedo_bb_planar_sw_band).setSamples(albedo_bb_planar_sw_data)
        target_tiles.get(self.albedo_bb_spherical_sw_band).setSamples(albedo_bb_spherical_sw_data)
        target_tiles.get(self.factor_band).setSamples(factor_data)
        target_tiles.get(self.o3_sice_band).setSamples(o3_sice_data)
        target_tiles.get(self.cv1_band).setSamples(cv1_data)
        target_tiles.get(self.cv2_band).setSamples(cv2_data)
        target_tiles.get(self.impurity_load_band).setSamples(impurity_load_data)
        target_tiles.get(self.pol_type_band).setSamples(pol_type_data)

        for key in sice2_constants.snow_type_flags_map:
            isnow_data[np.where(isnow_data == float(key))] = sice2_constants.snow_type_flags_map[key]
        target_tiles.get(self.isnow_band).setSamples(isnow_data)

        if self.write_spectral_albedos:
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
        """

        :param data:
        :param width:
        :param height:
        :param long_name:
        :param unit:
        :return:
        """
        data_da = xr.DataArray(
            data=data.reshape(width, height),
            dims=["x", "y"],
            attrs={'long_name': long_name,
                   'unit': unit}
        )
        return data_da.compute().stack(xy=("x", "y"))

    def setup_valid_expression_filter(self, context, rect):
        """

        :param context:
        :param rect:
        :return:
        """
        if self.cloud_product is not None:
            if self.cloud_product.containsBand(sice2_constants.SCDA_FLAG_BAND_NAME):
                # SCDA
                flagname = sice2_constants.SCDA_FLAG_BAND_NAME
                expr = sice2_constants.DEFAULT_SCDA_VALID_PIXEL_EXPR
                flag_coding_dict = sice2_constants.SCDA_BITMASK_FLAG_CODING_DICT
            else:
                flagname = sice2_constants.IDEPIX_FLAG_BAND_NAME
                expr = sice2_constants.DEFAULT_IDEPIX_VALID_PIXEL_EXPR
                flag_coding_dict = sice2_constants.IDEPIX_BITMASK_FLAG_CODING_DICT

            print('flagname: ' + flagname)
            print('expr: ' + expr)
            condition = sice2_v21_utils.get_condition_from_valid_pixel_expr(flagname, expr, flag_coding_dict)
            print('condition: ' + condition)
            # identify which of the product variables are in the valid expr...
            variables_in_expr_dict = {}
            for name in self.input_products_all_band_names:
                if name in expr:
                    if self.source_product.containsBand(name):
                        _band = self.source_product.getBand(name)
                    else:
                        _band = self.cloud_product.getBand(name)
                    _tile = context.getSourceTile(_band, rect)
                    if ProductData.isIntType(_band.getDataType()):
                        variables_in_expr_dict[name] = np.array(_tile.getSamplesInt(), dtype=np.int32)
                    elif ProductData.isFloatingPointType(_band.getDataType()):
                        variables_in_expr_dict[name] = np.array(_tile.getSamplesFloat(), dtype=np.float32)
                    variables_in_expr_dict[name] = np.reshape(variables_in_expr_dict[name], (rect.width, rect.height))

            return sice2_v21_utils.get_valid_expression_filter_array(condition, variables_in_expr_dict,
                                                                     rect.width, rect.height).flatten()
        else:
            # no filtering
            return np.full((rect.width * rect.height), True)