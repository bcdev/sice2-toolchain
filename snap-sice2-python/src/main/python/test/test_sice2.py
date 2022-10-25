import os
import sys
import unittest
import numpy as np
import xarray as xr
from sys import platform
from datetime import datetime

import sice2_v21_io
import sice2_v21_utils

# noinspection PyUnresolvedReferences
class TestSice2(unittest.TestCase):
    # def setUp(self):
    #     print('Platform: ', platform)
    #     # parent_dir = os.path.dirname(os.path.normpath(os.path.dirname(__file__)))
    #
    #     resource_root = os.path.dirname(__file__)
    #
    #     sys.path.append(resource_root)
    #
    #     SICE2_HOME = os.path.dirname(os.path.abspath(__file__))
    #     sys.path.append(SICE2_HOME)
    #
    #     import sice2_algo
    #     self.sice2algo = sice2_algo.Sice2Algo(-0.2, 0.8)

    # @unittest.skip("skipping test...")
    # def test_sice2algo(self):
    #     print('*** Compute NDSI ***')
    #     upper_data = 0.2
    #     lower_data = 0.1
    #
    #     ndsi = self.sice2algo.compute_ndsi(lower_data, upper_data)
    #
    #     print('NDSI: ' + str(ndsi))

    def test_sice2_io_nparrays_to_olci_scene_ds(self):
        print('*** NPARRAYS to OLCI_SCENE_DS ***')

        # 2*3 test numpy arrays...
        width = 2
        height = 3
        # width = 4864
        # height = 4085

        sza_data = np.random.rand(width*height) * 100.
        saa_data = np.random.rand(width*height) + 20.
        vaa_data = np.random.rand(width*height) + 100.
        vza_data = np.random.rand(width*height) + 50.
        ozone_data = np.random.rand(width*height) + 300.
        elevation_data = np.random.rand(width*height) + 1000.
        lon_data = np.random.rand(width*height) + 120.
        lat_data = np.random.rand(width*height) + 60.

        rad_data = np.empty([21, width*height])
        for i in range(rad_data.shape[0]):
            rad_data[i] = np.random.rand(width*height) * (i+1)

        variables = {
            'solar_zenith_angle': ['sza', 'deg', sza_data],
            'solar_azimuth_angle': ['saa', 'deg', saa_data],
            'satellite_azimuth_angle': ['vaa','deg', vaa_data],
            'satellite_zenith_angle': ['vza','deg', vza_data],
            'total_ozone': ['ozone','DU', ozone_data],
            'altitude': ['elevation', 'm', elevation_data]
        }

        olci_scene = xr.Dataset()

        for variable in variables:
            var_name = variables[variable][0]
            var_unit = variables[variable][1]
            var_data = variables[variable][2]
            olci_scene[var_name] = self.get_var(var_data, width, height, var_data, var_unit)

        olci_scene = olci_scene.assign_coords(
            longitude=self.get_var(lon_data, width, height, 'longitude', 'deg'),
            latitude=self.get_var(lat_data, width, height, 'latitude', 'deg'))

        coef = 1 / np.cos(np.deg2rad(olci_scene['sza'])) / 100.
        bands = [f'Oa{i:02}' for i in range(1, rad_data.shape[0] + 1)]

        toa = []
        for i in range(rad_data.shape[0]):
            toa.append(np.clip(self.get_var(rad_data[i], width, height, bands[i], 'dl') * coef, 0, 1))
        olci_scene['toa'] = xr.concat(toa, dim='band')

        print('done')

    def get_var(self, data, width, height, long_name, unit):
        data_da = xr.DataArray(
            data=data.reshape(width, height),
            dims=["x", "y"],
            attrs={'long_name': long_name,
                   'unit': unit}
        )
        return data_da.compute().stack(xy=("x", "y"))


    def test_check_if_tile_was_processed(self):
        width = 4865
        height = 4091
        tile_width = 1000
        tile_height = 1000
        tiles_processed = []
        rect_x = 1000
        rect_y = 3000
        tile_to_process = [rect_x, rect_y]
        was_processed = sice2_v21_io.Sice2V21Io.check_if_tile_was_processed(width, height, tile_width, tile_height,
                                                            tile_to_process, tiles_processed)
        self.assertFalse(was_processed)

        tiles_processed.append(tile_to_process)
        was_processed = sice2_v21_io.Sice2V21Io.check_if_tile_was_processed(width, height, tile_width, tile_height,
                                                                            tile_to_process, tiles_processed)
        self.assertTrue(was_processed)

    def test_rad_to_refl(self):
        rad = 39.03614
        sza = 82.5072
        flux = 1472.06
        refl = sice2_v21_utils.Sice2V21Utils.rad_to_refl(rad, sza, flux)
        self.assertAlmostEqual(0.6388, refl, 3)


# suite = unittest.TestLoader().loadTestsFromTestCase(TestSice2)
# unittest.TextTestRunner(verbosity=2).run(suite)
