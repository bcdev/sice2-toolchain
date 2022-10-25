# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""

import os
from glob import glob

import numpy as np
import xarray as xr


class Sice2V21TifIo(object):

    def __init__(self, dirname):
        """

        :param dirname:
        """

        self.dirname = dirname
        # self.open = self.open_tif


    def open_tif(self, x0=0, y0=0, width=None, height=None):
        """

        :param x0:
        :param y0:
        :param width:
        :param height:
        :return:
        """

        # # %% ========= input tif ===============
        if not os.path.isdir(self.dirname):
            raise Exception("dirname must be a directory")

        print('hier 1')
        # import rasterio as rio
        print('hier 2')
        # from rasterio import open as rioopen
        # from rasterio.transform import Affine as riotransformaffime

        # self.meta = rio.open(os.path.join(self.dirname, 'r_TOA_01.tif')).meta
        # self.meta = rioopen(os.path.join(self.dirname, 'r_TOA_01.tif')).meta

        # self.original_width = self.meta['width']
        # self.original_height = self.meta['height']

        def read_tif(filename):
            chunks = 'auto'
            data = xr.open_rasterio(os.path.join(self.dirname, filename), chunks=chunks).squeeze(dim='band', drop=True)
            # data = rioxarray.open_rasterio(os.path.join(self.dirname, filename), chunks=chunks).squeeze(dim='band', drop=True)

            if width is not None:
                data = data.isel(x=slice(x0, x0 + width))
            if height is not None:
                data = data.isel(y=slice(y0, y0 + height))

            return data.stack(xy=("x", "y")).compute()

        # self.meta['transform'] = rio.transform.Affine(1.0, 0.0, 0.0, 0.0, -1.0, 0.0)  # to improve. This is invalid in some cases
        # self.meta['transform'] = riotransformaffime(1.0, 0.0, 0.0, 0.0, -1.0, 0.0)  # to improve. This is invalid in some cases
        # self.meta.update(compress='DEFLATE')
        self.toa = []
        for i in range(21):
            try:
                dat = read_tif(f'r_TOA_{i + 1:02}.tif')
            except:
                if i in [16, 20]:
                    raise Exception('Missing the necessary bands')
                else:
                    print('Cannot load ', 'r_TOA_'+str(i+1).zfill(2)+'.tif, replacing by nans')
                    dat = xr.full_like(self.toa[0], fill_value=np.nan)
            self.toa.append(dat)
        self.olci_scene = xr.Dataset()
        self.olci_scene['toa'] = xr.concat(self.toa, dim='band')
        # print("toa=", self.toa.coords)

        self.olci_scene['ozone'] = read_tif('O3.tif')
        # self.water = read_tif('WV.tif')  # save memory, it is not used
        self.olci_scene['sza'] = read_tif('SZA.tif')
        self.olci_scene['saa'] = read_tif('SAA.tif')
        self.olci_scene['vza'] = read_tif('OZA.tif')
        self.olci_scene['vaa'] = read_tif('OAA.tif')
        self.olci_scene['elevation'] = read_tif('height.tif').astype(np.float64)

        mask = ~np.isnan(self.olci_scene.toa.sel(band=0))
        self.olci_scene = self.olci_scene.where(mask)

        t = self.olci_scene.elevation.unstack('xy')
        print('hier 3')
        # self.meta['width'] = len(t.x)
        # self.meta['height'] = len(t.y)
        self.original_width = len(t.x)
        self.original_height = len(t.x)
        print('hier 4')


    def write_output(self, snow, OutputFolder):
        """

        :param snow:
        :param OutputFolder:
        :return:
        """
        file_name_list = {
            "BXXX": "O3_SICE",
            "diameter": "grain_diameter",
            "area": "snow_specific_area",
            "al": "al",
            "r0": "r0",
            "isnow": "isnow",
            "conc": "conc",
            "rp3": "albedo_bb_planar_sw",
            "rs3": "albedo_bb_spherical_sw",
            "factor": "factor",
        }
        print('Printing out:')
        for var in ["diameter", "area", "rp3", "rs3", "isnow", "r0", "al"]:
            print(var)
            snow[var].unstack(dim="xy").transpose("y", "x").rio.to_raster(
                os.path.join(OutputFolder, file_name_list[var] + ".tif")
            )
        snow.alb_sph.sel(band=0).unstack(dim="xy").transpose("y", "x").rio.to_raster(OutputFolder+'/alb_sph_01_solved.tif')
        snow.rp.sel(band=0).unstack(dim="xy").transpose("y", "x").rio.to_raster(OutputFolder+'/alb_pl_01_solved.tif')


# class Sice2V21SafeIo(object):
#     def __init__(self, dirname):
#         """
#
#         :param dirname:
#         """
#
#         self.dirname = dirname
#         self.myopen = self.open_satpy
#
#     def test_open_satpy(self):
#         import satpy
#
#     def open_satpy(self, x0=0, y0=0, width=None, height=None, with_geom=True):
#         """
#
#         :param x0:
#         :param y0:
#         :param width:
#         :param height:
#         :param with_geom:
#         :return:
#         """
#
#         print('hier 1')
#         import satpy  # this is not good practice but avoid satpy to be a compulsary dependence
#         print('hier 2')
#
#         filenames = glob(os.path.join(self.dirname, "*.nc"))
#
#         scene = satpy.Scene(reader="olci_l1b", filenames=filenames)
#         print('hier 3')
#
#         variables = {
#             'solar_azimuth_angle': 'saa',
#             'solar_zenith_angle': 'sza',
#             'satellite_azimuth_angle': 'vaa',
#             'satellite_zenith_angle': 'vza',
#             'total_ozone': 'ozone',
#             'altitude': 'elevation'
#         }
#
#         scene.load(list(variables.keys()))
#         print('hier 4')
#
#         width = scene.to_xarray_dataset().dims['x']
#         height = scene.to_xarray_dataset().dims['y']
#
#         islice = {}
#         if width is not None:
#             islice['x'] = slice(x0, x0 + width)
#         if height is not None:
#             islice['y'] = slice(y0, y0 + height)
#
#         print('hier 5')
#         def get_var(variable):
#             # return the variable and remove what needs to be remove
#             data = scene[variable].isel(islice).compute().stack(xy=("x", "y"))
#             data.attrs = {}  # remove attributes, due to some conflict with tà_zarr being unable to serialize datatime
#             if 'crs' in data.coords:
#                 del data.coords['crs']  # idem. zarr complains
#             return data
#
#         self.olci_scene = xr.Dataset()
#         print('hier 6')
#         if with_geom:
#             scene.load(['longitude', 'latitude'])
#             self.olci_scene = self.olci_scene.assign_coords(longitude=get_var('longitude'),
#                                                             latitude=get_var('latitude'))
#         for variable in variables:
#             self.olci_scene[variables[variable]] = get_var(variable)
#
#         scene.unload()  # maybe useless
#
#         coef = 1 / np.cos(np.deg2rad(self.olci_scene['sza'])) / 100.
#
#         bands = [f'Oa{i:02}' for i in range(1, 22)]
#         scene.load(bands)
#
#         scene.load([satpy.DataQuery(name=band, calibration='reflectance') for band in bands])
#         toa = []
#         for band in bands:
#             toa.append(np.clip(get_var(band) * coef, 0, 1))
#         self.olci_scene['toa'] = xr.concat(toa, dim='band')
#
#         if 'crs' in self.olci_scene['toa'].coords:
#             del self.olci_scene['toa'].coords['crs']  # idem. zarr complains
#
#         scene.unload()  # probably useless



