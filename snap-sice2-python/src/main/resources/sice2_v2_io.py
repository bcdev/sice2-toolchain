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
import xarray as xr
# import netCDF4
import numpy as np
import rioxarray
import rasterio as rio


class SiceV2Io(object):
    def __init__(self, dirname):
        """
        SiceV2Io initialize

        :param dirname: input directory name (tif, zarr or SEN3 folder)
        """
        self.dirname = dirname

        if dirname.endswith(".zarr"):
            self._get_size_zarr()
            self.open = self.open_zarr
        elif os.path.exists(os.path.join(dirname, "Oa01_radiance.nc")):
            self._get_size_satpy()
            self.open = self.open_satpy

        elif os.path.exists(os.path.join(dirname, "r_TOA_01.tif")):
            self._get_size_tif()

            self.open = self.open_tif

    def _open_tif(self, filename):
        """
        Opens tif file

        :param filename:
        :return: rio object
        """
        return rio.open(os.path.join(self.dirname, filename))

    def _get_size_tif(self):
        """
        Extracts size of tif file

        :return:
        """
        self.meta = self._open_tif("r_TOA_01.tif").meta
        self.original_width = self.meta["width"]
        self.original_height = self.meta["height"]

    def open_tif(self, x0=0, y0=0, width=None, height=None):
        """
        Opens tif file.

        :param x0:
        :param y0:
        :param width:
        :param height:
        :return:
        """

        # if os.path.isfile(self.dirname):
        #     raise NotImplementedError("this part needs to be cleaned")
        #     InputFolder = os.path.dirname(os.path.dirname(dirname_or_filename)) + '/'
        #     print('Text file input')
        #     # data_in = pd.read_csv(sys.argv[1])
        #     data_in = pd.read_csv(sys.argv[1])
        #     self.toa = np.expand_dims(data_in[[c for c in data_in.columns if c.find('reflec') >= 0]].to_numpy().transpose(), axis=2)

        #     self.ozone = np.expand_dims(data_in['total_ozone'], axis=1)
        #     self.water = np.expand_dims(data_in['total_columnar_water_vapour'], axis=1)
        #     self.sza = np.expand_dims(data_in['sza'], axis=1)
        #     self.saa = np.expand_dims(data_in['saa'], axis=1)
        #     self.vza = np.expand_dims(data_in['vza'], axis=1)
        #     self.vaa = np.expand_dims(data_in['vaa'], axis=1)
        #     self.height = np.expand_dims(data_in['altitude'], axis=1)

        #     self.sza[np.isnan(toa[0, :, :])] = np.nan
        #     self.saa[np.isnan(toa[0, :, :])] = np.nan
        #     self.vza[np.isnan(toa[0, :, :])] = np.nan
        #     self.vaa[np.isnan(toa[0, :, :])] = np.nan

        # # %% ========= input tif ===============
        if not os.path.isdir(self.dirname):
            raise Exception("dirname must be a directory")

        def read_tif(filename):
            chunks = None
            chunks = "auto"
            # data = xr.open_rasterio(os.path.join(self.dirname, filename), chunks=chunks).squeeze(dim='band', drop=True)
            data = rioxarray.open_rasterio(
                os.path.join(self.dirname, filename), chunks=chunks
            ).squeeze(dim="band", drop=True)

            if width is not None:
                data = data.isel(x=slice(x0, x0 + width))
            if height is not None:
                data = data.isel(y=slice(y0, y0 + height))

            return data.stack(xy=("x", "y")).compute()

        self.meta["transform"] = rio.transform.Affine(
            1.0, 0.0, 0.0, 0.0, -1.0, 0.0
        )  # to improve. This is invalid in some cases
        self.meta.update(compress="DEFLATE")
        self.toa = []
        for i in range(21):
            try:
                dat = read_tif(f"r_TOA_{i + 1:02}.tif")
            except:
                if i in [16, 20]:
                    raise Exception("Missing the necessary bands")
                else:
                    print(
                        "Cannot load ",
                        "r_TOA_" + str(i + 1).zfill(2) + ".tif, replacing by nans",
                        )
                    dat = xr.full_like(self.toa[0], fill_value=np.nan)
            self.toa.append(dat)
        self.olci_scene = xr.Dataset()
        self.olci_scene["toa"] = xr.concat(self.toa, dim="band")
        # print("toa=", self.toa.coords)

        self.olci_scene["ozone"] = read_tif("O3.tif")
        # self.water = read_tif('WV.tif')  # save memory, it is not used
        self.olci_scene["sza"] = read_tif("SZA.tif")
        self.olci_scene["saa"] = read_tif("SAA.tif")
        self.olci_scene["vza"] = read_tif("OZA.tif")
        self.olci_scene["vaa"] = read_tif("OAA.tif")
        self.olci_scene["elevation"] = read_tif("height.tif").astype(np.float64)

        mask = ~np.isnan(self.olci_scene.toa.sel(band=0))
        self.olci_scene = self.olci_scene.where(mask)

        t = self.olci_scene.elevation.unstack("xy")
        self.meta["width"] = len(t.x)
        self.meta["height"] = len(t.y)

    def _get_size_satpy(self):
        """
        Extracts size of OLCI L1b SEN3 file

        :return:
        """
        filename = os.path.join(self.dirname, "Oa01_radiance.nc")
        rootgrp = netCDF4.Dataset(filename, "r")
        self.original_width = rootgrp.dimensions["columns"].size
        self.original_height = rootgrp.dimensions["rows"].size

    def open_satpy(self, x0=0, y0=0, width=None, height=None):
        """
        Opens OLCI L1b SEN3 file.

        :param x0:
        :param y0:
        :param width:
        :param height:
        :return:
        """
        import satpy  # this is not good practice but avoid satpy to be a compulsary dependence

        filenames = glob(os.path.join(self.dirname, "*.nc"))

        scene = satpy.Scene(reader="olci_l1b", filenames=filenames)

        variables = {
            "solar_azimuth_angle": "saa",
            "solar_zenith_angle": "sza",
            "satellite_azimuth_angle": "vaa",
            "satellite_zenith_angle": "vza",
            "total_ozone": "ozone",
            "altitude": "elevation"
            # 'longitude': 'longitude',
            # 'latitude': 'latitude'
        }
        scene.load(list(variables.keys()))

        islice = {}
        if width is not None:
            islice["x"] = slice(x0, x0 + width)
        if height is not None:
            islice["y"] = slice(y0, y0 + height)

        def get_var(variable):
            # return the variable and remove what needs to be remove
            data = scene[variable].isel(islice).compute().stack(xy=("x", "y"))
            data.attrs = (
                {}
            )  # remove attributes, due to some conflict with tà_zarr being unable to serialize datatime
            if "crs" in data.coords:
                del data.coords["crs"]  # idem. zarr complains
            return data

        for variable in variables:
            setattr(self, variables[variable], get_var(variable))
        scene.unload()  # maybe useless
        coef = 1 / np.cos(np.deg2rad(self.sza)) / 100.0

        bands = [f"Oa{i:02}" for i in range(1, 22)]
        scene.load(bands)

        scene.load(
            [satpy.DataQuery(name=band, calibration="reflectance") for band in bands]
        )
        self.toa = []
        for band in bands:
            self.toa.append(np.clip(get_var(band) * coef, 0, 1))
        self.toa = xr.concat(self.toa, dim="band")
        if "crs" in self.toa.coords:
            del self.toa.coords["crs"]  # idem. zarr complains

        scene.unload()  # probably useless

    def _get_size_zarr(self):
        """
        Extracts size of zarr file

        :return:
        """
        ds = xr.open_zarr(self.dirname)
        self.original_width = len(ds.x)
        self.original_height = len(ds.y)

    def open_zarr(self, x0=0, y0=0, width=None, height=None):
        """
        Opens zarr file.

        :param x0:
        :param y0:
        :param width:
        :param height:
        :return:
        """

        variables = {
            "solar_azimuth_angle": "saa",
            "solar_zenith_angle": "sza",
            "satellite_azimuth_angle": "vaa",
            "satellite_zenith_angle": "vza",
            "total_ozone": "ozone",
            "altitude": "elevation"
            # 'longitude': 'longitude',
            # 'latitude': 'latitude'
        }

        ds = xr.open_zarr(self.dirname)

        islice = {}
        if width is not None:
            islice["x"] = slice(x0, x0 + width)
        if height is not None:
            islice["y"] = slice(y0, y0 + height)

        def get_var(variable):
            # return the variable and remove what needs to be remove
            return ds[variable].isel(islice).stack(xy=("x", "y")).compute()

        for variable in variables:
            setattr(self, variables[variable], get_var(variable))

        bands = [f"Oa{i:02}" for i in range(1, 22)]
        self.toa = []
        for band in bands:
            self.toa.append(get_var(band))
        self.toa = xr.concat(self.toa, dim="band")


    def write_output(self, snow, output_folder):
        """
        Write snow retrievals to tif files (one file per variable)

        :param snow: xarray Dataset with snow retrievals
        :param output_folder: output folder with tif files
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

        def da_to_tif(da, file_path):
            da = da.unstack(dim="xy").transpose("y", "x")
            da.reindex(y=list(reversed(da.y))).rio.to_raster(file_path,
                                                             dtype='float32',compress='DEFLATE')

        print('Printing out:')
        for var in ["diameter", "area", "rp3", "rs3", "isnow", "r0", "al", 'factor']:
            print(var)
            da_to_tif(snow[var],
                      os.path.join(output_folder, file_name_list[var] + ".tif"))

        da_to_tif(snow.alb_sph.sel(band=0), output_folder + '/alb_sph_01_solved.tif')
        da_to_tif(snow.alb_sph_direct.sel(band=0), output_folder + '/alb_sph_01.tif')
        da_to_tif(snow.rp.sel(band=0), output_folder + '/alb_pl_01_solved.tif')
        da_to_tif(snow.rp_direct.sel(band=0), output_folder + '/alb_pl_01.tif')
        # for var in ['BXXX', ]:
        #     var = OLCI_scene[var].unstack(dim='xy').transpose('y', 'x').rio.to_raster(os.path.join(OutputFolder, file_name_list[var] + '.tif'))
