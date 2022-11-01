# -*- coding: utf-8 -*-
"""
@author: olafd

"""
import os
from math import ceil
import numpy as np

import sice2_constants


class Sice2V21Utils(object):
    def __init__(self):
        pass

    @staticmethod
    def check_if_tile_was_processed(width, height, tile_width, tile_height, tile_to_process, tiles_processed):
        num_tiles_x = ceil(width/tile_width)
        num_tiles_y = ceil(height/tile_height)
        tiles_x = []
        for i in range(num_tiles_x):
            tiles_x.append((i+1)*tile_width)
        tiles_y = []
        for i in range(num_tiles_y):
            tiles_y.append((i+1)*tile_height)

        return tile_to_process in tiles_processed

    @staticmethod
    def rad_to_refl(rad, sza, flux):
        # (float) ((rad * Math.PI) / (e0 * Math.cos(sza * RAD_PER_DEG)));
        return (rad * np.pi) / (flux * np.cos(np.deg2rad(sza)))

    @staticmethod
    def get_tif_source_product_paths(tif_input_directory):
        """

        :param tif_input_directory:
        :return:
        """

        # mandatory tif files are sice2_constants.mandatory_tif_inputs_for_retrieval:
        # (more r_TOA_xx optional)
        #

        # check if all inputs exist:
        for src in sice2_constants.mandatory_tif_inputs_for_retrieval:
            if not os.path.exists(tif_input_directory + os.sep + src):
                raise Exception("Missing necessary band " + src + " - cannot proceed.")

        tif_source_product_list = []
        for src in os.listdir(tif_input_directory):
            if src.lower().endswith('.tif'):
                if src in sice2_constants.mandatory_tif_inputs_for_retrieval or src.startswith('r_TOA'):
                    tif_source_product_list.append(tif_input_directory + os.sep + src)

        return tif_source_product_list
