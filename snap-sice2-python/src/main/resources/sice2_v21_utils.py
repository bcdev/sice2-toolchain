# -*- coding: utf-8 -*-
"""
@author: olafd

"""

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
