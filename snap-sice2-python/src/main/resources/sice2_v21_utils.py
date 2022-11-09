# -*- coding: utf-8 -*-
"""
@author: olafd

"""
import os
import re
from math import ceil
import numpy as np

import sice2_constants


def check_if_tile_was_processed(width, height, tile_width, tile_height, tile_to_process, tiles_processed):
    num_tiles_x = ceil(width / tile_width)
    num_tiles_y = ceil(height / tile_height)
    tiles_x = []
    for i in range(num_tiles_x):
        tiles_x.append((i + 1) * tile_width)
    tiles_y = []
    for i in range(num_tiles_y):
        tiles_y.append((i + 1) * tile_height)

    return tile_to_process in tiles_processed


def rad_to_refl(rad, sza, flux):
    # (float) ((rad * Math.PI) / (e0 * Math.cos(sza * RAD_PER_DEG)));
    return (rad * np.pi) / (flux * np.cos(np.deg2rad(sza)))


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


def get_condition_from_valid_pixel_expr(flagname, expr, flag_coding_dict):
    """
    Taken from OMAPS.

    :param flagname:
    :param expr:
    :param flag_coding_dict:
    :return:
    """
    # consider Python logical and comparison operators (https://www.w3schools.com/python/python_operators.asp)
    # e.g. valid_pixel_expr = '(WQSF.WATER  and Oa10_reflectance < 0.13) or (WQSF.LAND and not WQSF.CLOUD)'
    # extract strings (e.g. WQSF.WATER, WQSF.LAND, and WQSF.CLOUD) from expr
    flagname_starts = re.finditer(flagname + '.', expr)
    flagname_start_indices = [match.start() for match in flagname_starts]
    flag_substrings = []
    for index in flagname_start_indices:
        # search for blank
        flag_substring_end_index = expr[index:].find(' ')
        if flag_substring_end_index < 0:
            # search for close bracket
            flag_substring_end_index = expr[index:].find(')')
            if flag_substring_end_index < 0:
                # must be the end of expression
                flag_substring_end_index = len(expr)
        flag_substring = expr[index:index + flag_substring_end_index]
        flag_substrings.append(flag_substring)

    # look up flag coding values (bit indices) from flag_coding_dict:
    flag_substring_bit_indices = []
    for flag_substring in flag_substrings:
        flag_substring_bit_index = flag_coding_dict[flag_substring]
        flag_substring_bit_indices.append(int(np.log2(flag_substring_bit_index)))

    # transform the flag terms into a bitwise exploration
    # e.g. valid_pixel_expr =
    #        '(WQSF.WATER  and Oa10_reflectance < 0.13) or (WQSF.LAND and not WQSF.CLOUD)'
    # -->    'bool((wqsf & (1 << 1) and Oa10_reflectance < 0.13) or (wqsf & (1 << 2) and not wqsf & (1 << 3)))'
    for i in range(len(flag_substrings)):
        replacement = flagname + ' & (1 << ' + str(flag_substring_bit_indices[i]) + ')'
        expr_new = expr.replace(flag_substrings[i], replacement, 1)
        expr = expr_new

    # return string with boolean condition which can e.g. be further evaluated like: if eval(condition): ...
    return 'bool(' + expr + ')'


def get_valid_expression_filter_array(condition, variables_in_expr_dict, width, height):
    """

    :param condition:
    :param variables_in_expr_dict:
    :param width:
    :param height:
    :return:
    """
    # method extracted for test purpose
    valid_expr_arr = np.full((width, height), True, dtype=bool)

    # e.g. valid_pixel_expr = '(WQSF.WATER  and Oa10_reflectance < 0.13) or (WQSF.LAND and not WQSF.CLOUD)'
    # bool((wqsf & (1 << 1) and Oa10_reflectance < 0.13) or (wqsf & (1 << 2) and not wqsf & (1 << 3)))

    # or e.g. valid_pixel_expr = '(bitmask.LAND and Rrs_412 < 0.0013) or (bitmask.CLOUD_BASE and not bitmask.CASE2)'
    # bool((bitmask & (1 << 0) and Rrs_412 < 0.0013) or (bitmask & (1 << 1) and not bitmask & (1 << 10)))

    for i in range(width):
        for j in range(height):
            for name, var_arr in variables_in_expr_dict.items():
                command_str = name + ' = ' + str(var_arr[i][j])
                exec(command_str)
            valid_expr_arr[i][j] = eval(condition)

    print('leave get_valid_expression_filter_array... ')
    return valid_expr_arr
