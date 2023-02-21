# -*- coding: utf-8 -*-
"""
This script contains I/O utility functions needed by the SICE2 SNAP processors.

@author: Olaf Danne, BC (Brockmann Consult)
"""

from esa_snappy import FlagCoding

from esa_snappy import jpy

BitSetter = jpy.get_type('org.esa.snap.core.util.BitSetter')
IndexCoding = jpy.get_type('org.esa.snap.core.datamodel.IndexCoding')
ColorPaletteDef = jpy.get_type('org.esa.snap.core.datamodel.ColorPaletteDef')
Point = jpy.get_type('org.esa.snap.core.datamodel.ColorPaletteDef$Point')
ImageInfo = jpy.get_type('org.esa.snap.core.datamodel.ImageInfo')
BandMathsType = jpy.get_type('org.esa.snap.core.datamodel.Mask$BandMathsType')
Color = jpy.get_type('java.awt.Color')


def create_pol_type_index_coding(flag_id):
    """
    Creates index coding for snow pollution type flag.

    :param flag_id: the flag id
    :return: pol_type_index_coding, image_info
    """
    pol_type_index_coding = IndexCoding(flag_id)
    points = []
    pol_type_flags = ['SOOT', 'DUST', 'OTHER', 'MIXTURE']
    pol_type_descr = ['Soot', 'Dust', 'Other', 'Mixture']
    pol_type_colors = [Color.RED, Color.PINK, Color.YELLOW, Color.ORANGE]
    for i in range(len(pol_type_flags)):
        pol_type_index_coding.addIndex(pol_type_flags[i], i, pol_type_descr[i])
        points.append(Point(i, pol_type_colors[i], pol_type_descr[i]))

    cpd = ColorPaletteDef(points)
    image_info = ImageInfo(cpd)

    return pol_type_index_coding, image_info


def create_snow_retrieval_flag_coding(flag_id):
    """
    Creates flag coding for snow retrieval flag.

    :param flag_id: the flag id
    :return: snow_retrieval_flag_coding
    """
    snow_retrieval_flag_coding = FlagCoding(flag_id)
    snow_retrieval_flag_coding.addFlag("CLEAN_SNOW", BitSetter.setFlag(0, 0), "Clean snow")
    snow_retrieval_flag_coding.addFlag("POLLUTED_SNOW", BitSetter.setFlag(0, 1), "Polluted snow")
    snow_retrieval_flag_coding.addFlag("PARTIALLY_SNOW_COVERED", BitSetter.setFlag(0, 2),
                                       "Partially snow covered pixel")
    snow_retrieval_flag_coding.addFlag("SZA_OOR", BitSetter.setFlag(0, 3),
                                       "SZA out of range (< 75 deg), no retrival")
    snow_retrieval_flag_coding.addFlag("RTOA_01_OOR", BitSetter.setFlag(0, 4),
                                       "TOA reflectance at band 1 < 0.1, no retrieval")
    snow_retrieval_flag_coding.addFlag("RTOA_21_OOR", BitSetter.setFlag(0, 5),
                                       "TOA reflectance at band 21 < 0.2, no retrieval")
    snow_retrieval_flag_coding.addFlag("GRAIN_DIAMETER_OOR", BitSetter.setFlag(0, 6),
                                       "grain_diameter < 0.1, no retrieval (potential cloud flag)")
    snow_retrieval_flag_coding.addFlag("SPH_ALB_NEG", BitSetter.setFlag(0, 7),
                                       "Retrieved spherical albedo negative in band 1, 2 or 3")
    snow_retrieval_flag_coding.addFlag("SPH_ALB_NO_SOLUTION", BitSetter.setFlag(0, 8),
                                       "Impossible to solve snow spherical albedo equation")

    return snow_retrieval_flag_coding


def create_snow_retrieval_bitmask(snow_product):
    """
    Creates a bit mask for the snow retrieval flag.

    :param snow_product: the snow retrieval product
    :return: void
    """
    index = 0
    w = snow_product.getSceneRasterWidth()
    h = snow_product.getSceneRasterHeight()

    mask = BandMathsType.create("CLEAN_SNOW", "Clean snow", w, h,
                                "isnow.CLEAN_SNOW", Color.CYAN, 0.5)
    snow_product.getMaskGroup().add(index, mask)
    index = index + 1

    mask = BandMathsType.create("POLLUTED_SNOW", "Polluted snow", w, h,
                                "isnow.POLLUTED_SNOW", Color.PINK, 0.5)
    snow_product.getMaskGroup().add(index, mask)
    index = index + 1

    mask = BandMathsType.create("PARTIALLY_SNOW_COVERED", "Partially snow covered pixel", w, h,
                                "isnow.PARTIALLY_SNOW_COVERED", Color.ORANGE, 0.5)
    snow_product.getMaskGroup().add(index, mask)
    index = index + 1

    mask = BandMathsType.create("SZA_OOR", "SZA out of range (< 75 deg), no retrival", w, h,
                                "isnow.SZA_OOR", Color.BLUE, 0.5)
    snow_product.getMaskGroup().add(index, mask)
    index = index + 1

    mask = BandMathsType.create("RTOA_01_OOR", "TOA reflectance at band 21 < 0.1, no retrieval", w, h,
                                "isnow.RTOA_01_OOR", Color.GREEN, 0.5)
    snow_product.getMaskGroup().add(index, mask)
    index = index + 1

    mask = BandMathsType.create("RTOA_21_OOR", "TOA reflectance at band 1 < 0.2, no retrieval", w, h,
                                "isnow.RTOA_21_OOR", Color(50, 150, 150), 0.5)
    snow_product.getMaskGroup().add(index, mask)
    index = index + 1

    mask = BandMathsType.create("GRAIN_DIAMETER_OOR", "grain_diameter < 0.1, no retrieval (potential cloud flag)", w, h,
                                "isnow.GRAIN_DIAMETER_OOR", Color.YELLOW, 0.5)
    snow_product.getMaskGroup().add(index, mask)
    index = index + 1

    mask = BandMathsType.create("SPH_ALB_NEG", "Retrieved spherical albedo negative in band 1, 2 or 3", w, h,
                                "isnow.SPH_ALB_NEG", Color(150, 50, 50), 0.5)
    snow_product.getMaskGroup().add(index, mask)
    index = index + 1

    mask = BandMathsType.create("SPH_ALB_NO_SOLUTION", "Impossible to solve snow spherical albedo equation", w, h,
                                "isnow.SPH_ALB_NO_SOLUTION", Color.MAGENTA, 0.5)
    snow_product.getMaskGroup().add(index, mask)


def create_scda_flag_coding(flag_id):
    """
    Creates flag coding for SCDA cloud flag.

    :param flag_id: the flag id
    :return: scda_flag_coding
    """
    scda_flag_coding = FlagCoding(flag_id)
    scda_flag_coding.addFlag("SCDA_INVALID", BitSetter.setFlag(0, 0), "Invalid")
    scda_flag_coding.addFlag("SCDA_CLEAR", BitSetter.setFlag(0, 1), "Clear")
    scda_flag_coding.addFlag("SCDA_CLOUDY", BitSetter.setFlag(0, 2), "Cloudy")
    return scda_flag_coding


def create_scda_bitmask(scda_product):
    """
    Creates a bit mask for the SCDA cloud flag.

    :param scda_product: the SCDA product
    :return: void
    """
    index = 0
    w = scda_product.getSceneRasterWidth()
    h = scda_product.getSceneRasterHeight()

    mask = BandMathsType.create("SCDA_INVALID", "Invalid", w, h, "scda_cloud_mask.SCDA_INVALID", Color.red, 0.5)
    scda_product.getMaskGroup().add(index, mask)
    index = index + 1

    mask = BandMathsType.create("SCDA_CLEAR", "Clear", w, h, "scda_cloud_mask.SCDA_CLEAR", Color.CYAN, 0.5)
    scda_product.getMaskGroup().add(index, mask)
    index = index + 1

    mask = BandMathsType.create("SCDA_CLOUDY", "Cloudy", w, h, "scda_cloud_mask.SCDA_CLOUDY", Color.yellow, 0.5)
    scda_product.getMaskGroup().add(index, mask)


def create_idepix_bitmask(snow_product):
    """
    Creates a bit mask for the Idepix classification flag.

    :param snow_product: the snow retrieval product
    :return: void
    """
    index = 0
    w = snow_product.getSceneRasterWidth()
    h = snow_product.getSceneRasterHeight()

    mask = BandMathsType.create("IDEPIX_INVALID", "Invalid pixels", w, h,
                                "pixel_classif_flags.IDEPIX_INVALID", Color.MAGENTA, 0.5)
    snow_product.getMaskGroup().add(index, mask)
    index = index + 1

    mask = BandMathsType.create("IDEPIX_CLOUD", "Pixels which are either cloud_sure or cloud_ambiguous", w, h,
                                "pixel_classif_flags.IDEPIX_CLOUD", Color.RED, 0.5)
    snow_product.getMaskGroup().add(index, mask)
    index = index + 1

    mask = BandMathsType.create("IDEPIX_CLOUD_BUFFER",
                                "A buffer of n pixels around a cloud. n is a user supplied parameter. Applied to "
                                "pixels masked as 'cloud'", w, h,
                                "pixel_classif_flags.IDEPIX_CLOUD_BUFFER", Color.YELLOW, 0.5)
    snow_product.getMaskGroup().add(index, mask)
    index = index + 1

    mask = BandMathsType.create("IDEPIX_CLOUD_SHADOW", "Pixel is affected by a cloud shadow", w, h,
                                "pixel_classif_flags.IDEPIX_CLOUD_SHADOW", Color.ORANGE, 0.5)
    snow_product.getMaskGroup().add(index, mask)
    index = index + 1

    mask = BandMathsType.create("IDEPIX_SNOW_ICE", "Clear snow/ice pixels", w, h,
                                "pixel_classif_flags.IDEPIX_SNOW_ICE", Color.CYAN, 0.5)
    snow_product.getMaskGroup().add(index, mask)
    index = index + 1

    mask = BandMathsType.create("IDEPIX_COASTLINE", "Pixels at a coastline", w, h,
                                "pixel_classif_flags.IDEPIX_COASTLINE", Color.PINK, 0.5)
    snow_product.getMaskGroup().add(index, mask)
    index = index + 1

    mask = BandMathsType.create("IDEPIX_LAND", "Land pixels", w, h,
                                "pixel_classif_flags.IDEPIX_LAND", Color.GREEN, 0.5)
    snow_product.getMaskGroup().add(index, mask)
