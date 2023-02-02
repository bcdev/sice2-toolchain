import sys
import os

sys.path.append(os.path.expanduser('~') + os.sep + '.snap' + os.sep + 'snap-python')
import esa_snappy
# If a Java type is needed which is not imported by snappy by default it can be retrieved manually.
# First import jpy
from esa_snappy import jpy

HashMap = jpy.get_type('java.util.HashMap')

from esa_snappy import ProductIO
from esa_snappy import GPF

import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) != 2:
    print('Usage: test_scda.py </path/to/slstr_l1_product')
    sys.exit(1)

p = ProductIO.readProduct(sys.argv[1])
rad550_band = p.getBand('S1_radiance_an')
bt37_band = p.getBand('S7_BT_in')
bt11000_band = p.getBand('S8_BT_in')
bt12000_band = p.getBand('S9_BT_in')
wRad = rad550_band.getRasterWidth()
hRad = rad550_band.getRasterHeight()
rad13_data = np.zeros(wRad * hRad, np.float32)
rad550_band.readPixels(0, 0, wRad, hRad, rad13_data)

wBT = bt37_band.getRasterWidth()
hBT = bt37_band.getRasterHeight()

# band subset for BT bands:
parameters = HashMap()
parameters.put('bandNames', 'S7_BT_in,S8_BT_in,S9_BT_in')
operator_name = 'Subset'
bt_subset_product = GPF.createProduct(operator_name, parameters, p)
for subset_band in bt_subset_product.getBands():
    print('subset band: ' + subset_band.getName())
bt37_subset_band = bt_subset_product.getBand('S7_BT_in')
bt11000_subset_band = bt_subset_product.getBand('S8_BT_in')
bt12000_subset_band = bt_subset_product.getBand('S9_BT_in')

p.dispose()

# resample BT bands:
parameters = HashMap()
parameters.put('targetWidth', wBT * 2)
parameters.put('targetHeight', hBT * 2)
operator_name = 'Resample'
bt_subset_product_resampled = GPF.createProduct(operator_name, parameters, bt_subset_product)
for subset_band_resampled in bt_subset_product_resampled.getBands():
    print('resampled subset band: ' + subset_band_resampled.getName())
bt37_subset_band_resampled = bt_subset_product_resampled.getBand('S7_BT_in')
bt11000_subset_band_resampled = bt_subset_product_resampled.getBand('S8_BT_in')
bt12000_subset_band_resampled = bt_subset_product_resampled.getBand('S9_BT_in')

wBTResampled = bt37_subset_band_resampled.getRasterWidth()
hBTResampled = bt37_subset_band_resampled.getRasterHeight()

rad13_data.shape = hRad, wRad
imgplot = plt.imshow(rad13_data)
imgplot.write_png('S1_radiance_an.png')
