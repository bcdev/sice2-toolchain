import sys
import os

# for SNAPPY usage see
# https://senbox.atlassian.net/wiki/spaces/SNAP/pages/19300362/How+to+use+the+SNAP+API+from+Python

# Append esa_snappy installation dir to path:
sys.path.append(os.path.expanduser('~') + os.sep + '.snap' + os.sep + 'snap-python')

# import esa_snappy
# If a Java type is needed which is not imported by snappy by default it can be retrieved manually.
# First import jpy
from esa_snappy import jpy

HashMap = jpy.get_type('java.util.HashMap')

from esa_snappy import Product
from esa_snappy import ProductIO
from esa_snappy import ProductUtils
from esa_snappy import GPF

import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) != 3:
    print('Usage: test_scda.py </path/to/olci_l1_subset_product> </path/to/slstr_l1_product>')
    sys.exit(1)

print('sys.argv[1]: ' + sys.argv[1])
print('sys.argv[2]: ' + sys.argv[2])

olci_rad_product = ProductIO.readProduct(sys.argv[1])
# band subset for one OLCI band:
#parameters = HashMap()
#parameters.put('bandNames', 'Oa10_radiance')
#operator_name = 'Subset'
#olci_rad_subset_product = GPF.createProduct(operator_name, parameters, p_olci_l1)

rad10_band = olci_rad_product.getBand('Oa10_radiance')
wRad = rad10_band.getRasterWidth()
hRad = rad10_band.getRasterHeight()
olci_rad_subset_product = Product('Rad subset', 'Rad subset', wRad, hRad)
ProductUtils.copyGeoCoding(olci_rad_product, olci_rad_subset_product)
ProductUtils.copyBand('Oa10_radiance', olci_rad_product, olci_rad_subset_product, True)
for subset_band in olci_rad_subset_product.getBands():
    print('rad subset band: ' + subset_band.getName())

olci_rad_product.dispose()

p_slstr_l1 = ProductIO.readProduct(sys.argv[2])
rad550_band = p_slstr_l1.getBand('S1_radiance_an')
bt37_band = p_slstr_l1.getBand('S7_BT_in')
bt11000_band = p_slstr_l1.getBand('S8_BT_in')
bt12000_band = p_slstr_l1.getBand('S9_BT_in')
wRad = rad550_band.getRasterWidth()
hRad = rad550_band.getRasterHeight()
rad13_data = np.zeros(wRad * hRad, np.float32)
rad550_band.readPixels(0, 0, wRad, hRad, rad13_data)


# band subset for BT bands:
# parameters = HashMap()
# parameters.put('bandNames', 'S7_BT_in,S8_BT_in,S9_BT_in')
# parameters.put('tiePointGridNames', '')
# operator_name = 'Subset'
# bt_subset_product = GPF.createProduct(operator_name, parameters, p_slstr_l1)

wBT = bt37_band.getRasterWidth()
hBT = bt37_band.getRasterHeight()
slstr_bt_subset_product = Product('BT subset', 'BT subset', wBT, hBT)
ProductUtils.copyGeoCoding(p_slstr_l1, slstr_bt_subset_product)
ProductUtils.copyBand('S7_BT_in', p_slstr_l1, slstr_bt_subset_product, True)
ProductUtils.copyBand('S8_BT_in', p_slstr_l1, slstr_bt_subset_product, True)
ProductUtils.copyBand('S9_BT_in', p_slstr_l1, slstr_bt_subset_product, True)

p_slstr_l1.dispose()

for subset_band in slstr_bt_subset_product.getBands():
    print('bt subset band: ' + subset_band.getName())
for subset_tpg in slstr_bt_subset_product.getTiePointGrids():
    print('bt subset tpg: ' + subset_tpg.getName())
bt37_subset_band = slstr_bt_subset_product.getBand('S7_BT_in')
bt11000_subset_band = slstr_bt_subset_product.getBand('S8_BT_in')
bt12000_subset_band = slstr_bt_subset_product.getBand('S9_BT_in')

# collocate the two subsets, OLCI as master:
input_products = HashMap()
input_products.put('master', olci_rad_subset_product)
input_products.put('slave', slstr_bt_subset_product)
parameters = HashMap()
parameters.put('slaveComponentPattern', '${ORIGINAL_NAME}')
operator_name = 'Collocate'
collocate_product = GPF.createProduct(operator_name, parameters, input_products)
for collocate_band in collocate_product.getBands():
    print('collocate band: ' + collocate_band.getName())
    # print('collocate band width: ' + str(collocate_band.getRasterWidth()))
    # print('collocate band height: ' + str(collocate_band.getRasterHeight()))

# resample BT bands:
parameters = HashMap()
parameters.put('targetWidth', wBT * 2)
parameters.put('targetHeight', hBT * 2)
operator_name = 'Resample'
bt_subset_product_resampled = GPF.createProduct(operator_name, parameters, slstr_bt_subset_product)
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
