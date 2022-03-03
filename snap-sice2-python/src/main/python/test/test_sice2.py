import os
import sys
import unittest
from sys import platform


# noinspection PyUnresolvedReferences
class TestSice2(unittest.TestCase):
    def setUp(self):
        print('Platform: ', platform)
        # parent_dir = os.path.dirname(os.path.normpath(os.path.dirname(__file__)))

        resource_root = os.path.dirname(__file__)

        sys.path.append(resource_root)

        SICE2_HOME = os.path.dirname(os.path.abspath(__file__))
        sys.path.append(SICE2_HOME)

        import sice2_algo
        self.sice2algo = sice2_algo.Sice2Algo(-0.2, 0.8)

    # @unittest.skip("skipping test...")
    def test_sice2algo(self):
        print('*** Compute NDSI ***')
        upper_data = 0.2
        lower_data = 0.1

        ndsi = self.sice2algo.compute_ndsi(lower_data, upper_data)

        print('NDSI: ' + str(ndsi))


# suite = unittest.TestLoader().loadTestsFromTestCase(TestSice2)
# unittest.TextTestRunner(verbosity=2).run(suite)
