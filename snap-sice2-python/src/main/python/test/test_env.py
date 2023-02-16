import os
import sys

import unittest

from sys import platform


# noinspection PyUnresolvedReferences
class TestEnv(unittest.TestCase):
    def setUp(self):
        print('Platform: ', platform)
        parent_dir = os.path.dirname(os.path.normpath(os.path.dirname(__file__)))
        sys.path.append(parent_dir)

        # for a in sys.path:
        #     print(a)

    # @unittest.skip("skipping test...")
    def test_snappy_available(self):
        import snappy
        print('snappy import ok...')
        import os
        print('os import ok...')

# suite = unittest.TestLoader().loadTestsFromTestCase(TestEnv)
# unittest.TextTestRunner(verbosity=2).run(suite)
