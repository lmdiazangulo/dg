import unittest
import numpy as np

import sys, os
sys.path.insert(0, os.path.abspath('..'))

from spatial_discretization.tesselation_1d import Tesselation1D

str_tesselation_1d = \
    {
        "dimension": 1,
        "type": "dg",
        "flux_type": "centered",
        "basis": {
            "type": "nodal_lgl",
            "order": 1
        },
        "grid": [
            {
                "elem_id": 0,
                "steps": 0.025,
                "bounds": ["pec", "pec"]
            }
        ]
    }

class TestTesselation1D(unittest.TestCase):

    def test_use_simple_ctor(self):
        try:
            t = Tesselation1D(str_tesselation_1d)
        except:
            raise ValueError("Unable to build")


if __name__ == '__main__':
    unittest.main()
