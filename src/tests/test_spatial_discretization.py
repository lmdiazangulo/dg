import unittest
import numpy as np

import sys, os
sys.path.insert(0, os.path.abspath('..'))

from spatial_discretization.tesselation_1d import Tesselation1D

str_tesselation_1d = \
    {
        "dimension": 1,
        "type": "dg",
        "fields": [
            "E",
            "H"
        ],
        "flux_type": "upwind",
        "basis": {
            "type": "nodal_lgl",
            "order": 1
        },
        "grid": {
            "box": [0.0 , 1.0],
            "steps": 0.025,
            "bounds": ["pec", "pec"]
        }
    }

class TestTesselation1D(unittest.TestCase):

    def test_basic_features(self):
        try:
            t = Tesselation1D(str_tesselation_1d)
        except:
            raise ValueError("Unable to build")
        
    def test_sizes(self):
        t = Tesselation1D(str_tesselation_1d)

        order = str_tesselation_1d["basis"]["order"]
        box   = str_tesselation_1d["grid"]["box"]
        step  = str_tesselation_1d["grid"]["steps"]
        number_of_elements = (box[1] - box[0]) / step
        num_vars = (order + 1) * number_of_elements

        self.assertEqual(t.field("E").size, num_vars)
        self.assertEqual(t.field("E").size, t.curl( "E").size)
        self.assertEqual(t.field("E").size, t.flux( "E").size)

        self.assertEqual(t.field("E").size, t.field("H").size)
        self.assertEqual(t.field("H").size, t.curl( "H").size)
        self.assertEqual(t.field("H").size, t.flux( "H").size)



if __name__ == '__main__':
    unittest.main()
