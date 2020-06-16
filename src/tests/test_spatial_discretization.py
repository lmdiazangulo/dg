import unittest
import numpy as np

import sys, os
sys.path.insert(0, os.path.abspath('..'))

from spatial_discretization.tesselation_1d import Tesselation1D

tesselation_request = \
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

    order = tesselation_request["basis"]["order"]
    number_of_nodes = order + 1
    box   = tesselation_request["grid"]["box"]
    step  = tesselation_request["grid"]["steps"]
    number_of_elements = (box[1] - box[0]) / step
    num_vars = number_of_nodes * number_of_elements

    def test_basic_features(self):
        try:
            t = Tesselation1D(tesselation_request)
        except:
            raise ValueError("Unable to build")
        
    def test_sizes(self):
        t = Tesselation1D(tesselation_request)

        self.assertEqual(t.field("E").size, self.num_vars)
        self.assertEqual(t.field("E").size, t.curl( "E").size)
        self.assertEqual(t.field("E").size, t.flux( "E").size)

        self.assertEqual(t.field("E").size, t.field("H").size)
        self.assertEqual(t.field("H").size, t.curl( "H").size)
        self.assertEqual(t.field("H").size, t.flux( "H").size)

    def test_smallest_distance(self):
        t = Tesselation1D(tesselation_request)
        min_dist = self.step/(self.number_of_nodes-1)
        self.assertEqual(t.get_smallest_distance(), min_dist)

if __name__ == '__main__':
    unittest.main()
