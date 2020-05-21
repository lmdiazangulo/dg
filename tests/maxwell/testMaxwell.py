import unittest
import json

import maxwell

class TestMaxwell(unittest.TestCase):

    default_cases_path = "cases/"

    def test_maxwell_1D(self):
        case_file_name = "cavity_1d.case.json"
        case = json.load(open(case_file_name))
        
        solver = Maxwell(case)
        solver.solve()
        
        viewer = Viewer(solver.getProbes()[0])

if __name__ == '__main__':
    unittest.main()
