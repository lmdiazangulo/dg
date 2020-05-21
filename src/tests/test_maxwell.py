import unittest
import json

import sys, os
sys.path.insert(0, os.path.abspath('..'))

import maxwell

class TestMaxwell(unittest.TestCase):

    default_cases_path = "src/tests/maxwell_cases/"

    def test_basic_input_file(self):
        case_file_name = "cavity_1d.case.json"
        case = json.load(open(self.default_cases_path + case_file_name))
        try:
            solver = maxwell.Maxwell(case)
        except:
            self.fail("Problem parsing file")
        

if __name__ == '__main__':
    unittest.main()
