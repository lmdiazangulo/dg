from spatial_discretization.tesselator import Tesselator
from time_integrator.lserk4 import LSERK4

class Maxwell:
    """
        Maxwell is a EM Solver whichs uses a spatial discretization and a time 
        integrator.
    """

    @staticmethod
    def are_valid_inputs(case):
        return \
            Maxwell.are_valid_opts(case["solver"]) and \
            "geometry" in case and \
            "electromagnetics" in case

    @staticmethod
    def are_valid_opts(opts):
        return \
            opts["type"] == "maxwell" \
            "spatial_discretization" in opts and \
            "time_integrator" in opts

    def __init__(self, case):

        if not Maxwell.are_valid_inputs(case):
            raise ValueError("Invalid input case")
        
        # Creates spatial discretization
        try:
            tesselator = Tesselator(case["solver"]["spatial_discretization"])
            self.spatial_discretization = tesselator.discretize(case)
        except:
            raise ValueError("Invalid spatial discretization")

        # Creates time integrator.
        if case["solver"]["time_integrator"]["type"] == "lserk4":
            self.time_integrator = LSERK4(case, self.spatial_discretization)
        else:
            raise ValueError("Invalid time integrator")

