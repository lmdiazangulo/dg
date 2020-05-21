import dg
import lserk4

class Maxwell:

    def __init__(self, case):

        if case["solver"]["type"] != "maxwell":
            raise ValueError("Invalid solver type")
        
        # Creates spatial discretization.
        if case["solver"]["spatial_discretization"]["type"] == "dg":
            self.spatial_discretization = dg.DG(case)
        else:
            raise ValueError("Invalid spatial discretization")

        # Creates time integrator.
        if case["solver"]["time_integrator"]["type"] == "lserk4":
            self.time_integrator = \
                lserk4.LSERK4(case, self.spatial_discretization)
        else:
            raise ValueError("Invalid time integrator")

