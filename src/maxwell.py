from spatial_discretization.spatial_discretization import SpatialDiscretization
from time_integrator import LSERK4

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
            opts["type"] == "maxwell" and \
            "spatial_discretization" in opts and \
            "time_integrator" in opts


    def __init__(self, case):
        if not Maxwell.are_valid_inputs(case):
            raise ValueError("Invalid input case")
        
        self.case = case

        # Creates spatial discretization
        try:
            self.spatial_discretization = \
                SpatialDiscretization(case["solver"]["spatial_discretization"])
        except:
            raise ValueError("Invalid spatial discretization")

        # Creates time integrator.
        if case["solver"]["time_integrator"]["type"] == "lserk4":
            self.time_integrator = LSERK4(self.spatial_discretization)
        else:
            raise ValueError("Invalid time integrator type")

    @staticmethod
    def global_max_time_step(time_opts, spatial_discretization):
        if not "cfl" in time_opts:
            cfl = 1.0
        else:
            cfl = time_opts["cfl"] 
        return cfl * spatial_discretization.get_smallest_distance()

    def solve(self):
        time_opts = self.case["solver"]["time_integrator"]
        dt = self.global_max_time_step(time_opts, self.spatial_discretization)

        if "final_time" in time_opts:
            number_of_steps = int (time_opts["final_time"] / dt)
        elif "number_of_steps" in time_opts:
            number_of_steps = time_opts["number_of_steps"]
        else:
            raise ValueError("Ending condition is not defined")

        self.time_integrator.integrate(dt, number_of_steps)


