from spatial_discretization.spatial_discretization import SpatialDiscretization
from spatial_discretization.equation import Equation
from time_integrator import LSERK4

class MaxwellEquations(SpatialDiscretization, Equation):

    def __init__(self, opts):
        SpatialDiscretization.__init__(self, opts)


    def rhs(self, t=None):
        res = []
        for tess in self.tesselations:
            res.append( - tess.curl("H") + tess.flux("E") )
            res.append( - tess.curl("E") + tess.flux("H") )


class Maxwell:
    """
        Maxwell is an EM Solver which uses a spatial discretization and a time 
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

        # Creates spatial discretization.
        try:
            self.spatial_discretization = \
                MaxwellEquations(case["solver"]["spatial_discretization"])
        except:
            raise ValueError("Invalid spatial discretization")

        # Initializes probes.
        for probe in case["electromagnetics"]["probes"]:
            self._init_probe(probe)

        # Sets sources if exist.
        for source in case["electromagnetics"]["sources"]:
            self._init_source(source)

        # Creates time integrator.
        if case["solver"]["time_integrator"]["type"] == "lserk4":
            self.time_integrator = LSERK4(self.spatial_discretization)
        else:
            raise ValueError("Invalid time integrator type")


    def _global_max_time_step(self, time_opts):
        if not "cfl" in time_opts:
            cfl = 1.0
        else:
            cfl = time_opts["cfl"] 
        return cfl * self.spatial_discretization.get_smallest_distance()


    def _init_probe(self, probe_request):
        print("TBD") #TODO


    def _init_source(self, source_request):
        print("TBD") #TODO


    def solve(self):
        time_opts = self.case["solver"]["time_integrator"]
        dt = self._global_max_time_step(time_opts)

        if "final_time" in time_opts:
            number_of_steps = int (time_opts["final_time"] / dt)
        elif "number_of_steps" in time_opts:
            number_of_steps = time_opts["number_of_steps"]
        else:
            raise ValueError("Ending condition is not defined")
        
        t = 0.0
        for n in range(number_of_steps):
            self.time_integrator.integrate(dt, 1, t)
            t += dt
            self._process_probes(t)


    def _process_probes(self, t):
        print("TBD") #TODO




