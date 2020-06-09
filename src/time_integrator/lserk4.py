import numpy as np
import time
import math

rk4a = np.array([                             0.0, \
                  -567301805773.0/1357537059087.0, \
                 -2404267990393.0/2016746695238.0, \
                 -3550918686646.0/2091501179385.0, \
                 -1275806237668.0/842570457699.0])
rk4b = np.array([ 1432997174477.0/9575080441755.0, \
                 5161836677717.0/13612068292357.0, \
                  1720146321549.0/2090206949498.0, \
                  3134564353537.0/4481467310338.0, \
                 2277821191437.0/14882151754819.0])
rk4c = np.array([                             0.0, \
                  1432997174477.0/9575080441755.0, \
                  2526269341429.0/6820363962896.0, \
                  2006345519317.0/3224310063776.0, \
                  2802321613138.0/2924317926251.0])

class LSERK4:

    __timeStepPrint = 10

    @staticmethod
    def are_valid_options(opts):
        return \
            "type" in opts and opts["type"] == "lserk4" and \
            "cfl" in opts and \
            "final_time" in opts


    def __init__(self, case, spatial_discretization):

        self.opts = case["solver"]["time_integrator"]
        if not LSERK4.are_valid_options(self.opts):
            raise ValueError("LSERK4 are not valid")

        self.spatial_discretization = spatial_discretization
        self.residue = spatial_discretization.init_fields()


    def maximum_time_step(self):
        return self.opts["cfl"] * \
               self.spatial_discretization.get_smallest_space_step()


    def integrate(self):
        dt = self.maximum_time_step()

        if "final_time" in self.opts:
            number_of_steps = int (self.opts["final_time"] / dt)
        elif "number_of_steps" in self.opts:
            number_of_steps = self.opts["number_of_steps"]
        else:
            raise ValueError("Ending condition is not defined")

        tic = time.time()
        for n in range(number_of_steps):
            for i in range(5):
                self.residue *= rk4a[i]
                self.residue += dt * self.spatial_discretization.compute_rhs()
                
                self.spatial_discretization.fields += rk4b[i] * self.residue

            if n % self.__timeStepPrint == 0 or n+1 == number_of_steps:
                remaining = (time.time() - tic) *  (number_of_steps-n) / (n+1)
                min = math.floor(remaining / 60.0)
                sec = remaining % 60.0
                print("  Step: %6d of %6d. Remaining: %2.0f:%02.0f"% \
                    (n, number_of_steps-1, min, sec))
        
        print("    CPU Time: %f [s]" % (time.time() - tic))