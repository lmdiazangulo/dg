import numpy as np
import math
from time_integrator import TimeIntegrator

class LSERK4(TimeIntegrator):
    a = np.array([                             0.0, \
                   -567301805773.0/1357537059087.0, \
                  -2404267990393.0/2016746695238.0, \
                  -3550918686646.0/2091501179385.0, \
                   -1275806237668.0/842570457699.0])
    b = np.array([ 1432997174477.0/9575080441755.0, \
                  5161836677717.0/13612068292357.0, \
                   1720146321549.0/2090206949498.0, \
                   3134564353537.0/4481467310338.0, \
                  2277821191437.0/14882151754819.0])
    c = np.array([                             0.0, \
                   1432997174477.0/9575080441755.0, \
                   2526269341429.0/6820363962896.0, \
                   2006345519317.0/3224310063776.0, \
                   2802321613138.0/2924317926251.0])

    @staticmethod
    def are_valid_options(opts):
        return \
            "type" in opts and opts["type"] == "lserk4" and \
            "final_time" in opts

    def __init__(self, opts, equation):
        
        self.opts = opts
        self.equation = equation
        self.residue = np.zeros(equation.fields.size)

    def integrate(self, dt, number_of_steps):
        TimeIntegrator.start_timer()
        
        t = 0.0
        for n in range(number_of_steps):
            for i in range(5):
                self.residue *= self.a[i]
                self.residue += dt * self.equation.rhs(t)
                self.equation.fields += self.b[i] * self.residue

            t += dt
            TimeIntegrator.print_progress(n, number_of_steps)

        TimeIntegrator.print_cpu_time()