import unittest
import numpy as np

import sys, os
sys.path.insert(0, os.path.abspath('..'))

from time_integrator import LSERK4
from spatial_discretization.equation import Equation

class CapacitorCharge(Equation):
    """
        Ordinary differential equation example:
            dq(t)/dt = V/R - q(t)/(R*C)
    """
    def __init__(self, R, C, V):
        self.R = R
        self.V = V
        self.C = C
        self.vars = [np.zeros(1)]

    def rhs(self, t):
        return self.V / self.R - self.vars[0] / (self.R * self.C)
    
    def exact_solution(self, t):
        return self.V * self.C * (1 - np.exp(-t/(self.R*self.C)) )


class TestLSERK4(unittest.TestCase):

    def test_single_variable_time_integration(self):
        eq = CapacitorCharge(R=10, C=0.1, V=10)
        integrator = LSERK4(eq)
        
        dt = 0.01
        number_of_steps = 1000
        t = np.linspace(0, dt*number_of_steps, number_of_steps)
        
        integrated = np.zeros((number_of_steps))
        exact      = np.zeros((number_of_steps))
        for n in range(number_of_steps):

            integrated[n] = eq.vars[0]
            integrator.integrate(dt, 1)

            exact[n] = eq.exact_solution(t[n])

        error = np.linalg.norm(integrated - exact)/number_of_steps
        assert(error < 1e-5)


if __name__ == '__main__':
    unittest.main()
