import unittest
import numpy as np

import sys, os
sys.path.insert(0, os.path.abspath('..'))

from time_integrator import *
from time_integrator.lserk4 import *

class ODE_capacitor_charge:
    """
        Ordinary differential equation example:
            dq(t)/dt = V/R - q(t)/(R*C)
    """
    def __init__(self, R, C, V):
        self.R = R
        self.V = V
        self.C = C
        self.vars = np.zeros(1)

    def rhs(self, t):
        return self.V / self.R - self.vars / (self.R * self.C)
    
    def exact_solution(self, t):
        return self.V * self.C * (1 - np.exp(-t/(self.R*self.C)) )


class TestLSERK4(unittest.TestCase):

    def test_single_variable_time_integration(self):
        eq = ODE_capacitor_charge(R=10, C=0.1, V=10)
        integrator = LSERK4(eq)
        dt = 0.01
        number_of_steps = 1000
        t = np.linspace(0, dt*number_of_steps)
        integrated_var = np.zeros((1,number_of_steps))
        exact_var      = np.zeros((1,number_of_steps))
        for n in range(number_of_steps):
            integrated_var[n] = integrator.integrate(dt, 1)
            exact_var[n] = eq.exact_solution(t[n])

        

if __name__ == '__main__':
    unittest.main()
