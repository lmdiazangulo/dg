import time
import math
import numpy as np



class TimeIntegrator:
    """
        Time integrators perform the evolution of a set of equations of the form
            d{f} / dt = [A] {f}(t) + {g}(t)
        with {f} being a vector of fields, [A] is a matrix of coefficients, and
        {g} is a forcing function. 
    """
    __timeStepPrint = 10

    def start_timer(self):
        self.tic = time.time()

    def print_progress(self, n, number_of_steps):
        if n % self.__timeStepPrint == 0 or n+1 == number_of_steps:
            remaining = (time.time() - self.tic) *  (number_of_steps-n) / (n+1)
            min = math.floor(remaining / 60.0)
            sec = remaining % 60.0
            print("  Step: %6d of %6d. Remaining: %2.0f:%02.0f"% \
                (n, number_of_steps-1, min, sec))


    def print_cpu_time(self):
        print("    CPU Time: %f [s]" % (time.time() - self.tic))


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

    def __init__(self, equation):        
        self.equation = equation
        self.residue = \
            list(map(lambda x: np.zeros(x.shape), self.equation.vars))

    def integrate(self, dt, number_of_steps):
        t = 0.0
        for _ in range(number_of_steps):
            for i in range(5):
                for res in self.residue: 
                    res *= self.a[i]
                    res += dt * self.equation.rhs(t)
                    
                for j in range(len(self.equation.vars)):
                    self.equation.vars[j] += self.b[i] * self.residue[j]

            t += dt