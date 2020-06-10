import time
import math

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

