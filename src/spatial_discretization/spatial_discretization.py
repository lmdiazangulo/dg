from .tesselation_1d import Tesselation1D
import numpy as np

class SpatialDiscretization:
    """
        SpatialDiscretization manages a group of tesselations. It can be also 
        regarded as a system of equations of the form 
            d{f} / dt = [A] {f}(t) + {g}(t)
        with {f} being a vector of fields, [A] is a matrix of coefficients, and
        {g} is a forcing function. 
    """


    @staticmethod
    def are_valid_opts(opts):
        return ("tesselations" in opts)

    def __init__(self, opts):
        self.vars = [np.array(1)] #TODO

        if not SpatialDiscretization.are_valid_opts(opts):
            raise ValueError("DG has invalid options")

        self.tesselations = []
        for tesselation_request in opts["tesselations"]:
            if Tesselation1D.are_valid_opts(tesselation_request):
                self.tesselations.append(Tesselation1D(tesselation_request))
            else:
                raise ValueError("Tesselation options not supported")

    def init_fields(self):
        print("TBD") #TODO

    def rhs(self, t=None):
        print("TBD") #TODO

if __name__ == '__main__':
    import doctest
    doctest.testmod()