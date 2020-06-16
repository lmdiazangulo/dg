from .tesselation_1d import Tesselation1D
import numpy as np

class SpatialDiscretization:
    """
        SpatialDiscretization manages a group of tesselations managing the 
        variables in a system of equations of the form 
            d{f} / dt = [A] {f}(t) + {g}(t)
        with {f} being a vector of fields, [A] is a matrix of coefficients, and
        {g} is a forcing function. 
    """


    @staticmethod
    def are_valid_opts(opts):
        return ("tesselations" in opts)

    def __init__(self, opts):
        
        if not SpatialDiscretization.are_valid_opts(opts):
            raise ValueError( \
                "Spatial Discretization has invalid options")

        self.tesselations = []
        for tess_request in opts["tesselations"]:
            if Tesselation1D.are_valid_opts(tess_request):
                self.tesselations.append(Tesselation1D(tess_request))
            else:
                raise ValueError("Tesselation options not supported")

        self.vars = []
        for tess in self.tesselations:
            for field_type in tess_request["fields"]:
                self.vars.append(tess.field(field_type))

    def get_smallest_distance(self) -> float:
        res = float("inf")
        for tess in self.tesselations:
            aux = tess.get_smallest_distance()
            if aux < res:
                res = aux
        return res

if __name__ == '__main__':
    import doctest
    doctest.testmod()