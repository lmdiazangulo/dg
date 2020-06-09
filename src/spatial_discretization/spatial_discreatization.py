from tesselation.tesselation_1d import Tesselation1D

class SpatialDiscretization:
    """
        SpatialDiscretization manages a group of tesselations.
    """


    @staticmethod
    def are_valid_opts(opts):
        return ("tesselations" in opts)

    def __init__(self, opts):

        if not SpatialDiscretization.are_valid_opts(opts):
            raise ValueError("DG has invalid options")

        self.tesselations = []
        for tesselation in opts["tesselations"]:
            if Tesselation1D.are_valid_opts(tesselation):
                self.tesselations.append(Tesselation1D(tesselation))
            else:
                raise ValueError("Tesselation options not supported")

    def init_fields(self):
        print("TBD") #TODO

    def compute_rhs(self):
        print("TBD") #TODO

if __name__ == '__main__':
    import doctest
    doctest.testmod()