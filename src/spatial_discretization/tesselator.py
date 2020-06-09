from tesselation.tesselation_1d import Tesselation_1D

class Tesselator:
    
    @staticmethod
    def are_valid_opts(opts):
        return ("tesselations" in opts)

    def __init__(self, opts):

        if not Tesselator.are_valid_options(opts):
            raise ValueError("DG has invalid options")

        self.tesselations = []
        for tesselation in opts["tesselations"]:
            if Tesselation_1D.are_valid_opts(tesselation):
                self.tesselations.append(Tesselation_1D(tesselation))
            else:
                raise ValueError("Tesselation options not supported")

if __name__ == '__main__':
    import doctest
    doctest.testmod()