
class Equation:
    """
        Equation represents a system of equations of the form
            d{f} / dt = [A] {f}(t) + {g}(t)
        with {f} being a vector of fields, [A] is a matrix of coefficients, and
        {g} is a forcing function. 
    """

    def rhs(self, t=None):
        raise NotImplementedError