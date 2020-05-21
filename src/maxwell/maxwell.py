import dg

class Maxwell:

    def __init__(self, case):

        if case["solver"]["type"] != "maxwell":
            raise ValueError("Invalid solver type")
        
        dg = DiscontinuousGalerkin(case)
        