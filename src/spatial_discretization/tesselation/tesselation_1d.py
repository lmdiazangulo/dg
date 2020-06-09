

class Tesselation1D:

    @staticmethod
    def are_valid_opts(opts):
        return \
            opts["dimension"] == 1 and \
            opts["type"]      == "dg" and \
            opts["flux_type"] == "centered" and \
            opts["basis"]     == {"type": "nodal_lgl", "order": 1} and \
            "grid" in opts

    def __init__(self, opts):
        if not Tesselation1D.are_valid_opts(opts):
            raise ValueError("Tesselation parameters are not supported")

        print("TBD") #TODO