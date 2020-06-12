

class Tesselation1D:

    @staticmethod
    def are_valid_opts(opts):
        return \
            opts["dimension"] == 1 and \
            opts["type"]      == "dg" and \
            "flux_type" in opts and \
            opts["basis"]     == {"type": "nodal_lgl", "order": 1} and \
            "grid" in opts


    def __init__(self, opts):
        if not Tesselation1D.are_valid_opts(opts):
            raise ValueError("Tesselation parameters are not supported")

        print("TBD") #TODO


    def field(self, field_type):
        """
            field returns a one dimensional numpy array containing the 
            fields of type field_type
        """
        print("TBD") #TODO
    
    def curl(self, field_type):
        """
            curl returns a one dimensional numpy array containing the 
            discrete curl of fields of type field_type
        """
        print("TBD") #TODO

    def flux(self, field_type):
        """
            flux returns a one dimensional numpy array containing the 
            numerical flux of fields of type field_type
        """
        print("TBD") #TODO
