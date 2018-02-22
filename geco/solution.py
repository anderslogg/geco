class Solution:

    def __init__(self, NU, BB, MU, WW,
                     RHO, P00, P11, P33, P03, RMD,
                     data):

        # Store fields
        self.NU = NU
        self.BB = BB
        self.MU = MU
        self.WW = WW
        self.RHO = RHO
        self.P00 = P00
        self.P11 = P11
        self.P33 = P33
        self.P03 = P03
        self.RMD = RMD # rest mass density
        self.data = data

        # Extract function space and mesh
        self.V= NU.function_space()
        self.mesh = self.V.mesh()


    def geometry_fields_list(self):
        " Return geometry fields in list "
        
        return [self.NU, self.BB, self.MU, self.WW] 
