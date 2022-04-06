"""
-----------
solution.py
-----------
This module defines a solution class.

Copyright 2019 Anders Logg, Ellery Ames, Håkan Andréasson

This file is part of GECo.
GECo is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

GECo is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with GECo. If not, see <https://www.gnu.org/licenses/>.
"""


class Solution:
    def __init__(self, NU, BB, MU, WW, RHO, P00, P11, P33, P03, RMD, data):

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
        self.RMD = RMD  # rest mass density
        self.data = data

        # Extract function space and mesh
        self.V = NU.function_space()
        self.mesh = self.V.mesh()

    def geometry_fields_list(self):
        "Return geometry fields in list"

        return [self.NU, self.BB, self.MU, self.WW]
