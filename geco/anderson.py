"""
This module implements the Anderson algorithm for acceleration of
general fixed-point iterations. Often, it is enough to set m = 3.

  Solve x = g(x) with Anderson accelerated fixed-point iteration

  anderson = Anderson(3, x0)
  x = x0
  while not_converged:
    gx = g(x)
    x = anderson.update(gx)
"""

import numpy.linalg

class Anderson:

    def __init__(self, m, x0):
        "Initialize acceleration with depth m and initial value x0"
        self._m = m
        self._G = []
        if isinstance(x0, (list, tuple)):
            self._X = [[x0i.copy() for x0i in x0]]
        else:
            self._X = [x0.copy()]

    def update(self, gx):
        "Compute next iteration"

        # Special case: depth = 0
        if self._m == 0: return gx

        # Store new g(x) value
        self._G = self._G[-self._m:] + [gx.copy()]

        # Extract data
        X = self._X
        G = self._G
        mk = len(X) - 1

        # Wait until we have enough values
        if len(self._X) < 2:
            self._X = self._X[-self._m:] + [gx.copy()]
            return gx

        # Check that data makes sense
        if len(X) != len(G):
            raise RuntimeError("Number of vectors don't match.")
        
        # Compute vectors q_j = F_0 - F_j
        Q = [G[0].copy() for j in range(mk)]
        for j in range(mk):
            Q[j].axpy(-1.0, X[0])
            Q[j].axpy(-1.0, G[j + 1])
            Q[j].axpy(+1.0, X[j + 1])

        # Compute normal equations
        A = numpy.zeros((mk, mk))
        b = numpy.zeros((mk, 1))
        for i in range(mk):
            for j in range(mk):
                A[i, j] = Q[i].inner(Q[j])
            b[i] = Q[i].inner(G[0]) - Q[i].inner(X[0])
            
        # Solve linear system
        alpha = numpy.linalg.solve(A, b)[:, 0]

        # Compute weighted sum of iterates
        x = gx.copy()
        x.zero()
        x.axpy(1.0 - sum(alpha), G[0])
        for j in range(mk):
            x.axpy(alpha[j], G[j + 1])

        # Store new x value
        self._X = self._X[-self._m:] + [x.copy()]

        return x

    def update_system(self, gx):
        """Compute next iteration (system version). For this version
        the initial value and gx should be lists of vectors."""

        # Special case: depth = 0
        if self._m == 0: return gx

        # Store new g(x) value
        self._G = self._G[-self._m:] + [[gI.copy() for gI in gx]]

        # Extract data
        X = self._X
        G = self._G
        mk = len(X) - 1

        # Wait until we have enough values
        if len(self._X) < 2:
            self._X = self._X[-self._m:] + [[gI.copy() for gI in gx]]
            return gx

        # Check that data makes sense
        if len(X) != len(G):
            raise RuntimeError("Number of vectors don't match.")
        if len(X[0]) != len(G[0]):
            raise RuntimeError("Size of system does not match.")
        
        # Get size of system
        N = len(X[0])
        
        # Compute vectors q_j = F_0 - F_j
        Q = [[gI.copy() for gI in G[0]] for j in range(mk)]
        for I in range(N):
            for j in range(mk):
                Q[j][I].axpy(-1.0, X[0][I])  
                Q[j][I].axpy(-1.0, G[j + 1][I]) 
                Q[j][I].axpy(+1.0, X[j + 1][I])

        # Compute normal equations
        A = numpy.zeros((mk, mk))
        b = numpy.zeros((mk, 1))
        for I in range(N):
            for i in range(mk):
                for j in range(mk):
                    A[i, j] += Q[i][I].inner(Q[j][I])
                b[i] += Q[i][I].inner(G[0][I]) - Q[i][I].inner(X[0][I])
            
        # Solve linear system
        alpha = numpy.linalg.solve(A, b)[:, 0]

        # Compute weighted sum of iterates
        x = [gI.copy() for gI in gx]
        for I in range(N):
            x[I].zero()
            x[I].axpy(1.0 - sum(alpha), G[0][I])
            for j in range(mk):
                x[I].axpy(alpha[j], G[j + 1][I])

        # Store new x value
        self._X = self._X[-self._m:] + [[xI.copy() for xI in x]]

        return x
