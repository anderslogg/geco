"""
This module implements the Anderson algorithm for acceleration of
general fixed-point iterations. Often, it is enough to set m = 3.

Usage:

anderson = Anderson(3)
x = x0
while not_converged:

  gx = g(x)
  anderson.update_solution(gx)

  x = anderson.next()
  
  fx = f(x)
  anderson.update_residual(fx)
"""

import numpy.linalg

class Anderson:

    def __init__(self, m, x0):
        "Initialize acceleration with depth m and initial value x0"
        self._m = m
        self._X = [x0.copy()]
        self._G = []

    def update(self, gx):

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
            raise RuntimeError, "Number of vectors don't match."
        
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
