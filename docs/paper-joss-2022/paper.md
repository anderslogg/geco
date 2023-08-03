---
title: 'GECo: A collection of solvers for the Einstein-Vlasov equations'
tags:
  - Einstein-Vlasov
authors:
  - name: Ellery Ames
	orcid: FIXME
    equal-contrib: true
    affiliation: "1"
  - name: Anders Logg
    orcid: 0000-0002-1547-4773
    equal-contrib: true
    affiliation: "2"
affiliations:
 - name: FIXME
   index: 1
 - name: Chalmers University of Technology
   index: 2
date: 21 September 2022
bibliography: paper.bib
---

# Summary

Gothenburgh Einstein solver Collection (GECo) is a collection of solvers for stationary self-gravitating collisionless kinetic matter. The gravitational interaction may be either Newtonian or general relativistic. The code makes use a reduction method, finite elements, and adaptive mesh refinement.  

- High level answer to "what is GECO?"
- 

# Statement of need

FIXME
- unlike spherical symmetry, in axisymmetry the solution outside of the matter distribution is not canonical, and far-field boundary conditions must be applied.

We can cite stuff like this: [@amesAxisymmetricStationarySolutions2016] and [@amesCosmicStringBlack2019].

# Functionality

FIXME (examples of where and for what GECO has been used.)

# Method and implementation

FIXME
- reduction based on method of [NAME] by which the distribution is a function of the conserved quantities. This results in a semilinear system of equations for the gravitational field.
- FEM (say something about implementation and limitations of the scheme (eg which FE are used, and does the code provide hte user with any flexibility.))
- Mass-preserving fixed point scheme to solve the nonlinear system of equations with Anderson acceleration. 
- Mesh refinement 

# Documentation

The documentation for GECo is published on the
[GECo GitLab pages](https://gitlab.com/alogg/geco).

# Limitations and future work

FIXME
- extension to spherical symmetry for completeness
- allow different particle species to have different properties such as mass
- extension to the Einstein-Vlasov-Maxwell system. 

# Acknowledgements

FIXME

# References
