---
title: 'GECo: A collection of solvers for the self-gravitating Vlasov equations'
tags:
  - Einstein-Vlasov
authors:
  - name: Ellery Ames
    orcid: 0000-0001-9444-585X
    equal-contrib: true
    affiliation: "1"
  - name: Anders Logg
    orcid: 0000-0002-1547-4773
    equal-contrib: true
    affiliation: "2"
affiliations:
 - name: California Polytechnic University, Humboldt
   index: 1
 - name: Chalmers University of Technology
   index: 2
date: 21 September 2022
bibliography: paper.bib
---

# Summary

Gothenburgh Einstein solver Collection (GECo) is a collection of solvers for
stationary self-gravitating collisionless kinetic (Vlasov) matter. The
gravitational interaction may be taken to be either Newtonian or general
relativistic. GECo is focused on the solutions which are axisymmetric, meaning
that the gravitational and matter fields have a rotational symmetry. In this
setting stationary solutions may be generated with the choice of a particular
ansatz function for the Vlasov distribution function. GECo allows users to
easily introduce new ansatz functions and explore the properties of the
resulting stationary solutions.

# Statement of need

In understanding a physical model one usually starts with a simplified setting,
such as by imposing symmetry assumptions. In the case of self-gravitating
kinetic matter, stationary solutions in the spherically symmetric setting are
well understood [@BinneyTremaine:2008; @Andreasson:2011]. However, many of the
physical systems of interest such as accretion disks, galaxies, galaxy clusters,
and so on require models beyond spherical symmetry. When going beyond spherical
symmetry, the coupled and nonlinear PDE systems in high dimensions such as the
self-gravitating Vlasov equations are difficult to investigate analytically, and
numerical approaches are essential to understand behavior of solutions and to
answer questions of physical and mathematical interest. The GECo code started
with the desire to understand properties of stationary and axisymmetric
solutions of the Einstein-Vlasov system. 

# Method and implementation

To construct stationary solutions the code relies on a reduction method in which
the distribution function for the matter is assumed to depend on the position
and momentum phase-space coordinates solely through conserved quantities, such
as the particle energy and angular momentum about the axis of symmetry. With
this ansatz the Einstein--Vlasov or Vlasov--Poisson system (depending on the
gravitational model used) forms a semi-linear integro-differential system of
equations. In GECo, the form of the ansatz is called a `MaterialModel` and
several different choices are implemented as subclasses of the dolfin Expression
class. The semi-linear integro-differential system is solved via a
mass-preserving fixed point scheme using Anderson acceleration [@Walker:2011].
At each step of the fixed point method the linear system of equations is solved
using finite elements implemented with the FEniCS toolkit [@Logg:2012]. The
computational domain is taken to be the half-meridional plane $\{(r,z): r>0, z>0
\}$ in cylindrical coordinates, with a semi-circular outer boundary
\autoref{fig:Solution}. Details of the mathematical formulation and
implementation can be found in [@Ames:2016]

# Functionality

The entrypoint for GECo is a run script written in python. In this file the user
selects the solver class (`EinsteinVlasovSolver` or `VlasovPoisson`) that
specifies the model for the gravitational interaction, a `MaterialModel` to
specify the particular form of the reduction ansatz, and several parameters
related to the model and discretization. Calling the `solve` method within the
script invokes the solver to construct a stationary solution via the fixed point
scheme mentioned above, which has converged to the specified tolerance.
Gravitational fields and matter quantities are saved in XMDF and XML format that
can be consumed by Paraview as well as postprocessing scripts. Multi-component
solutions may be constructed from multiple `MaterialModel`s by combining models
in a weighted sum.

GECo includes several postprocessing routines that:

* Generate additional scalar data not computed during the fixed point iteration. 
* Represent the matter density as well as an ergoregion (if present) in
  $\mathbb{R}^2$ (i.e. reflected about the reflection plane and symmetry axis),
  as shown in \autoref{fig:2Ddensity}.
* Represent the matter density as well as an ergoregion (if present) as a volume in $\mathbb{R}^3$, facilitating visualization of contours, as shown in \autoref{fig:3Ddensity}. 
* Represent the density as a three-dimensional point cloud,as shown in \autoref{fig:PointCloud}. 
* Compute the Kretschmann curvature scalar.

![Torus spatial density on quarter plane computational domain \label{fig:Solution}](./figures/density_computational_domain.png){ width=25% }
![2D Density\label{fig:2Ddensity}](./figures/density_2d_density.png){ width=25% }
![3D Density\label{fig:3Ddensity}](./figures/density_3d_density.png){ width=25% }
![Pointcloud\label{fig:PointCloud}](./figures/density_3d_pointcloud.png){ width=25% }

# Documentation

The documentation for GECo is published on the
[GECo GitLab pages](https://gitlab.com/alogg/geco).

# Limitations and future work

We briefly list a few directions of interest for future work. 

- GECo currently uses a uniform mesh. However, in axisymmetry (unlike spherical
symmetry) the solution is not uniquely defined outside the support of the matter
and asymptotically flat boundary conditions must be applied sufficiently far
from the matter. An adaptive mesh refinement algorithm was developed and used in
[@Ames:2019] to investigate properties of extreme rotating toroidal solutions.
It remains however to integrate such an adaptive mesh refinement scheme into the
core of GECo. 
- Currently the particles only interact via the gravitational field generated by
the particle distribution. An exciting area at the frontier of astrophysics
currently is the study of accretion disks, where both central black holes and
electromagnetic fields play important roles. To lay groundwork for this area in
fundamental relativity, it is thus highly desirable to extend GECo to the
Einstein-Vlasov-Maxwell system and allow the inclusion of central black holes.
- While multi-species solutions can be generated in which the different species
follow different distribution ansatzes, the particle properties are otherwise
taken to be the same. Astrophysical systems however often consist of
particle-like entities with very different properties (such as stars and dust).
We thus propose to allow different particle species to have different particle
properties such as mass and charge.

# References
