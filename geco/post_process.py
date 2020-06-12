'''
Solutions class that holds useful information about solutions in the given directory
'''
import geco_utils as gu
import os
import numpy as np

class PostProcess:

    def __init__(self, subdir, resolution=250, radius=50, z=0):
      print(os.getcwd())

      self.subdir = subdir
      self.resolution = resolution
      self.radius = radius
      self.z = z

      mesh, parameters, density_components, self.U, V = gu.GatherFiles(subdir)

      #removes the combined rho.xdmf, we usually arent interested in it
      # might we be though?
      if len(density_components) > 1:
        density_components.remove(density_components[0])

      self.titles_for_plot = [gu.GetTitle(params) for params in parameters]
      self.parameter_dict = [gu.GetParametersDicts(params) for params in parameters]
      self.parameter_str = [gu.GetParametersStrings(params) for params in parameters]
      self.density_arrays = [gu.ToNumpyArray(comp, r_max=self.radius, res=self.resolution) for comp in density_components]
      self.radius_of_support = [gu.GetRadiusSupport(sup, res=self.resolution) for sup in density_components]

      self.radii, self.inv_sqrt_r,self.circular_velocity = [],[],[]
      rotation_curve_output= [gu.RotationCurve(self.U, sup, res=resolution,z=z) for sup in self.radius_of_support]

      for entry in rotation_curve_output:
        self.radii.append(entry[0])
        self.inv_sqrt_r.append(entry[1])
        self.circular_velocity.append(entry[2])


    def compare_with_observed(self,observed_radii,  observed_velocities,density_component=0):
      radius_vector = self.radii[density_component]
      support = self.radius_of_support[density_component]

      max_observed = np.max(observed_radii)

      scaled_radius = radius_vector/support*max_observed

      mltplr, error = gu.ModelToObsVelocity(observed_radii,observed_velocities,
                                 scaled_radius, self.circular_velocity[density_component])

      scaled_velocity = mltplr * self.circular_velocity[density_component]

      return error, scaled_velocity, scaled_radius

    def forward_abel_transform(self):

      transformed = [gu.ForwardAbelTransform(arr) for arr in self.density_arrays]

      return transformed
