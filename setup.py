from setuptools import setup
from os.path import join as pj
from os import listdir

cppfiles = [pj('geco', 'cppcode', f) for f in listdir('geco/cppcode')]
meshfiles = [pj('geco', 'meshes', f) for f in listdir('geco/meshes')]

setup(name='GECo',
      version='0.1dev',
      packages=['geco'],
      scripts=[pj('bin', 'geco'),
               pj('bin', 'geco-postprocess-data'),
               pj('bin', 'geco-postprocess-print'),
               pj('bin', 'geco-postprocess-volume'),
               pj('bin', 'geco-postprocess-torus'),
               pj('bin', 'geco-postprocess-pointcloud'),
               pj('bin', 'geco-postprocess-moments'),
               pj('bin', 'geco-postprocess-ergoregion'),
               pj('bin', 'geco-dict-plot')],
      data_files=[(pj('geco','cppcode'), cppfiles),
                  (pj('geco','meshes'), meshfiles)],
      license='To be decided...',
      long_description=open('README.rst').read(),
      include_package_data=True)
