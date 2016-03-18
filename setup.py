from setuptools import setup
from os.path import join as pj

setup(name='GECo',
      version='0.1dev',
      packages=['geco'],
      scripts=[pj('bin', 'geco'),
               pj('bin', 'geco-postprocess-volume'),
               pj('bin', 'geco-postprocess-moments'),
               pj('bin', 'geco-postprocess-pointcloud'),
               pj('bin', 'geco-postprocess-ergoregion'),               
               pj('bin', 'geco-dict-plot')],
      license='To be decided...',
      long_description=open('README.rst').read(),
      include_package_data=True)
