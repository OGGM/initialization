from setuptools import setup

req_packages = ['numpy',
                'scipy',
                'pandas',
                'matplotlib>=2.0.0',
                'shapely',
                'requests',
                'configobj',
                'netcdf4',
                'xarray',
                'salem',
                'oggm']

setup(
	  # Project info
	  name='initialization',
      version='1.0.0',
      description='initialization of past glacier states',
	  # Github repostiory link
      url='http://github.com/OGGM/initialization',
	  # Authors details
      author='OGGM Developers',
      author_email='jeis@uni-bremen.de',
	  # License
      license='BSD-3-Clause',
	  # What does the project relate to?
	  keywords=['geosciences', 'glaciers', 'initialization','reconstruction'],
      packages=['initialization'],
	  # Python 3 only 
	  python_requires='>=3.5',
	  # Install all dependencies --> same as for OGGM
	  install_requires=req_packages,
      zip_safe=False)

