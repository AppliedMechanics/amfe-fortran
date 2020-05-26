from setuptools import Extension, setup
from Cython.Build import cythonize

ext_package = 'amfefortran'
extensions = [Extension('shapefunction', ['amfefortran/shapefunction.pyx'],
                        libraries=['shapefunction', 'gfortran'],
                        library_dirs=['lib'],
                        )
              ]

setup(name='amfefortran',
      version='0.1.6',
      description='Fortran Accelerator Module for the Finite Element Code AMfe',
      author='Christian Meyer',
      author_email='christian.meyer@tum.de',
      url='www.am.mw.tum.de',
      ext_package=ext_package,
      ext_modules=cythonize(extensions),
      packages=['amfefortran'],
      include_package_data=True,
      )
