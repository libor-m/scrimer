#!/usr/bin/env python

from distutils.core import setup
import glob

with open('README') as f:
    long_description = f.read()

with open('LICENSE.txt') as f:
    license = f.read()

setup(name='scrimer',
      version='1.0',
      description='Primer designing pipeline',
      author='Libor Morkovsky',
      author_email='morkovsk@natur.cuni.cz',
      url='https://github.com/libor-m/scrimer',
      py_modules=[f.replace('.py', '') for f in glob.glob('modules/*.py')],
      scripts=glob.glob('scripts/*.py'),
      requires=['PyVCF', 'pybedtools', 'pysam'],
      data_files=[('demo', glob.glob('demo/*'))],
      long_description=long_description,
      # https://pypi.python.org/pypi?%3Aaction=list_classifiers
      classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)',
          'Natural Language :: English',
          'Operating System :: POSIX :: Linux',
          'Programming Language :: Python :: 2.7',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
      ],
      license=license,
     )
