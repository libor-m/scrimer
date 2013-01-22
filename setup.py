#!/usr/bin/env python

from distutils.core import setup
import glob

setup(name='Scrimer',
      version='1.0',
      description='Primer designing pipeline',
      author='Libor Morkovsky',
      author_email='morkovsk@natur.cuni.cz',
      url='https://github.com/libor-m/scrimer',
      py_modules=[f.replace('.py', '') for f in glob.glob('*.py')],
      requires=['PyVCF', 'pybedtools', 'pysam'],
     )
