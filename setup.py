#!/usr/bin/env python

try:
    from setuptools import setup
    have_setuptools = True
except ImportError:
    from distutils.core import setup
    have_setuptools = False

import glob

with open('README') as f:
    long_description = f.read()

with open('LICENSE.txt') as f:
    license = f.read()

setup(name='scrimer',
      version='1.1',
      description='Primer designing pipeline',
      author='Libor Morkovsky',
      author_email='morkovsk@natur.cuni.cz',
      url='https://github.com/libor-m/scrimer',
      py_modules=[f.replace('.py', '') for f in glob.glob('scrimer/*.py')],
      scripts=glob.glob('scripts/*.py'),
      # http://stackoverflow.com/questions/6947988/when-to-use-pip-requirements-file-versus-install-requires-in-setup-py
      #install_requires=['PyVCF (>=0.6)', 'pybedtools (>=0.6)', 'pysam (>=0.7)'],
      install_requires=['PyVCF', 'pybedtools', 'pysam'],
      # requires=['PyVCF (>=0.6)', 'pybedtools (>=0.6)', 'pysam (>=0.7)'],
      data_files=[('demo', glob.glob('demo/*'))],
      long_description=long_description,
      # https://pypi.python.org/pypi?%3Aaction=list_classifiers
      classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)',
          'Natural Language :: English',
          'Operating System :: POSIX',
          'Operating System :: POSIX :: Linux',
          'Programming Language :: Python',
          'Programming Language :: Python :: 2.7',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
      ],
      license=license,
      entry_points="""
      [vcf.filters]
      contrast-samples = scrimer.pyvcf_filters:DistinguishingVariants
      """,
     )
