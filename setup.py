#!/usr/bin/env python

from distutils.core import setup

setup(name = 'sv_tools',
      version = '0.1',
      description = 'Assorted tools for working with structural variants',
      author = 'Leon Di Stefano',
      author_email = 'distefano.l@wehi.edu.au',
      # url='',
      packages = ['sv_tools'],
      install_requires = ['matplotlib', 'networkx',
                          'pandas', 'numpy', 'itertools', 're', 'scipy']
     )
