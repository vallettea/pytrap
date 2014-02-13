#!/usr/bin/env python
# -*- coding: utf-8 -*-

import distutils.core
import sys


error_msg = "I'm sorry.  This package is for Python 2.5 and higher only."
try:
    if sys.version_info[:2] < (2, 5):
        print >> sys.stderr, error_msg
        sys.exit(1)
except AttributeError:  # sys.version_info was introduced in Python 2.0
    print >> sys.stderr, error_msg
    sys.exit(1)


distutils.core.setup(
    name='pytrap',
    version='2.0', 
    author='Alexandre Vallette',
    author_email='vallettea@gmail.com',
    url='http://packages.python.org/pytrap/',
      
    license='''\
This software can be used under one of the following two licenses: \
(1) The BSD license. \
(2) Any other license, as long as it is obtained from the original \
author.''',
      
    description=('Utilities for physicists using the '
                 ' Electrostatic Ion Beam Trap'),
    
    long_description=u'''\
====================================
Welcome to the PyTrap package
====================================

This package is destined to physicists using or studying the
(`Electron Ion Beam Trap <http://pra.aps.org/abstract/PRA/v55/i3/pR1577_1>`_).
The main features are:
   - a fast (analytic) calculation of the **potential inside the trap** depending on the set of potentials
   - a **stability map** indicating which potentials lead to stable trapping
   - trajectories simulations.
Most of the proofs can be found in this `article`_.

Examples are the best way to dive into PyTrap.
''',
      
    keywords=['physics', 'electrostatics',
              'trap'],
    
    classifiers=[
    'Development Status :: 2 - Pre-Alpha',
    'Intended Audience :: Education',
    'Intended Audience :: Other Audience',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Operating System :: OS Independent',
    'Programming Language :: Python',
    'Programming Language :: Python :: 2.5',
    'Programming Language :: Python :: 2.6',
    'Programming Language :: Python :: 2.7',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Mathematics',
    'Topic :: Scientific/Engineering :: Physics',
    'Topic :: Utilities'
    ],
    
    # Files are defined in MANIFEST
    packages=['pytrap'],
        
    )  # End of setup definition

