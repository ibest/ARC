#!/usr/bin/env python

from distutils.core import setup

setup(
    name='ARC',
    version='2.0.0',
    author='Institute for Bioinformatics and Evolutionary Studies',
    url='http://www.ibest.uidaho.edu',
    packages=['ARC', 'ARC.runners'],
    scripts=['bin/ARC']
)
