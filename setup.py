#!/usr/bin/env python

from distutils.core import setup

setup(
    name='ARC',
    version='1.1.0',
    author='Institute for Bioinformatics and Evolutionary Studies',
    url='http://www.ibest.uidaho.edu',
    packages=['ARC', 'ARC.runners'],
    scripts=['bin/ARC']
)
