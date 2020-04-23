#!/usr/bin/env python

from setuptools import setup
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name='angcorrwat',
    version='0.0.1',
    description='Calculate angular distributions and correlations',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/op3/angcorrwat', 
    author='O. Papst',
    author_email='opapst@ikp.tu-darmstadt.de',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    packages=['angcorrwat'],
    install_requires=['sympy'],
    setup_requires=['pytest-runner'],
    tests_require=['pytest', 'pytest-cov'],
)
