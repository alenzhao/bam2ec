#!/usr/bin/env python
# -*- coding: utf-8 -*-
from glob import glob
import os


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read().replace('.. :changelog:', '')

requirements = [
    # TODO: put package requirements here
]

test_requirements = [
    # TODO: put package test requirements here
]

on_rtd = os.environ.get('READTHEDOCS', None)
if not on_rtd:
    requirements.append('pysam>=0.8.1')
    requirements.append('numpy>1.8')
    requirements.append("emase>=0.9.8")


setup(
    name='bam2ec',
    version='0.1.0',
    description="Convert BAM and SAM files to binary format with Equivalence Classes",
    long_description=readme + '\n\n' + history,
    author="Matthew Vincent",
    author_email='mvincent@jax.org',
    url='https://github.com/mvincent/bam2ec',
    packages=[
        'bam2ec',
    ],
    package_dir={'bam2ec':
                 'bam2ec'},
    include_package_data=True,
    install_requires=requirements,
    license="ISCL",
    zip_safe=False,
    scripts=glob("bin/*"),
    keywords='bam2ec',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: ISC License (ISCL)',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
