#! /usr/bin/env python

import os
import sys
import subprocess

from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand

VERSION = '0.1.0'

#stolen from Guidebook
PKG_NAME = 'dgbio'
def include(*dirs):
    'Generate explicit mapping of asset -> destination for setuptools'
    results = []
    for src_dir in dirs:
        for root, _, files in os.walk(src_dir):
            results.append((os.path.join(PKG_NAME, root),
                            map(lambda f: os.path.join(root, f), files)))
    return results


class Tox(TestCommand):
    user_options = [('tox-args=', 'a', "Arguments to pass to tox")]
    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.tox_args = None
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True
    def run_tests(self):
        #import here, cause outside the eggs aren't loaded
        import tox
        import shlex
        errno = tox.cmdline(args=shlex.split(self.tox_args))
        sys.exit(errno)

setup(
    name='big-search',
    version=VERSION,
    url='https://github.com/DeskGen/big-search',
    author='YOUR NAME',
    author_email='YOUR EMAIL',
    maintainer='',
    maintainer_email='',
    description='Case Study in Python searching datasets larger than RAM.',
    packages=find_packages(),
    data_files=include('bin'),
    install_requires=[
        'pysam',
        'numpy',
    ],
    tests_require=[
        'virtualenv',
        'tox',
        'pytest',
        'mock',
    ],
    cmdclass={
        'test': Tox,
    },
)
