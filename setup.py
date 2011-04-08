#!/usr/bin/env python

import os
import sys

from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

def configuration(parent_package='', top_path=None):
    config = Configuration('osrefl', parent_package, top_path)
    config.set_options(quiet=True) # silence debug/informational messages

    # Add subpackages (top level name spaces) and data directories.
    # Note that subpackages may have their own setup.py to drill down further.
    config.add_subpackage('osrefl')
    config.add_data_dir('doc')
    config.add_data_dir('examples')
    config.add_data_dir('tests')

    config.add_data_files('LICENSE.txt')
    config.add_data_files('README.txt')

    return config


if __name__ == '__main__':
    if len(sys.argv) == 1: sys.argv.append('install')
    setup(**configuration(top_path='').todict())
