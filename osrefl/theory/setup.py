#!/usr/bin/env python

import os
import sys

from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

def configuration(parent_package='', top_path=None):
    config = Configuration('theory', parent_package, top_path)
    config.set_options(quiet=True) # silence debug/informational messages

    # Just treat lib as a data directory for now.
    config.add_data_dir('lib')

    return config


if __name__ == '__main__':
    if len(sys.argv) == 1: sys.argv.append('install')
    setup(**configuration(top_path='').todict())
