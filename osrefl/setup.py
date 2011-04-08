#!/usr/bin/env python

import os
import sys

from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

def configuration(parent_package='', top_path=None):
    config = Configuration('osrefl', parent_package, top_path)
    config.set_options(quiet=True) # silence debug/informational messages

    config.add_subpackage('loaders')
    config.add_subpackage('model')
    config.add_subpackage('theory')
    config.add_subpackage('viewers')

    return config


if __name__ == '__main__':
    if len(sys.argv) == 1: sys.argv.append('install')
    setup(**configuration(top_path='').todict())
