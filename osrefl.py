#!/usr/bin/env python

# Copyright (C) 2008-2011 University of Maryland
# All rights reserved.
# See LICENSE.txt for details.
# Author: Christopher Metting

"""
This script starts the OsRefl Off-Specular Reflectometry application.
"""
import os
import sys
import matplotlib
matplotlib.use('WXAgg')

# Normally the inversion package will be installed, but if it is not installed,
# augment sys.path to include the parent directory of the package.  This
# assures that the module search path will include the package namespace and
# allows the application to be run directly from the source tree, even if the
# package has not been installed.
try:
    import osrefl
except:
    this_dir_path = os.path.dirname(os.path.abspath(__file__))
    if os.path.basename(this_dir_path) == 'osrefl':
        sys.path.insert(1, (os.path.dirname(this_dir_path)))
    else:
        print """\
        *** To run this script, either install the inversion package or
        *** place this module in the top-level directory of the package."""

if __name__ == "__main__":
    if len(sys.argv) == 1:
        file = os.path.join('examples', 'AuFit.py')
    else:
        file = sys.argv[1]

    try:
        f = open(file)
        f.close()
    except:
        print "*** Script not found:", file
        sys.exit()

    print "Executing script", file, "...\n"
    execfile(file)
