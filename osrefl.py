#!/usr/bin/env python

# Copyright (C) 2008-2011 University of Maryland
# All rights reserved.
# See LICENSE.txt for details.
# Author: Christopher Metting

"""
This script starts the DiRefl Direct Inversion Reflectometry application.
"""

import os
import sys

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
    #import osrefl.examples.AuFit
    #osrefl.someplace.AuFit.main()
    print os.path.join('examples', 'AuFit.py')
    execfile(os.path.join('examples', 'AuFit.py'))
