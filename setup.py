#!/usr/bin/env python
import sys
import os

# Find the version number of the application
for line in open("osrefl/__init__.py").readlines():
    if line.startswith('__version__'):
        exec line

#from distutils.core import Extension
from setuptools import setup, find_packages, Extension
#import fix_setuptools_chmod

packages = find_packages()

if len(sys.argv) == 1:
    sys.argv.append('install')

# reflmodule extension
def reflmodule_config():
    srcpath = os.path.join('osrefl','loaders','reduction','lib')
    sources = [os.path.join(srcpath,f)
               for f in "reduction.cc","str2imat.c"]
    module = Extension('osrefl.loaders.reduction._reduction', sources=sources)
    return module

#TODO: write a proper dependency checker for packages which cannot be
# installed by easy_install
#dependency_check('numpy>=1.0', 'scipy>=0.6', 'matplotlib>=1.0', 'wx>=2.8.9')

dist = setup(
        name = 'osrefl',
        version = __version__,
        author='Christopher Metting',
        author_email='mettingc@umd.edu',
        url='http://www.reflectometry.org/danse/osrefl.html',
        description='3-D reflectometry modelling',
        long_description=open('README.txt').read(),
        classifiers=[
            'Development Status :: 4 - Beta',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: Public Domain',
            'Operating System :: OS Independent',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering :: Chemistry',
            'Topic :: Scientific/Engineering :: Physics',
            ],
        packages = packages,
        scripts = ['osrefl.py'],
        ext_modules = [reflmodule_config()],
        )

# End of file
