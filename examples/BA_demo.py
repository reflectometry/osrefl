# Copyright (C) 2008 University of Maryland
# All rights reserved.
# See LICENSE.txt for details.
# Author: Christopher Metting

#Starting Date:6/12/2009
from osrefl.model.sample_prep import *
from numpy import log,abs, min, max
from pylab import figure, show, subplot, imshow
from osrefl.loaders.andr_load import *
from osrefl.theory.scatter import *
from osrefl.theory import approximations, scatter
from osrefl.viewers import view

# Variables used to define the parameters of multiple Ellipse objects without
#redefining each individually.
shell_dim =[5500.0,3500.0,100.0]
core_dim = [2000.0,2000.0,100.0]

'''
each layer is defined by a set of parameters. These are the same parameters that
will be used for fitting. The layer need at minimum a Scattering Length density
and some dimension. 
'''
NiFe = Ellipse(SLD = 9.12e-6,dim = shell_dim, Ms = 2.162e-6)
NiFe_mag = Ellipse(SLD = 2.162e-6,dim = core_dim)
Cu = Ellipse(SLD = 6.54e-6,dim = shell_dim)
Co = Ellipse(SLD = 2.26e-6,dim = shell_dim)
Co_mag = Ellipse(SLD = 4.12e-6,dim = core_dim)
Au = Ellipse(SLD = 4.5e-6,dim = shell_dim)

#IrMn = Layer(SLD = -0.06585e-6,thickness_value = 40.0)
IrMn = Layer(SLD = 6.585e-6,thickness_value = 40.0)

'''
By using some built-in functions, we can specify the locations of shapes
relative to one another
'''
#NiFe.on_top_of(IrMn)
NiFe_mag.is_core_of(NiFe)
Cu.on_top_of(NiFe)
Co.on_top_of(Cu)
Co_mag.is_core_of(Co)
Au.on_top_of(Co)

'''
once we are satisfied with the shapes that we have, we can create a scene.
Although it has not been implemented yet, this will effectively tie the shapes
to each other and create one main shape that will be used for the fitting
'''
#mag_ellips_scene = Scene([IrMn,NiFe,NiFe_mag,Cu,Co,Co_mag,Au])
mag_ellips_scene = Scene([IrMn,NiFe,NiFe_mag,Cu])
#non_mag_ellips_scene = Scene([IrMn,NiFe,Cu,Co,Au])
non_mag_ellips_scene = Scene([IrMn,NiFe,Cu])
core_ellips_scene = Scene([NiFe_mag,Co_mag])
test = Scene([NiFe])
'''
now we can create the unit_cell. This class contains all of the required
information about a single unit cell.
'''

kathryns_non_mag_unit = GeomUnit(Dxyz = [9000.0,4500.0,None],
                                 n = [20,20,20], scene = test)

kathryns_non_mag_unit = kathryns_non_mag_unit.buildUnit()
kathryns_non_mag_unit.add_media()

q_space = Q_space([-.001,-0.002,0.0002],[.001,.002,0.3],[100,25,50])
lattice = Rectilinear([20,20,1],kathryns_non_mag_unit)

beam = Beam(5.0,None,None,0.05,None)

'''
This solves the structure
BA for a single, non-magnetic system
'''
sample = scatter.Calculator(lattice,beam,q_space,kathryns_non_mag_unit)
sample_two = scatter.Calculator(lattice,beam,q_space,kathryns_non_mag_unit)

sample.DWBA(refract = False)

sample.resolution_correction()

sample.viewCorUncor()

'''
********* INSERT SAVE FOR sample HERE**********
'''








