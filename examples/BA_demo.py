# Copyright (C) 2008 University of Maryland
# All rights reserved.
# See LICENSE.txt for details.
# Author: Christopher Metting

#Starting Date:6/12/2009
from osrefl.model.sample_prep import *
from numpy import log, abs, min, max
from pylab import figure, show, subplot, imshow
from osrefl.loaders.andr_load import *
from osrefl.theory.scatter import *
from osrefl.theory import approximations, scatter
from osrefl.viewers import view

# Variables used to define the parameters of multiple Ellipse objects without
#redefining each individually.
shell_dim =[3000.0, 3000.0, 400.0]

core_dim =[3000.0, 1000.0, 300.0]

MaterialB = Parallelapiped(SLD = 3.12e-6, dim = shell_dim, Ms = 1.0e-6)

MaterialA = []

offsetx = 0.0
offsety = -1250.0 #- (core_dim[0] / 2 + core_dim[1] / 2)
offsetz = 0.0

offset_iterator = [offsetx, offsety, offsetz]

for i in range(6):  
    MaterialA.append(TriangularPrism(SLD = 2.0e-6, dim = core_dim, Ms = 1.5e-6))
    MaterialA[i].is_core_of(MaterialB, offset_iterator)
    offset_iterator[1] += 500.0
    



scenelist = []
scenelist.append(MaterialB)
#scenelist.append(MaterialA)
for i in range(len(MaterialA)):
    scenelist.append(MaterialA[i])

test = Scene(scenelist)

'''
now we can create the unit_cell. This class contains all of the required
information about a single unit cell.
'''

kathryns_non_mag_unit = GeomUnit(Dxyz = [6500.0,4500.0,None],
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

sample.SMBAfft(refract = False)

sample.resolution_correction()

sample.viewCorUncor()

'''
********* INSERT SAVE FOR sample HERE**********
'''








