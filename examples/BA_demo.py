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


# Define the Samples
alternating = AlternatingSample(shell_dim = [3000.0, 3000.0, 600.0], 
                                core_dim = [3000.0, 200.0, 600.0],
                                x_increment = 0.0,
                                y_increment = 500.0,                          
                                offset = [0.0, -1400.0, 0.0])

triprism = TriPrismSample(shell_dim = [3000.0, 3000.0, 600.0], 
                                core_dim = [3000.0, 900.0, 300.0], # y seems to be off
                                x_increment = 0.0,
                                y_increment = 500.0,                          
                                offset = [0.0, -1250.0, 0.0])

cylinder = CylinderSample(shell_dim = [3000.0, 3000.0, 500.0], 
                          core_dim = [200.0, 200.0, 700.0],        # z seems to be off
                          x_increment = 1200.0,
                          y_increment = 600.0,
                          offset = [-1500.0, -1500.0, -300.0],
                          offset2 = [600.0, 300.0, 0])

# Build the Samples
alternating.Create()
triprism.Create()
cylinder.Create()

# Grab the scene object 
triprismscene = triprism.getScene()
altscene = alternating.getScene()
cylinderscene = cylinder.getScene()

# Define and build the Geometry Unit for the Cylinder sample
triprismunit = GeomUnit(Dxyz = [3000.0,3000.0,600.0], n = [150,150,150], scene = triprismscene)
triprismunit = triprismunit.buildUnit()
triprismunit.add_media()

altunit = GeomUnit(Dxyz = [3000.0,3000.0,600.0], n = [150,150,150], scene = altscene)
altunit = altunit.buildUnit()
altunit.add_media()

cylinderunit = GeomUnit(Dxyz = [3000.0,3000.0,600.0], n = [150,150,150], scene = cylinderscene)
cylinderunit = cylinderunit.buildUnit()
cylinderunit.add_media()

# View composition of the Geometry Unit
triprismunit.viewSlice()
altunit.viewSlice()
cylinderunit.viewSlice()

'''

# Define the Q space and Lattice structure
q_space = Q_space([-.001,-0.002,0.0002],[.001,.002,0.3],[100,25,50])
lattice = Rectilinear([20,20,1],altunit)

# Define the Beam parameters 
beam = Beam(5.0,None,None,0.05,None)

# Calculate the Scattering and display results 
sample1 = scatter.Calculator(lattice,beam,q_space,altunit)
sample1.SMBAfft(refract = False)
sample1.resolution_correction()
sample1.viewCorUncor()

'''








