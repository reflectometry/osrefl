# Copyright (C) 2008 University of Maryland
# All rights reserved.
# See LICENSE.txt for details.
# Author: Christopher Metting

#Starting Date:6/12/2009
from osrefl.model.sample_prep import *
from osrefl.model.samples import *
from numpy import log, abs, min, max
from pylab import figure, show, subplot, imshow
from osrefl.loaders.andr_load import *
from osrefl.theory.scatter import *
from osrefl.theory import approximations, scatter
from osrefl.viewers import view


# Define the Samples
alternating = AlternatingSample(shell_dim = [3000.0, 3000.0, 500.0], 
                                core_dim = [3000.0, 200.0, 500.0],
                                x_increment = 0.0,
                                y_increment = 500.0,                          
                                offset = [0.0, -1400.0, 0.0])

triprism = TriPrismSample(shell_dim = [3000.0, 3000.0, 400.0], 
                          core_dim = [3000.0, 500.0, 300.0], # z seems to be off
                          x_increment = 0.0,
                          y_increment = 500.0,                          
                          offset = [0.0, -1250.0, 0.0])

cylinder = CylinderSample(shell_dim = [3000.0, 3000.0, 500.0], 
                          core_dim = [200.0, 200.0, 300.0], 
                          x_increment = 1200.0,
                          y_increment = 600.0,
                          offset = [-1500.0, -1500.0, -100.0],
                          offset2 = [600.0, 300.0, 0.0])

# Build the Samples with different materials
print("\nBuilding Samples...")
alternating.Create(core_SLD = 2.0e-6, core_Ms = 2.162e-6, 
                   shell_SLD = 3.12e-6, shell_Ms = 1.0e-6)
triprism.Create(core_SLD = 2.0e-6, core_Ms = 2.162e-6, 
                shell_SLD = 3.12e-6, shell_Ms = 1.0e-6)
cylinder.Create(core_SLD = 2.0e-6, core_Ms = 2.162e-6, 
                shell_SLD = 3.12e-6, shell_Ms = 1.0e-6)
print("...")
# Grab the scene object 
altscene = alternating.getScene()
triprismscene = triprism.getScene()
cylinderscene = cylinder.getScene()

# Define Resolution for Slice Drawing
resolution = [200,200,200]

xyz = [9000.0, 9000.0, 600.0]

# Define and build the Geometry Unit for the different samples
print("Creating Geometry Units...")
altunit = GeomUnit(Dxyz = xyz, n = resolution, scene = altscene)
altunit = altunit.buildUnit()
altunit.add_media()
print("...")
triprismunit = GeomUnit(Dxyz = xyz, n = resolution, scene = triprismscene)
triprismunit = triprismunit.buildUnit()
triprismunit.add_media()
print("...")
cylinderunit = GeomUnit(Dxyz = xyz, n = resolution, scene = cylinderscene)
cylinderunit = cylinderunit.buildUnit()
cylinderunit.add_media()
print("...")
# Define the Q space
q_space = Q_space([-.001,-0.002,0.0002],[.001,.002,0.3],[100,25,50])

#define the lattice Structures of each unit
print("Creating Lattice Structures...")
altlattice = Rectilinear([20,20,1],altunit)
triprismlattice = Rectilinear([20,20,1],triprismunit)
cylinderlattice = Rectilinear([20,20,1],cylinderunit)

# Define the Beam parameters 
beam = Beam(5.0,None,None,0.05,None)

# Calculate the Scattering and display results 
print("Calculating Sample 1...")
sample1 = scatter.Calculator(altlattice,beam,q_space,altunit)
sample1.SMBAfft(refract = False)
sample1.resolution_correction()

print("\nCalculating Sample 2...")
sample2 = scatter.Calculator(triprismlattice,beam,q_space,triprismunit)
sample2.BA()
sample2.resolution_correction()

print("\nCalculating Sample 3...")
sample3 = scatter.Calculator(cylinderlattice,beam,q_space,cylinderunit)
sample3.DWBA(refract = False)
sample3.resolution_correction()

# View Results
print("\nViewing Results...")
altunit.viewSlice()
sample1.viewCorUncor()

triprismunit.viewSlice()
sample2.viewCorUncor()

cylinderunit.viewSlice()
sample3.viewCorUncor()

