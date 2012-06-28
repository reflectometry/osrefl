from osrefl.model.sample_prep import *
from osrefl.model.samples import *
from numpy import log, abs, min, max
from pylab import figure, show, subplot, imshow
from osrefl.loaders.andr_load import *
from osrefl.theory.scatter import *
from osrefl.theory import approximations, scatter
from osrefl.theory.approximations import *
from osrefl.viewers import view
import numpy as np
import matplotlib.pyplot as plt



# Define the Samples
alternating = AlternatingSample(shell_dim = [3000.0, 3000.0, 500.0], 
                                core_dim = [3000.0, 200.0, 500.0],
                                x_increment = 0.0,
                                y_increment = 500.0,                  
                                offset = [0.0, -1400.0, 0.0])


#alternating = AlternatingSample(shell_dim = [3000.0, 3000.0, 500.0], 
#                                core_dim = [200.0, 3000.0, 500.0],
#                                x_increment = 500.0,
#                                y_increment = 0.0,                  
#                                offset = [-1400, 0, 0.0])

triprism = TriPrismSample(shell_dim = [3000.0, 3000.0, 400.0], 
                          core_dim = [3000.0, 500.0, 300.0],
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
alternating.Create(core_SLD = 3.0e-6, core_Ms = 2.162e-6, 
                   shell_SLD = 1.0e-6, shell_Ms = 1.0e-6)
triprism.Create(core_SLD = 1e-6, core_Ms = 2.162e-6, 
                shell_SLD = 3.12e-6, shell_Ms = 1.0e-6)
cylinder.Create(core_SLD = 4.0e-6, core_Ms = 2.162e-6, 
                shell_SLD = 3.12e-6, shell_Ms = 1.0e-6)

# Grab the scene object 
altscene = alternating.getScene()
triprismscene = triprism.getScene()
cylinderscene = cylinder.getScene()

# Define Resolution for Slice Drawing
resolution = [ 100 , 100 , 100 ]
xyz = [ 3000 , 3000 , 600 ]

# Define and build the Geometry Unit for the different samples
print("Creating Geometry Units...")
altunit = GeomUnit(Dxyz = xyz, n = resolution, scene = altscene)
altunit = altunit.buildUnit()
#altunit.add_media()

triprismunit = GeomUnit(Dxyz = xyz, n = resolution, scene = triprismscene)
triprismunit = triprismunit.buildUnit()
#triprismunit.add_media()

cylinderunit = GeomUnit(Dxyz = xyz, n = resolution, scene = cylinderscene)
cylinderunit = cylinderunit.buildUnit()
#cylinderunit.add_media()

# Define the Q space
q_space = Q_space([ -0.1 , -0.1 , 0.002 ], [ 0.1 , 0.1 , .12 ], [ 300 , 300 , 300 ])

#define the lattice Structures of each unit
altlattice = Rectilinear([1,1,1],altunit)
triprismlattice = Rectilinear([1,1,1],triprismunit)
cylinderlattice = Rectilinear([1,1,1],cylinderunit)


# Define the Beam parameters 
beam = Beam(4.75, 0.02, None, 0.02, None)

# Calculate the Scattering and display results 
###############################################################################

print("Calculating Sample 1... (Alternating)")

sample1 = scatter.Calculator(None, beam, q_space, altunit)

raw_intensity1 = sample1.BA_FormFactor()

##############################################################################

print("Calculating Sample 2... (Saw Tooth)")

sample2 = scatter.Calculator(None, beam, q_space, triprismunit)

raw_intensity2 = sample2.BA_FormFactor()

###############################################################################

print("Calculating Sample 3... (Cylinder)")

sample3 = scatter.Calculator(None, beam, q_space, cylinderunit)

raw_intensity3 = sample3.BA_FormFactor()

###############################################################################

# View Angular Results

print("Viewing Sample 1... (Alternating)")
sample1.toAngular(0.25, raw_intensity1)
sample1.viewAngular()

print("Viewing Sample 1... (Saw Tooth)")
sample2.toAngular(0.25, raw_intensity2)
sample2.viewAngular()

print("Viewing Sample 1... (Cylinder)")
sample3.toAngular(0.25, raw_intensity3)
sample3.viewAngular()

###############################################################################



















