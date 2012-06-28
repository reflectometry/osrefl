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
                   shell_SLD = 2.12e-6, shell_Ms = 1.0e-6)
triprism.Create(core_SLD = 2.0e-6, core_Ms = 2.162e-6, 
                shell_SLD = 3.12e-6, shell_Ms = 1.0e-6)
cylinder.Create(core_SLD = 4.0e-6, core_Ms = 2.162e-6, 
                shell_SLD = 3.12e-6, shell_Ms = 1.0e-6)

# Grab the scene object 
altscene = alternating.getScene()
#triprismscene = triprism.getScene()
#cylinderscene = cylinder.getScene()

# Define Resolution for Slice Drawing
resolution = [100,100,50]

#xyz = [12000,12000,1200]
xyz = [3000,3000,600]

# Define and build the Geometry Unit for the different samples
print("Creating Geometry Units...")
altunit = GeomUnit(Dxyz = xyz, n = resolution, scene = altscene)
altunit = altunit.buildUnit()
#altunit.add_media()

#triprismunit = GeomUnit(Dxyz = xyz, n = resolution, scene = triprismscene)
#triprismunit = triprismunit.buildUnit()
#triprismunit.add_media()

#cylinderunit = GeomUnit(Dxyz = xyz, n = resolution, scene = cylinderscene)
#cylinderunit = cylinderunit.buildUnit()
#cylinderunit.add_media()

# Define the Q space
#150 150 150
q_space = Q_space([-0.01,-0.1,0.002],[0.01,0.1,.12],[300, 300, 300])
#q_space = Q_space([-10.0, -1.0, 0.0002],[10.0, 1.0, 14.0],[100,100,100])
#q_space = Q_space([-.10,-0.15,0.02],[.10,0.15,.14],[100,100,20])

#define the lattice Structures of each unit
altlattice = Rectilinear([1,1,1],altunit)
#triprismlattice = Rectilinear([20,20,1],triprismunit)
#cylinderlattice = Rectilinear([1,1,1],cylinderunit)


# Define the Beam parameters 
#beam = Beam(5.0,None,None,0.05,None)
beam = Beam(4.75, 0.02, None, 0.02, None)

# Calculate the Scattering and display results 
###############################################################################
'''
print("Calculating Sample 1...")

sample1 = scatter.Calculator(cylinderlattice, beam, q_space, cylinderunit)

raw_intensity = BA_FT(cylinderunit.unit, cylinderunit.step, q_space)
raw_intensity = abs(complete_formula(raw_intensity, cylinderunit.step, q_space))**2

qx_array = q_space.q_list[0].reshape(q_space.points[0],1,1)
qy_array = q_space.q_list[1].reshape(1,q_space.points[1], 1)
qz_array = q_space.q_list[2].reshape(1,1,q_space.points[2])

raw_intensity *= (4.0 * pi / (qz_array * qx_array * qy_array))**2
#raw_intensity *= (4.0 * pi / (qz_array))**2
#raw_intensity *= (4.0 * pi / (qz_array * altunit.Dxyz[0] * altunit.Dxyz[1]))**2 
 
sample1.results = sum(raw_intensity,1).astype('float64')
'''

##############################################################################

print("Calculating Sample 2...")

sample2 = scatter.Calculator(altlattice, beam, q_space, altunit)

raw_intensity = sample2.BA_FormFactor()

sample2.results = sum(raw_intensity, axis = 1)
sample2.viewUncor()

#sample2.DWBA()
#sample2.viewUncor()

sample2.toAngular(0.5, raw_intensity)
sample2.viewAngular()

###############################################################################

###############################################################################
'''
print("Calculating Sample 3...")

triprismunit.viewSlice()

sample3 = scatter.Calculator(None, beam, q_space, triprismunit)
raw_intensity = BA_FT(triprismunit.unit, triprismunit.step, q_space)
raw_intensity = abs(complete_formula(raw_intensity, triprismunit.step, q_space))**2

qx_array = q_space.q_list[0].reshape(q_space.points[0],1,1)
qy_array = q_space.q_list[1].reshape(1,q_space.points[1], 1)
qz_array = q_space.q_list[2].reshape(1,1,q_space.points[2])

raw_intensity *= (4.0 * pi / (qz_array * qx_array * qy_array))**2
    
sample3.results = sum(raw_intensity,axis=1).astype('float64')

# View Results
print("\nViewing Results for Sample 1 (cylinder)...")

sample1.viewUncor()

print("\nViewing Results for Sample 2 (alternating)...")

sample2.viewUncor()

print("\nViewing Results for Sample 3 (triprism)...")

sample3.viewUncor()

'''
