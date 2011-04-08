# Copyright (C) 2008 University of Maryland
# All rights reserved.
# See LICENSE.txt for details.
# Author: Christopher Metting

#Starting Date:6/12/2009

from sample_prep import *
#from . sample_prep import *
import approximations,view, scatter
from numpy import log,abs, min, max, amax,where, flipud
from pylab import figure, show, subplot, imshow,colorbar,plot,title,xlabel,ylabel


#Au = Parallelapiped(SLD = 4.506842e-5,dim=[5.0e4,5.0e4,550.0])
Au = Ellipse(SLD = 4.506842e-6,dim=[3.75e4,3.75e4,700.0])
Cr = Layer(SLD = 3.01e-6,thickness_value = 20.0)

Au.on_top_of(Cr)
scene = Scene([Au,Cr])

GeoUnit = GeomUnit(Dxyz = [10.0e4,10.0e4, 2000.0], n = [50,51,52],scene = scene)#, inc_sub = [4.506842e-6,2.7e-6])
unit = GeoUnit.buildUnit()
#unit.add_media()

#unit.viewSlice()
q_space = Q_space([-.0001,-0.001,0.00002],[.0001,0.001,0.04],[50,20,50])
lattice = Rectilinear([20,20,1],unit)

beam = Beam(5.0,None,None,0.05,None)
sample = scatter.Calculator(lattice,beam,q_space,unit)
sample2 = scatter.Calculator(lattice,beam,q_space,unit)
sample3 = scatter.Calculator(lattice,beam,q_space,unit)

import time


t0 = time.time()
sample.BA()
print "time of BA",time.time()-t0
sample.resolution_correction()
sample.viewCorUncor()
'''
sample.view_uncorrected()
sample.resolution_correction()
sample.view_corrected()
show()
'''
'''
t0 = time.time()
sample2.longBA()
print "time of float64",time.time()-t0
'''
'''
t0 = time.time()
sample3.BA()
print "time of BA",time.time()-t0
'''
'''
sub = abs(sample.results - sample2.results)/(abs(sample.results + 1.0e-20)).real 
print max(sub)

title('Relative Difference btwn Python64 and Cuda64')
ylabel('qz axis')
xlabel('qx axis')
extent = [q_space.minimums[0],q_space.maximums[0],q_space.minimums[2],q_space.maximums[2]]
imshow(flipud(log(sum(sub.real,axis = 1).T)),extent = extent, aspect = 'auto')

colorbar()
show()

sub = abs(sample.results - sample2.results)/((abs(sample.results + sample2.results)/2.0)).real 
print max(sub)
'''
'''
sliceone = sum(sample.results,axis = 1)
slicetwo = sum(sample2.results,axis = 1)
'''

