# Copyright (C) 2008 University of Maryland
# All rights reserved.
# See LICENSE.txt for details.
# Author: Christopher Metting

#Starting Date:9/21/2009
import os

import numpy
from numpy import log

from pylab import figure

from osrefl.model.sample_prep import *
from osrefl.loaders.andr_load import *
from osrefl.theory.scatter import *
from osrefl.theory import approximations, scatter
from osrefl.viewers import view

# The following line does not work on Windows - jak
#os.nice(10)

Au_measurments = Data()
#Au = Ellipse(SLD = 4.506842e-6,dim=[3.75e4,3.75e4,600.0])
Au = [None]*1
test_data = [None]*1
#shape test
stdDim=[3.9e4,3.75e4,550.0]
Au = RoundedParPip(SLD = 4.506842e-6,dim=[3.75e4,3.75e4,630.0], curve = .56)
Cr = Layer(SLD = 3.01e-6,thickness_value = 48.0)

Au.on_top_of(Cr)

scene = Scene([Au,Cr])
GeoUnit = GeomUnit(Dxyz = [10.0e4,10.0e4,700.0], n = [100,100,300],scene = scene, inc_sub = [0.0,2.7e-6])
unit = GeoUnit.buildUnit()
#unit.viewSlice()
print size(unit.unit[unit.unit==4.506842e-6])
unit.add_media()

lattice = Rectilinear([20.0,20.0,1.0],unit)

beam = Beam(5.0,.02,None,0.02,None)
scale  = 1.7e8
#q_space = Q_space([-.0002,-0.002,0.00002],[.0002,.002,0.03],[200,50,200])
q_space = Au_measurments.space


test_data = Calculator(lattice,beam,q_space,unit)
test_data.SMBAfft()
#test_data.results[test_data.results < 1.0e-15] = 1.0e-15
test_data.resolution_correction()
test_data.corrected_results *=scale
#test_data.results[Au_measurments.data==0.0]=0.0


#Masking
artre = min(test_data.corrected_results[nonzero(test_data.corrected_results)])
print artre

test_data.corrected_results[test_data.corrected_results == 0.0] = artre
test_data.corrected_results[test_data.corrected_results == NaN] = artre
test_data.corrected_results[Au_measurments.data==0.0]=0.0

'''
from numpy import save
save('mod.npy',test_data.corrected_results)
save('data.npy',Au_measurments.data)
'''
#extraData = [test_data[1].corrected_results,test_data[2].corrected_results,test_data[3].corrected_results,test_data[4].corrected_results]
#test_data[0].fitCompare(Au_measurments,extraData,titles = ['data','curve = 0','curve = 25','curve = 56','curve = 75','curve = 100'])
test_data.fitCompare(Au_measurments,titles = ['data','Au-56% Curve'])
#test_data.scale(Au_measurments)
