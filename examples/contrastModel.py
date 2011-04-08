# Copyright (C) 2008 University of Maryland
# All rights reserved.
# See LICENSE.txt for details.
# Author: Christopher Metting

#Starting Date:9/21/2009
import numpy
from sample_prep import *
from andrload import *
from scatter import *
import approximations,view, scatter
from numpy import log
from pylab import figure


Au_measurments = Data()
#Au = Ellipse(SLD = 4.506842e-6,dim=[3.75e4,3.75e4,600.0])
Au = [None]*1
test_data = [None]*1
#shape test
stdDim=[3.9e4,3.75e4,550.0]
fillMediaSLD = [4.506842e-6,5.506842e-6,6.506842e-6]
feHight = 630.0
liquid = [None]*3
models = [None]*3
print Au_measurments.type
liquid[0] = Layer(SLD = fillMediaSLD[0],thickness_value = feHight)
liquid[1] = Layer(SLD = fillMediaSLD[1],thickness_value = feHight)
liquid[2] = Layer(SLD = fillMediaSLD[2],thickness_value = feHight)

for i in range(size(liquid)):
    Au = RoundedParPip(SLD = 4.506842e-6,dim=[3.75e4,3.75e4,feHight], curve = .56)
    Cr = Layer(SLD = 3.01e-6,thickness_value = 48.0)
    
    liquid[i].on_top_of(Cr)
    Au.on_top_of(Cr)
    
    scene = Scene([liquid[i],Cr])
    GeoUnit = GeomUnit(Dxyz = [10.0e4,10.0e4,700.0], n = [100,100,300],scene = scene, inc_sub = [liquid[i].SLD,2.7e-6])
    unit = GeoUnit.buildUnit()
    unit.add_media()
    #unit.viewSlice()
    
    
    lattice = Rectilinear([20.0,20.0,1.0],unit)
    
    beam = Beam(5.0,.02,None,0.02,None)
    scale  = 1.7e4
    q_space = Q_space([-.0002,-0.002,0.00002],[.0002,.002,0.03],[200,50,200])
    #q_space = Au_measurments.space
    
    
    test_data = Calculator(lattice,beam,q_space,unit)
    test_data.BA()
    #test_data.results[test_data.results < 1.0e-15] = 1.0e-15 
    test_data.resolution_correction()
    test_data.corrected_results *=scale
    #test_data.results[Au_measurments.data==0.0]=0.0
    
        
    #Masking
    artre = min(test_data.corrected_results[nonzero(test_data.corrected_results)])
    
    
    test_data.corrected_results[test_data.corrected_results == 0.0] = artre
    test_data.corrected_results[test_data.corrected_results == NaN] = artre
    #test_data.corrected_results[Au_measurments.data==0.0]=0.0
    
    models[i] = test_data
'''
from numpy import save
save('mod.npy',test_data.corrected_results)
save('data.npy',Au_measurments.data)
'''
extraData = [models[2].corrected_results]
#extraData = [test_data[1].corrected_results,test_data[2].corrected_results,test_data[3].corrected_results,test_data[4].corrected_results]
#test_data[0].fitCompare(Au_measurments,extraData,titles = ['data','curve = 0','curve = 25','curve = 56','curve = 75','curve = 100'])
models[0].fitCompare(models[1],extraData,titles = ['data',str(fillMediaSLD[0]),str(fillMediaSLD[1]),str(fillMediaSLD[2])])
#test_data.scale(Au_measurments)