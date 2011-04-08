# Copyright (C) 2008 University of Maryland
# All rights reserved.
# See LICENSE.txt for details.
# Author: Christopher Metting

#Starting Date:6/5/2009
from numpy import size,array,shape,indices, searchsorted, linspace
from numpy import log, log10, abs, min, max, nonzero,isnan
from zoom_colorbar import *
import sys,copy

from numpy import roll

def feature_plot(unit,size,delta_xyz):
    '''
    A three dimensional plot the the voxilized feature. Does what spy_on but
    does a three dimensional plot rather then z slices.
    
    feature (type: Lego_Model): contains all of the information required to
    calculate the scattering off the feature
    '''
    from enthought.mayavi import mlab
    
    dumbie_unit = roll(unit,1,axis = 0)
    dumbie_unit = roll(dumbie_unit,1,axis = 1)
    dumbie_unit[:,0,:] = 0.0
    dumbie_unit[:,-1,:] = 0.0
    dumbie_unit[0,:,:] = 0.0
    dumbie_unit[-1,:,:] = 0.0

    xyz = indices((unit.shape), 'd')
    

    xyz[0] *= delta_xyz[0]
    xyz[1] *= delta_xyz[1] 
    xyz[2] *= delta_xyz[2]
    

    feature_size = shape(unit)
    mlab.figure(0)
    s = mlab.contour3d(xyz[0],xyz[1],xyz[2],dumbie_unit,opacity=.07,contours = 20)
    mlab.figure(1)
    t = mlab.contour3d(xyz[0],xyz[1],xyz[2]*10,dumbie_unit,opacity=.07,contours = 20)
    mlab.figure(2)
    u = mlab.contour3d(dumbie_unit,opacity=.05,contours = 20)
    mlab.show()
    return


def intensity_plot(intensity,mins,maxs, header = None, bar = True, 
                   vmin = None, vmax = None):
    '''
    creates a two three dimensional plot of the qz vs qx and the intensity
    of the scattering.
    
    This plotter can be used for both resolution corrected and uncorrected
    intensity plots. The intensity has the lowest non-zero value added to it
    to eliminated the limitations that exists when taking the log of zero
    '''
    from pylab import imshow,colorbar,show, title, xlabel, ylabel
    
    
    plotxmin = mins[0]
    plotxmax = maxs[0] 
    plotzmin = mins[-1]
    plotzmax = maxs[-1]
    print min(intensity)
    if vmax == None:
        vmax = max(log10(intensity))
    print max(log(intensity))
    if vmin == None:
        vmin =  max(log10(intensity)) - 15.0
        
    intensity[isnan(intensity)] = 0.0

    if size(abs(intensity[nonzero(intensity.real)])) == 0:
        lower_lim = 0.0
    else:
        lower_lim = min(abs(intensity[nonzero(intensity.real)]))

    plot_extent = (plotxmin,plotxmax,plotzmin,plotzmax)

    graph = imshow(log10(abs(intensity.T+lower_lim)),aspect='auto',
                   interpolation='nearest',extent=plot_extent,origin='lower', 
                   vmin = vmin, vmax = vmax)
    zoom_colorbar(graph)
       
    title(str(header))
    xlabel('qx(A^-1)')
    ylabel('qz(A^-1)')
    
    return graph


def linear_plot(intensity,mins,maxs, header = None, bar = True,
                 vmin = None, vmax = None):
    '''
    creates a two three dimensional plot of the qz vs qx and the intensity
    of the scattering.
    
    This plotter can be used for both resolution corrected and uncorrected
    intensity plots. The intensity has the lowest non-zero value added to it
    to eliminated the limitations that excists when taking the log of zero
    '''
    from pylab import imshow,colorbar,show, title, xlabel, ylabel
    
    plotxmin = mins[0]
    plotxmax = maxs[0] 
    plotzmin = mins[-1]
    plotzmax = maxs[-1]

    lower_lim = min(intensity[nonzero(intensity.real)])

    plot_extent = (plotxmin,plotxmax,plotzmin,plotzmax)

    graph = imshow((abs(intensity.T+lower_lim)),aspect='auto',
                   interpolation='nearest',
                   extent=plot_extent,origin='lower')
    colorbar()
       
    title(str(header))
    xlabel('qz(A^-1)')
    ylabel('qx(A^-1)')
    
    return graph

def qz_slice(intensity,mins,maxs,q_slice = 0.0,second_intensity = None):
    
    '''
    This takes a qz slice from the uncorrected intensity plot and, if the
    resolution corrected plot exists, will also show the qz slice for that data
    '''
    from pylab import plot, legend, title, xlabel, ylabel
    
    print shape(intensity)
    qz_array = linspace(mins[2],maxs[2],shape(intensity)[1])
    z_position = searchsorted(qz_array,q_slice)

    graph = plot(log(intensity[z_position,:]),
                 xdata = qz_array,label = 'Uncorrected')

    if (second_intensity != None):
        plot(log(second_intensity[z_position,:]),
             xdata = qz_array,label= 'Corrected' )
        
    legend()
    title('Qz Slice at '+ str(q_slice))
    xlabel('qz(A^-1)')
    ylabel('Normalized Intensity')

    return graph


def data_compare(intensity_one,intensity_two,mins,maxs):
    from pylab import imshow,colorbar,show, title, xlabel, ylabel,subplot
    
    plotxmin = mins[0]
    plotxmax = maxs[0] 
    plotzmin = mins[-1]
    plotzmax = maxs[-1]
    
    intensity_one[isnan(intensity_one)] = 0.0
    intensity_two[isnan(intensity_one)] = 0.0
    
    intensity_one += min(intensity_one[(intensity_one)>0])/2
    intensity_two += min(intensity_one[(intensity_one)>0])/2
    
    vmin = min(log(intensity_one))
    vmax = max(log(intensity_one))
    
    plot_extent = (plotxmin,plotxmax,plotzmin,plotzmax)
    
    subplot(311)
    
    imshow(log(abs(intensity_one.T)),
           aspect='auto',interpolation='nearest',
           extent=plot_extent,origin='lower', 
           vmin = vmin, vmax = vmax)
    colorbar()
    
    subplot(312)
    
    imshow(log(abs(intensity_two.T)),
           aspect='auto',interpolation='nearest',
           extent=plot_extent,origin='lower', 
           vmin = vmin, vmax = vmax)
    colorbar()
    
    subplot(313)
    
    imshow((abs(intensity_one - intensity_two)/intensity_one).T,aspect='auto',
           interpolation='nearest',extent=plot_extent,origin='lower',
            vmin = vmin, vmax = vmax)
    colorbar()
    
    
def test():
    '''
    this test was used to fix the plotters
    '''
    intensity = array([[1,2,3,4,5],[5,4,3,2,1],[1,5,2,4,3],[5,1,4,2,3],[5,1,4,2,3]])
    mins = array([-3,-3,-3])
    maxs = array([3,3,3])
    qz_slice(intensity,mins,maxs)
if __name__=="__main__":test()