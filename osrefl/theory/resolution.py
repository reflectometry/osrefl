# Copyright (C) 2008 NIST Center for Neutron Research
# All rights reserved.
# See LICENSE.txt for details.

#author: Brian Maranville
#Documented: Christopher Metting

from numpy import sin, cos, arctan, arcsin, pi, log, sqrt
from numpy import array, zeros,  ones, histogram2d, cumsum
from numpy import  arange, linspace, meshgrid

from scipy import ndimage, interpolate
from pylab import figure, imshow, colorbar

OVERSAMPLE_U = 3.0
OVERSAMPLE_TH = 2.0

class Resolution_correction():
    '''
    Overview:
        An object for doing a resolution correction. It takes the 2-d Q-space 
    map, with metadata (maximum and minimum values for both axes), delta_theta
    (in degrees),delta_wl, wl (neutron wavelength), and desired angle range
    (degrees) in intermediate map.default angle range is +/- 20 degrees, which
    is larger than most off-specular scans intermediate object is created: data
    is transformed into u and theta (Q) coordinates. Angle resolution of polar
    map is determined by angle resolution of rectilinear Q-space map at highest
    Qz. Q resolution is the same as Qz resolution of rectilinear map, and u is a
    coordinate very similar to Q but stretched so that delta_u is constant
    across the map.
          
    blur is applied (2d gaussian) in which Q resolution is: 
    sqrt((delta_theta/twotheta)^2 + (delta_lambda/lambda)^2)
    and theta resolution is delta_theta.
    
    polar Q-space map is transformed back to rectilinear coordinates for output.
    
    
    Parameters(__init__):
    
    map_in:(float,[]|reflectivity) = The two dimensional Q-space array that that
    is to have the resolution correction applied to.
    
    map_qx_range:(float[2]|angstroms^-1) = The range of qx values that are
    represented by the map_in attribute. This attribute is structured as:
    [minimum qx, maximum qx]
    
    map_qz_range:(float[2]|angstroms^-1) = The range of qz values that are
    represented by the map_in attribute. This attribute is structured as:
    [minimum qz, maximum qz]
    
    delta_theta:(float|degrees) = The beam's angular divergence. When the beam
    is emitted onto the sample, it is not a point but rather a finite spread of 
    probing beam approaching at some finite range of angles. This parameter 
    represents the spread of the incoming beam +/- the measured angle.
    This is an instrument specific parameter.
    
    delta_wl:(float|Angstroms) = The wavelength(energy) divergence of the
    beam.This is an instrument specific parameter. 
    
    wl:(float|Angstroms) = The wavelength(energy) of the incoming beam.This is 
    an instrument specific parameter.
    
    angle_range:(float,[2]|degrees) = Gives a range of angles for the u plot
    to be solved for.
    
    
    Parameters(Class):
    
    u_th_map:(float,[]|map) = The 'map_in' attribute represented in polar
    coordinates.
    
    th_u_blurred:(float,[]|map) = The 'u_th_map' attribute after a convolution
    is applied in the q space and the theta direction of the polar coordinate
    plot. This is an intermediate plot that allows the covelutions in the x and
    y direction of the array to be applied separately.
    
    map_out:(float,[]|map) =  The final convolution of 'map_in'.
    
    Note: 
    - Steps to convolution:
        1) Convert inputed Qx-Qz array to polar coordinates with a theta and qz
        axis.
        
        2) Apply a simple gaussian convolution along the theta axis of the polar
        coordinate map.
        
        3) Apply a simple gaussian convolution along the q axis of the polar
        coordinate map.
        
    - The convolution applied includes wavelength divergence and beam angle
    divergence.

    '''
  
    def __init__(self, map_in, map_qx_range, map_qz_range, delta_theta = None,
                 delta_wl = None, wl = None, angle_range = [-20.0, 20.0]):
        
        self.map_in = map_in
        self.map_qx_range = map_qx_range
        self.map_qz_range = map_qz_range
        
        plot_data = 0
        if delta_theta == None:self.delta_theta = 0.049
        else: self.delta_theta = delta_theta
        
        if delta_wl == None: self.delta_wl = 0.05
        else: self.delta_wl = delta_wl
          
        if wl == None: self.wl = 5.000
        else: self.wl = wl

        self.angle_range = angle_range
        self.polar_map = None
        self.polar_extent = None
        
        self.generate_u_th_map()

        # then blur it in theta (fast)
        self.th_blurred = self.blur_theta(self.u_th_map)

        # blur the resulting map along q
        self.th_u_blurred = self.blur_q(self.th_blurred)

        # now transform back into linear Qx, Qz coordinates:
        self.map_out = self.convert_back(self.th_u_blurred)
        
        extent = array([self.map_qx_range, self.map_qz_range]).flatten()
        if plot_data:
            self.blurred_back = figure()
            imshow(log(self.map_out.T + 1), origin='lower', aspect='auto', 
                   extent = extent)
            colorbar()
            
        return
    
    
    def generate_u_th_map(self):
        
        '''
        Overview:
            This method evaluates spacing at max Qz and produces a scattering 
        map in Q|theta (polar) coordinates where theta is the angle between Q 
        and the Qz axis (theta = 0 when Q = Qz)
        
        
        Notes:
        -Coordinate u has constant error: 
            delta_u = self.delta_q.max() * pixels_per_q
        
        -To transform back to Qx, Qz, divide by delta_q, multiply by sin and cos
        '''
        
        qz_max = array(self.map_qz_range).max()
        qz_min = array(self.map_qz_range).min()
        qx_max = array(self.map_qx_range).max()
        qx_min = array(self.map_qx_range).min()
        
        #protect against divide by zero
        if qz_min == 0 :
            qz_min = 1e-9
        
        self.qx_stepsize = (qx_max - qx_min) / self.map_in.shape[0]
        self.qz_stepsize = (qz_max - qz_min) / self.map_in.shape[1]
        
        x0 = -qx_min / self.qx_stepsize
        z0 = -qz_min / self.qz_stepsize
        
        # making coordinate grid for the polar plot
        self.th = arange(self.angle_range[0], self.angle_range[1], 
                         (1.0/OVERSAMPLE_TH)*
                         arctan(self.qx_stepsize/qz_max) * 180./pi)
        
        self.q = linspace(qz_min, qz_max, self.map_in.shape[1])
        
        self.twoth = arcsin( self.wl * self.q / (4. * pi) ) * 2.0 * 180./ pi
        
        self.delta_q = sqrt((self.delta_wl / self.wl)**2 + 
                            (self.delta_theta / self.twoth)**2) * self.q
                            

        self.u_steps = self.delta_q.max() / self.delta_q
        self.delta_u = self.delta_q.max()
        
        q_u_coords = cumsum(self.u_steps) - self.u_steps[0]
        # generate lookup function to find q(u)
        f = interpolate.interp1d(q_u_coords, self.q)
        u_range = q_u_coords[-1]
        self.u = linspace(0,int(u_range), int(u_range) * OVERSAMPLE_U)
    
        # now we have an expanded set of q values (spaced according to error)
        self.new_q = f(self.u)
        xy = meshgrid(self.th,self.new_q)

        qzp = xy[1] * cos(xy[0]*pi/180.)
        qxp = xy[1] * sin(xy[0]*pi/180.)
        
        self.qxpp = qxp / self.qx_stepsize + x0
        self.qzpp = qzp / self.qz_stepsize + z0
        
        self.qxpp.clip(0)
        self.qzpp.clip(0)

        newim = ndimage.map_coordinates(
                        self.map_in, array([self.qxpp,self.qzpp]), order = 1)
        
        self.u_th_extent = array([self.angle_range[0], self.angle_range[1], 
                                  self.u.min(), self.u.max()])
        
        
        self.u_th_map = newim
        self.fig_u_th = figure()
        
        return
    

    
    
    def blur_theta(self, data):
        '''
        Overview:
            The delta_theta was obtained from the input parameters. Now in polar
        coordinates.This is a simple 1-d gaussian convolution.
        
        '''
        
        conv_factor = 2 * sqrt(2 * log(2))
        sigma_th = self.delta_theta / conv_factor
        
        # what is the width in pixel space?  angle-to-pixel conversion:
        pixels_per_degree = (self.u_th_map.shape[1] / 
        ( array(self.angle_range).max() - array(self.angle_range).min() ))
        
        sigma_x = sigma_th * pixels_per_degree
        
        th_blurred_map = ndimage.gaussian_filter1d(data, sigma_x, axis=1)
        
        return th_blurred_map
  
  
    def blur_q(self, data):
        '''
        Overview:
            Applies the 1-d gaussian convolution in the q space direction. This
        is applied to the attribute that is in polar coordinates.
        '''
        conv_factor = 2 * sqrt(2 * log(2))
        sigma_u = self.delta_u / conv_factor
        
        # what is the width in pixel space?  u-to-pixel conversion:
        pixels_per_u = (self.u_th_map.shape[0] / 
        ( array(self.map_qz_range).max() - array(self.map_qz_range).min() ))
        
        sigma_y = sigma_u * pixels_per_u
        
        q_blurred_map = ndimage.gaussian_filter1d(data, sigma_y, axis=0)
        
        return q_blurred_map



    def convert_back(self, data):
        '''
        Overview:
            Converts the polar coordinate space attribute back to the starting
        rectilinear coordinates which is the the coordinate system used by the
        inputed data.
        
        Notes:
        -histogram2d returns the histogram as well as the edge values for the
        x and y axis. Only the histogram result is used in this conversion.
        '''

        # here's the coordinates in rectilinear of the polar stuff
        qxpp = self.qxpp
        qzpp = self.qzpp
        select_matrix = ( (qxpp < self.map_in.shape[0]) * (qxpp >= 0) * 
                          (qzpp < self.map_in.shape[1]) * (qzpp >= 0) )
        
        xpp = qxpp[select_matrix]
        zpp = qzpp[select_matrix]
        xBin = xpp.astype(int)
        zBin = zpp.astype(int)
        intens = data[select_matrix]
        
        outshape = self.map_in.shape
        output_bins = zeros((self.map_in.shape[0], self.map_in.shape[1], 3))
        
        if (len(xBin) > 0):
            
            hist2d = histogram2d(
                     xBin,zBin, bins = (outshape[0],outshape[1]), 
                     range=((0,outshape[0]),(0,outshape[1])), weights=intens)
            
            output_bins[:,:,0] += hist2d[0]
            
            hist2d = histogram2d(
                     xBin,zBin, bins = (outshape[0],outshape[1]), 
                     range=((0,outshape[0]),(0,outshape[1])), 
                     weights=ones(len(intens)) )
            
            output_bins[:,:,1] += hist2d[0]
            

        nonzero_mask = (output_bins[:,:,1] != 0)
        output_bins[:,:,2][nonzero_mask] = ( output_bins[:,:,0][nonzero_mask] / 
                                             output_bins[:,:,1][nonzero_mask] )

        return output_bins[:,:,2]
    
        
    def display(self, data):
        '''
        Overview:
            A simply way to view any of the array maps produced at any step in
        in the convolution process. This is available for testing purposes and
        will not generally be used by consumers of this module.
        
        
        '''
        self.fig_polar = figure()
        imshow(log(data + 1.), aspect = 'auto', origin='lower', 
               extent = self.polar_extent)
        return
def main():
    import numpy,pylab
    
    #Point count can not be altered. From loaded data file.
    points = numpy.array([50,60])
    mins = numpy.array([-.001,0.0002])
    maxs = numpy.array([.001,0.3])

    data = numpy.fromfile('BA_data')
    data = data.reshape(50,60)
    
    qx = numpy.linspace(mins[0],maxs[0],points[0])
    qz = numpy.linspace(mins[1],maxs[1],points[1])
    
    qx = numpy.reshape(qx,[points[0],1])
    qz = numpy.reshape(qz,[1,points[1]])
    
    pyRC = Resolution_correction(data, [mins[0],maxs[0]], [mins[1],maxs[1]], 
                                 delta_theta = None,
                 delta_wl = None, wl = None, angle_range = [-20.0, 20.0])
    figure(0)
    pylab.imshow(numpy.log(numpy.fliplr(data).T),extent = [mins[0],maxs[0],mins[1],maxs[1]],aspect='auto')
    figure(1)
    pylab.imshow(numpy.log(numpy.fliplr(pyRC.map_out).T),extent = [mins[0],maxs[0],mins[1],maxs[1]],aspect='auto')
    pylab.colorbar()
    pylab.show()
    
    return
if __name__ == "__main__":
    main()
