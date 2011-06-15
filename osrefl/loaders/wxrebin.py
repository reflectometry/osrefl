#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# check to see if all parameters are set on the command line
# if they are, then don't open GUI interface
from math import *
import sys
import os, pickle
#from Tkinter import *
#from .reduction import *
from .reduction import load as icp_load
from numpy import *
from pylab import imshow,cm,colorbar,hot,show,xlabel,ylabel
#from qxqz_data_array_class import *
from matplotlib.widgets import  RectangleSelector

# Note that the following import can be done three ways: as as standard import,
# as a relative import, or as an import that explicitly references the root
# directory of the package.
#from binned_data import plottable_2d_data
from .binned_data import plottable_2d_data
#from osrefl.loaders.binned_data import plottable_2d_data

import __main__
#global rebinned_data_objects
#rebinned_data_objects = []
import wx
app = wx.App()
app.MainLoop()

__main__._rebinned_data_objects = []

class Supervisor():
    """ class to hold rebinned_data objects and increment their reference count """
    def __init__(self):
        self.refcount = 0
        self.rebinned_data_objects = []

    def __iadd__(self, new_rebinned_data_object):
        if isinstance(new_rebinned_data_object, rebinned_data):
            self.refcount += 1
            self.rebinned_data_objects.append(new_rebinned_data_object)
        print('adding to myself')
        return self

supervisor = Supervisor()

ID_EXIT=102


class rebinned_data():
    """ class to import ICP files, and rebin them according to
    user-supplied parameters.  A GUI is opened when an instance
    of the class is created, and parameters are entered there.
    Add, subtract and divide functions work with other instances,
    as long as the array dimensions are equal, so use with care. """

    def __init__( self, filenames = [], plot_data = True, normalization = 'monitor' ):
        self.detectorMaxXPixel = 608.0
        self.pixelsPerDegree = 80.0
        self.inFileNames = filenames
        self.inFileObjects = []
        self.theta_offset = 0.0
        self.twoth_zero_pixel = 309.0
        self.normalization = normalization
        self.qxqz_2d_data = None
        self.th2th_data = None
        self.th_in_th_out_data = None
        self.supervisor = None


        # ask which filenames to work with
        #self.inFileNames = self.ChooseInputFiles()
        if (len(filenames) < 1):
            self.get_load_params_gui(None, -1, 'Select Data Files', self)
        print 'names: ' + str(self.inFileNames)
        #self.get_filenames_gui(None, -1, 'Get filenames to work with')

        self.number = 0
        global supervisor
        print supervisor
        if supervisor:
            supervisor += self
            self.number = supervisor.refcount
            self.supervisor = supervisor

        for dataFileName in self.inFileNames:
            self.inFileObjects.append(icp_load(dataFileName))
            print('datafile: ' + dataFileName + ' loaded.')

        #grab wavelength from first file:
        self.wavelength = self.inFileObjects[0].detector.wavelength
        #grab description from first file:
        self.description = self.inFileObjects[0].description
        # then make those filenames into 2th-th map
        self.make_2th_th_map()
        if plot_data:
            self.th2th_data.wxplot()

    def motors(self, v):
        start = double(v[0])
        stop = double(v[-1])
        step = double((stop-start)/(len(v)-1))
        return start,step,stop

    def getTwoTheta (self, pixel, angle4 ):
        return ( ( (self.twoth_zero_pixel - pixel ) / self.pixelsPerDegree ) + angle4 )

    def twoThetaFromPixel(self, pixel, angle):
        return_angle = (self.params['qZeroPixel'] - pixel) / self.pixelsPerDegree
        return_angle += (angle - self.params['qZeroAngle'])
        return return_angle

    def setThetaOffset(self, offset):
        self.theta_offset = offset
        self.make_2th_th_map()
        if self.qxqz_2d_data:
            self.runConversion(self.qxqz_2d_data.params)

    def correctionFromPixel(self, xpixel, subDim=1):
        pixelsPerDegree = 80.0 * subDim
        wiggleAmplitude = 0.15
        pixelCorrection = ( (32.0 / (2.0 * pi) ) * wiggleAmplitude * sin( 2.0 * pi * xpixel / 32.0 ) )
        widthCorrection = ( wiggleAmplitude * cos( 2.0 * pi * xpixel / 32.0 ) )
        return [widthCorrection, pixelCorrection]

    def getQxQzI (self, A3, A4, x, rawI ):
        qLength = 2.0 * pi / self.wavelength
        rawTwoTheta = self.getTwoTheta( x, A4 )
        widthCorrection,angleCorrection = self.correctionFromPixel(x)
        correctedI = rawI / ( 1.0 + widthCorrection )
        twoTheta = rawTwoTheta + angleCorrection
        tilt = A3 - ( twoTheta / 2.0 )
        dq = 2.0 * qLength * sin( ( pi / 180.0 ) * ( twoTheta / 2.0 ) )
        qxOut = dq * sin( pi * tilt / 180.0 )
        qzOut = dq * cos( pi * tilt / 180.0 )
        return [qxOut, qzOut, correctedI]

    def runConversion (self, params_in={}, plot_result=True):
        """convert 2th-th dataset to Qx,Qz dataset"""

        # initialize the parameters that follow the dataset around:
        print 'beginning rebin:\n'
        params = {
          'description': '',
          'x_max': 0.003,
          'x_min': -0.003,
          'x_steps': 200,
          'y_max': 0.10,
          'y_min': 0.0,
          'y_steps': 200,
          'x_units': 'Qx (inv Angstroms)',
          'y_units': 'Qz (inv Angstroms)',
          'qZeroAngle': 0.0,
          'qZeroPixel': 309.0
          }

        # override default parameters with saved parameters from last time
        stored_name = '.stored_rebin_params'
        if os.access(stored_name, os.R_OK):
            stored_file = open(stored_name, 'r')
            stored_params = pickle.load(stored_file)
            stored_file.close()
            for label in stored_params.keys():
                if params.has_key(label):
                    params[label] = stored_params[label]

        for label in params_in.keys():
            if params.has_key(label):
                params[label] = params_in[label]

        if (len(params_in) < len(params)):
            s = self.get_params_gui(None,-1,'input Qx-Qz conversion parameters',params)
            params['x_steps'] = int(round(params['x_steps']))
            params['y_steps'] = int(round(params['y_steps']))

        #self.bin_data = zeros((params['x_steps'],params['y_steps'],4))
        #self.qxqz_data_array = plottable_2d_data(self.bin_data, params, self)
        #qZMax = self.params['qZMax']
        #qZMin = self.params['qZMin']
        qZSteps = params['y_steps']
        #qXMax = self.params['qXMax']
        #qXMin = self.params['qXMin']
        qXSteps = params['x_steps']
        #qzStepSize = ( qZMax - qZMin ) / qZSteps
        #qxStepSize = ( qXMax - qXMin ) / qXSteps

        #self.qxqz_array = self.convert_to_qxqz(self.twoth_th_array)
        self.qxqz_arrays = []

        #for dataFileObject in self.inFileObjects:
            #qxqz_array = self.convert_to_qxqz(dataFileObject)
            #self.qxqz_arrays.append(qxqz_array)
            #print('loaded: ' + dataFileObject.name)

        if (params['description'] == ''):
            params['description'] = self.description
            #params['description'] = self.inFileObjects[0].description

        qZSteps = params['y_steps']
        qXSteps = params['x_steps']
        self.bin_data = zeros((params['x_steps'],params['y_steps'],4))

        qxqz_arrays = [self.convert_2thth_qxqz(self.th2th_data)]

        #for dataFileObject in self.inFileObjects:
            #qxqz_array = self.convert_to_qxqz(dataFileObject)
            #self.qxqz_arrays.append(qxqz_array)
            #print('converted: ' + dataFileObject.name)

        for qxqz_array in qxqz_arrays:
            self.do_rebinning(qxqz_array, params, self.bin_data)
            print('rebinned.')

        #self.do_rebinning()

        #self.minimum_intensity = inf
        for j in range(qZSteps):
            for i in range(qXSteps):
                pixelCount = self.bin_data[i,j,1]
                if ( pixelCount > 0 ):
                    monitorTotal = self.bin_data[i,j,2]
                    avg = self.bin_data[i,j,0] / double(monitorTotal)
                    self.bin_data[i,j,3] = avg
                    #if (avg < self.minimum_intensity and avg > 0):
                      #self.minimum_intensity = avg

        qxqz_2d_params = params

        #self.make_2th_th_map()
        zero_pixels = (self.bin_data[:,:,1]==0)
        self.bin_data[zero_pixels] = self.reverse_lookup(params)[zero_pixels]
        # make sure that these pixels still show up as empty for sums:
        self.bin_data[zero_pixels][:,1] = 0
        title = 'dataset ' + str(self.number) + ': qxqz_2d_data'
        if self.qxqz_2d_data == None:
            self.qxqz_2d_data = plottable_2d_data(self.bin_data, qxqz_2d_params, self, title = title)
        else:
            self.qxqz_2d_data.bin_data = self.bin_data.copy()
            self.qxqz_2d_data.params = qxqz_2d_params
            self.qxqz_2d_data.title = title
        if plot_result:
            self.qxqz_2d_data.wxplot()

    def interpolate_vec(self,vec,n=1):
        if n==1:
            return vec
        else:
            new_indices = linspace(0.,float(len(vec))-1,(len(vec) * n)-1)
            mod = fmod(new_indices,1)
            new_vec = vec[floor(new_indices).astype(int)] * (1-mod) + vec[ceil(new_indices).astype(int)] * mod
            return new_vec

    def expand_arr(self,arr,n=1):
        if n==1:
            return arr
        else:
            new_indices_1 = array([arange(arr.shape[0]*n -1)]).T / n
            new_indices_2 = array([arange(arr.shape[1]*n -1)]) / n
            new_arr = arr[new_indices_1,new_indices_2]
            return new_arr

    def convert_2thth_qxqz(self, twoth_th_2d_data):
        qLength = 2.0 * pi / self.wavelength
        data = twoth_th_2d_data
        corrected_I = data.bin_data[:,:,0]
        monitor_array = data.bin_data[:,:,2]
        twoth_stepsize = (data.params['x_max'] - data.params['x_min']) / (data.params['x_steps'] - 1)
        twoTheta_array = arange(data.params['x_steps']) * twoth_stepsize + data.params['x_min']
        twoTheta_array = array([twoTheta_array]).T
        th_stepsize = (data.params['y_max'] - data.params['y_min']) / (data.params['y_steps'] - 1)
        th_array = arange(data.params['y_steps']) * th_stepsize + data.params['y_min']
        th_array = array([th_array])
        tilt_array = th_array - ( twoTheta_array / 2.0 )
        self.twoTheta_array = twoTheta_array
        self.tilt_array = tilt_array
        qxOut = 2.0 * qLength * sin( ( pi / 180.0 ) * ( twoTheta_array / 2.0 ) ) * sin( pi * tilt_array / 180.0 )
        qzOut = 2.0 * qLength * sin( ( pi / 180.0 ) * ( twoTheta_array / 2.0 ) ) * cos( pi * tilt_array / 180.0 )
        self.qxOut = qxOut
        self.qzOut = qzOut
        qxqz_data = array([qxOut,qzOut,corrected_I,monitor_array])

        return qxqz_data

    def convert_to_qxqz(self, inData):
        monitor_counts = inData.display_monitor
        self.wavelength = inData.detector.wavelength
        qLength = 2.0 * pi / self.wavelength
        A3_start, A3_step, A3_stop = self.motors(inData.sample.angle_x + self.theta_offset)
        A4_start, A4_step, A4_stop = self.motors(inData.detector.angle_x)

        #qZeroPixel = self.params['qZeroPixel']
        #qZeroAngle = self.params['qZeroAngle']
        #subDim = self.params['subDim']
        subDim = 1

        pixel_vector = self.interpolate_vec(arange(608.0),subDim)
        intens_correction, pixel_correction = self.correctionFromPixel(pixel_vector)
        corrected_pixel = pixel_vector + pixel_correction
        A4_vector = array([self.interpolate_vec(inData.detector.angle_x,subDim)]).T
        A3_vector = array([self.interpolate_vec(inData.sample.angle_x + self.theta_offset,subDim)]).T
        monitor_vector = array([self.interpolate_vec(inData.monitor.counts,subDim)]).T
        twoTheta_array = self.getTwoTheta(corrected_pixel,A4_vector)

        intens = self.expand_arr(inData.detector.counts,subDim)
        print 'intens.shape=', intens.shape
        print 'intens_correction.shape=',intens_correction.shape
        corrected_I = intens / (1.0 + intens_correction)
        tilt_array = A3_vector - ( twoTheta_array / 2.0 )
        qLength = 2.0 * pi / self.wavelength
        # dq = 2.0 * qLength * sin( ( pi / 180.0 ) * ( twoTheta / 2.0 ) )
        qxOut = 2.0 * qLength * sin( ( pi / 180.0 ) * ( twoTheta_array / 2.0 ) ) * sin( pi * tilt_array / 180.0 )
        qzOut = 2.0 * qLength * sin( ( pi / 180.0 ) * ( twoTheta_array / 2.0 ) ) * cos( pi * tilt_array / 180.0 )
        self.qxOut = qxOut
        self.qzOut = qzOut
        monitor_array = monitor_vector + zeros(shape(corrected_I))
        qxqz_data = array([qxOut,qzOut,corrected_I,monitor_array])

        return qxqz_data

    def convert_to_th2th(self, qx, qz):
        Q = sqrt(qx * qx + qz * qz) * sign(qz)
        twoth = arcsin( self.wavelength * Q / (4. * pi) ) * 2.0 * 180./ pi
        tilt = arctan( qx / qz ) * 180. / pi
        th = tilt + twoth / 2.0
        self.twoth = twoth
        self.th = th
        return th, twoth

    def do_rebinning(self, qxqz_data, params, output_bins):
        print 'rebinning, oy\n'
        qZMax = params['y_max']
        qZMin = params['y_min']
        qZSteps = params['y_steps']
        qXMax = params['x_max']
        qXMin = params['x_min']
        qXSteps = params['x_steps']
        qzStepSize = ( qZMax - qZMin ) / qZSteps
        qxStepSize = ( qXMax - qXMin ) / qXSteps
        subDim = 1

        d = shape(qxqz_data)

        select_matrix = (qxqz_data[0] < qXMax) * (qxqz_data[0] >= qXMin) * (qxqz_data[1] < qZMax) * (qxqz_data[1] >= qZMin)

        qx = qxqz_data[0][select_matrix]
        self.data_qx_list = qx.copy()
        qz = qxqz_data[1][select_matrix]
        self.data_qz_list = qz.copy()
        intens = qxqz_data[2][select_matrix]
        monitor = qxqz_data[3][select_matrix]
        xBin = floor( (qx - qXMin ) / qxStepSize ).astype(int32)
        zBin = floor( (qz - qZMin ) / qzStepSize ).astype(int32)

        outshape = output_bins.shape
        if (len(xBin) > 0):
            hist2d, xedges, yedges = histogram2d(xBin,zBin, bins = (outshape[0],outshape[1]), range=((0,outshape[0]),(0,outshape[1])), weights=intens)
            output_bins[:,:,0] += hist2d
            hist2d, xedges, yedges = histogram2d(xBin,zBin, bins = (outshape[0],outshape[1]), range=((0,outshape[0]),(0,outshape[1])), weights=ones(len(intens)) )
            output_bins[:,:,1] += hist2d
            hist2d, xedges, yedges = histogram2d(xBin,zBin, bins = (outshape[0],outshape[1]), range=((0,outshape[0]),(0,outshape[1])), weights=monitor)
            output_bins[:,:,2] += hist2d

        #for i in range(len(xBin)):
            #output_bins[xBin[i],zBin[i],0] += intens[i]
            #output_bins[xBin[i],zBin[i],1] += 1
            #output_bins[xBin[i],zBin[i],2] += monitor[i]

        return

    def reverse_lookup(self, params):
        base = zeros((params['x_steps'], params['y_steps']))
        qx = linspace(params['x_min'], params['x_max'], params['x_steps'])
        qx.shape = (qx.shape[0], 1)
        self.qx = qx + base
        qz = linspace(params['y_min'], params['y_max'], params['y_steps'])
        qz.shape = (1, qz.shape[0])
        self.qz = qz + base
        print self.qx.shape, self.qz.shape
        th, twoth = self.convert_to_th2th(self.qx, self.qz)
        #self.th = th
        #self.twoth = twoth
        print th.shape, twoth.shape
        #select_matrix = (th >= self.th_min_min and th <= self.th_max_max and twoth >= self.twoth_min_min and twoth <= self.twoth_max_max)
        th_index = (floor((th - self.th_min_min)/self.th_stepsize)).astype(int32)
        twoth_index = (floor((twoth - self.twoth_min_min)/self.twoth_stepsize)).astype(int32)

        th_index_min = 0
        th_index_max = self.th_2th_array[:,:,3].shape[0]
        twoth_index_min = 0
        twoth_index_max = self.th_2th_array[:,:,3].shape[1]

        new_th_index = th_index.copy()
        # map all indices that are outside the proper range to 0 index
        new_th_index[new_th_index < th_index_min] = 0
        new_th_index[new_th_index >= th_index_max] = 0
        new_twoth_index = twoth_index.copy()
        # map all indices that are outside the proper range to 0 index
        new_twoth_index[new_twoth_index < twoth_index_min] = 0
        new_twoth_index[new_twoth_index >= twoth_index_max] = 0
        #self.th_index = new_th_index
        #self.twoth_index = new_twoth_index

        # the action happens here:  the output array is filled with data with lookup indices
        # specified by th_index and twoth_index
        out_array = self.th_2th_array[new_th_index, new_twoth_index]

        # here we set all the data that was outside the mapped region in th,twoth to zero
        out_array[th_index < 0] = 0
        out_array[th_index >= th_index_max] = 0
        out_array[twoth_index < 0] = 0
        out_array[twoth_index >= twoth_index_max] = 0

        return out_array

    def make_2th_th_map(self):
        """ makes a composite array that contains all the data objects
        in a regular theta-2theta grid.  For theta, find stepsize in all
        data sets, then use largest step size as grid spacing.  In 2theta,
        somewhat simpler: use the pixel size as the grid spacing (1/80 deg)
        """
        from reduction import rebin as reb
        instr_resolution = 0.00001 # difference that is negligible for all motors

        num_files = len( self.inFileObjects )
        th_min = zeros(num_files)
        th_max = zeros(num_files)
        th_len = zeros(num_files)
        th_step = zeros(num_files)
        twoth_min = zeros(num_files)
        twoth_max = zeros(num_files)
        twoth_len = zeros(num_files)
        twoth_step = zeros(num_files)

        pixel_vector = arange(609.0)
        intens_correction, pixel_correction = self.correctionFromPixel(pixel_vector)
        corrected_pixel = pixel_vector + pixel_correction

        for i in range(num_files):
            th_min[i] = self.inFileObjects[i].sample.angle_x.min() + self.theta_offset
            th_max[i] = self.inFileObjects[i].sample.angle_x.max() + self.theta_offset
            th_len[i] = self.inFileObjects[i].sample.angle_x.shape[0]

            A4_vector = self.inFileObjects[i].detector.angle_x
            twoth_vector = self.getTwoTheta(corrected_pixel, array([A4_vector]).T)
            twoth_max[i] = twoth_vector.max()
            twoth_min[i] = twoth_vector.min()

            if th_len[i] > 1:
                th_step[i] = float(th_max[i] - th_min[i]) / th_len[i]
            else:
                th_step[i] = 0.0

            #if twoth_len[i] > 1:
                #twoth_step[i] = float(twoth_max[i] - twoth_min[i]) / twoth_len[i]
            #else:
                #twoth_step[i] = 0.0

        # take maximum A3 spacing as grid spacing for theta:
        th_stepsize = th_step.max()
        self.th_stepsize = th_stepsize
        # take the extrema as the limits of our new grid:
        th_min_min = th_min.min()
        self.th_min_min = th_min_min
        th_max_max = th_max.max()
        self.th_max_max = th_max_max
        th_steps = int(float(th_max_max - th_min_min) / th_stepsize)
        self.th_steps = th_steps

        # take 1/80 deg as grid spacing for 2theta:
        twoth_stepsize = 1.0/80.0
        self.twoth_stepsize = twoth_stepsize
        # take the extrema as the limits of our new grid:
        twoth_min_min = twoth_min.min()
        self.twoth_min_min = twoth_min_min
        twoth_max_max = twoth_max.max()
        self.twoth_max_max = twoth_max_max
        twoth_steps = int(round(float( twoth_max_max - twoth_min_min ) / twoth_stepsize))
        self.twoth_steps = twoth_steps

        th_twoth_data_extent = (twoth_min_min,twoth_max_max,th_min_min,th_max_max)
        #th_twoth_data = zeros((th_steps, twoth_steps, 4))
        print th_min_min, th_max_max, th_steps
        print twoth_min_min, twoth_max_max, twoth_steps

        twoth_bin_edges = (arange(twoth_steps + 1)*twoth_stepsize) + twoth_min_min
        #self.twoth_bin_edges = twoth_bin_edges
        th_bin_edges = (arange(th_steps + 1)*th_stepsize) + th_min_min
        #self.th_bin_edges = th_bin_edges

        th_2th_array = zeros((th_steps, twoth_steps, 4))
        for i in range(num_files):
            inData = self.inFileObjects[i]
            #pixel_vector = arange(609.0,0.0,-1.0)
            pixel_vector = arange(609.)
            intens_correction, pixel_correction = self.correctionFromPixel(pixel_vector)
            corrected_pixel = pixel_vector + pixel_correction
            intens = inData.detector.counts
            #cut off the last bit of intens_correction (which is larger by one than intens, because
            #we're making edges
            intens_correction = intens_correction[:-1]
            corrected_I = intens / (1.0 + intens_correction)
            #corrected_I = intens
            th = inData.sample.angle_x + self.theta_offset
            # need one more bin edge than data points:
            th_vector = zeros(len(th) + 1)
            th_vector[:len(th)] = th
            # this fails if there's only 1 theta value for the data file
            th_vector[-1] = th_vector[-2] + (th_vector[-2] - th_vector[-3])
            print th_vector
            th_vector_comp = (arange(th_len[i] + 1.0)*th_step[i]) + th_min[i]
            print th_vector_comp
            #self.th_vector = th_vector
            #th_array = array([th_vector]).T + zeros(shape(intens))
            # if the normalization is chosen to be vs. time instead of vs. monitor counts, the change is
            # made here:
            if self.normalization == 'time':
                monitor_vector = inData.monitor.count_time
            else: # self.normalization = monitor by default
                monitor_vector = inData.monitor.counts

            monitor_array = array([monitor_vector]).T + zeros(shape(intens))
            pixel_array = ones(shape(intens))

            A4_vector = inData.detector.angle_x
            if ( (A4_vector.max() - A4_vector.min() ) < instr_resolution ):
                #then the detector is fixed and we can pass a single 2theta vector to rebin2d
                twoth_vector = self.getTwoTheta(corrected_pixel,A4_vector.min())
                #self.twoth_vector = twoth_vector

                new_data = reb.rebin2d(th_vector,twoth_vector,corrected_I.astype('float64'),th_bin_edges,twoth_bin_edges)
                #print th_vector.shape, twoth_vector.shape, corrected_I.shape, th_bin_edges.shape, twoth_bin_edges.shape
                #new_data = reb.rebin2d(twoth_vector,th_vector,intens,twoth_bin_edges,th_bin_edges)
                #print new_data
                new_pixelcount = reb.rebin2d(th_vector,twoth_vector,pixel_array.astype('float64'),th_bin_edges,twoth_bin_edges)
                new_mon = reb.rebin2d(th_vector,twoth_vector,monitor_array.astype('float64'),th_bin_edges,twoth_bin_edges)
                #new_pixelcount = reb.rebin2d(twoth_vector,th_vector,pixel_array,twoth_bin_edges,th_bin_edges)
                #new_mon = reb.rebin2d(twoth_vector,th_vector,monitor_array,twoth_bin_edges,th_bin_edges)
                th_2th_array[:,:,0] += new_data
                th_2th_array[:,:,1] += new_pixelcount
                th_2th_array[:,:,2] += new_mon
            else:
                #then the detector is not fixed, and we have to pass in each A4 value at a time to rebin
                for j in range(th_vector.shape[0]-1):
                    twoth_vector = self.getTwoTheta(corrected_pixel,A4_vector[j])
                    #print th_vector.shape, twoth_vector.shape, corrected_I.shape, th_bin_edges.shape, twoth_bin_edges.shape
                    new_data = reb.rebin(twoth_vector,corrected_I[j].astype('float64'),twoth_bin_edges)
                    #new_data = reb.rebin2d(twoth_vector,th_vector[j],corrected_I[:,j],twoth_bin_edges,th_bin_edges)
                    new_pixelcount = reb.rebin(twoth_vector,pixel_array[j].astype('float64'),twoth_bin_edges)
                    new_mon = reb.rebin(twoth_vector,monitor_array[j].astype('float64'),twoth_bin_edges)
                    print j
                    th_2th_array[j,:,0] += new_data
                    th_2th_array[j,:,1] += new_pixelcount
                    th_2th_array[j,:,2] += new_mon
            print('loaded: ' + inData.name)


        nonzero_pixels = (th_2th_array[:,:,1] > 0)
        monitor_total = th_2th_array[:,:,2][nonzero_pixels]
        avg = th_2th_array[:,:,0][nonzero_pixels] / double(monitor_total)
        th_2th_array[:,:,3][nonzero_pixels] = avg
        self.th_2th_array = th_2th_array


        th2th_params = {
          'description': self.inFileObjects[0].description,
          'x_max': twoth_min_min,
          'x_min': twoth_max_max,
          'x_steps': twoth_steps,
          'y_max': th_max_max,
          'y_min': th_min_min,
          'y_steps': th_steps,
          'x_units': '$2\\theta ({}^{\circ})$',
          'y_units': '$\\theta ({}^{\circ})$'
          }

        twoth_th_array = flipud(th_2th_array.swapaxes(0,1))
        # cut off stuff that should really be zero - bad fp stuff
        twoth_th_array[:,:,3][twoth_th_array[:,:,3] < 1e-16] = 0.
        title = 'dataset ' + str(self.number) + ':  th2th_data'
        self.th2th_data = plottable_2d_data(twoth_th_array, th2th_params, self, title = title)
        return

    def make_th_in_th_out_map(self):
        if (not self.th2th_data):
            print('make_2th_th_map() must be performed first')
            return
        th2th_params = self.th2th_data.params
        th_steps = th2th_params['y_steps']
        th_max = th2th_params['y_max']
        th_min = th2th_params['y_min']
        th_stepsize = float(th_max - th_min)/th_steps
        th_in = arange(th_steps, dtype = 'float') * th_stepsize + th_min
        twoth_steps = th2th_params['x_steps']
        twoth_max = th2th_params['x_max']
        twoth_min = th2th_params['x_min']
        twoth_stepsize = float(twoth_max - twoth_min)/twoth_steps
        twoth = arange(twoth_steps, dtype = 'float') * twoth_stepsize + twoth_min
        #twoth = arange(self.twoth_steps, dtype = 'float') * self.twoth_stepsize + self.twoth_min_min
        th_out_max = twoth_max - th_min
        th_out_min = twoth_min - th_max
        from scipy import ndimage as nd
        tthdata = th2th_data.bin_data
        affine_transform = array([[1.0, -th_stepsize / twoth_stepsize],[0.,1.]])
        th_out_steps = int((th_max - th_min) / twoth_stepsize + twoth_steps)
        th_in_th_out = zeros((th_out_steps, th_steps,4))

        th_in_th_out[:,:,0] = nd.affine_transform(tthdata[:,:,0], affine_transform, offset = 0.0, output_shape=[th_out_steps,th_steps] )
        th_in_th_out[:,:,1] = nd.affine_transform(tthdata[:,:,1], affine_transform, offset = 0.0, output_shape=[th_out_steps,th_steps] )
        th_in_th_out[:,:,2] = nd.affine_transform(tthdata[:,:,2], affine_transform, offset = 0.0, output_shape=[th_out_steps,th_steps] )
        th_in_th_out[:,:,3] = nd.affine_transform(tthdata[:,:,3], affine_transform, offset = 0.0, output_shape=[th_out_steps,th_steps] )
        print th_in_th_out.shape
        print th_out_max, th_out_min
        th_in_th_out_params = {
          'description': self.description,
          'x_max': th_out_max,
          'x_min': th_out_min,
          'x_steps': th_out_steps,
          'y_max': th_max,
          'y_min': th_min,
          'y_steps': th_steps,
          'x_units': '$\\theta_{\\rm{out}} ({}^{\circ})$',
          'y_units': '$\\theta_{\\rm{in}} ({}^{\circ})$'
          }

        th_in_th_out = flipud(self.th_in_th_out)
        # cut off stuff that should really be zero - bad fp stuff
        th_in_th_out[:,:,3][th_in_th_out[:,:,3] < 1e-16] = 0.
        self.th_in_th_out_data = plottable_2d_data(th_in_th_out, th_in_th_out_params, self)
        return

    def do_born_calc(self, params, density_array, density_coords ):
        # density_coords is list [[xmin, xmax], [zmin, zmax]]
        [[xmin,xmax],[zmin,zmax]] = density_coords
        xStep = float((xmax-xmin)) / shape(density_array)[0]
        zStep = float((zmax-zmin)) / shape(density_array)[1]
        # make sure the coordinate matrix is in floating-point, not integers
        x,z = indices(shape(density_array)).astype(float)
        # put the marker in the center of the pixel
        x += 0.5
        z += 0.5
        # fill the coordinate matrix with x and z values
        x *= xStep
        z *= zStep
        # make offset arrays:
        x_sub1 = x.copy()
        x_sub1[1:] = x[:-1]
        x_sub1[0] = x[-1]

        z_sub1 = z.copy()
        z_sub1.T[1:] = z.T[:-1]
        z_sub1.T[0] = z.T[-1]

        qZMax = params['qZMax']
        qZMin = params['qZMin']
        qZSteps = params['qZSteps']
        qXMax = params['qXMax']
        qXMin = params['qXMin']
        qXSteps = params['qXSteps']
        qzStepSize = ( qZMax - qZMin ) / qZSteps
        qxStepSize = ( qXMax - qXMin ) / qXSteps

        Qx,Qz = indices((qXSteps,qZSteps)).astype(float)
        Qx += 0.5
        Qz += 0.5
        Qx *= qxStepSize
        Qx += qXMin
        Qz *= qzStepSize
        Qz += qZMin
        # now have complete matrices of Q values...

        def fourier_sum(qx,qz):
            array_to_sum = density_array * (exp(1j * qz * z) - exp(1j * qz * z_sub1)) * \
            (exp(1j * qx * x) - exp(1j * qx * x_sub1) )
            return array_to_sum.sum()

        reflectivity_matrix = fourier_sum(Qx,Qz)

        prefactor = exp(0.5j * ((Mz-1.0)*Qz*zmax + (Mx-1.0)*Qx*xmax) )
        structure_factor = ( (sin(Mz*Qz*zmax/2.0) / sin(Qz*zmax/2.0)) * (sin(Mx*Qx*xmax/2.0)/sin(Qx*xmax/2.0)) )

        return reflectivity_matrix * prefactor * structure_factor


    def ChooseInputFiles(self,event = None):
        filenames_out = []
        #dlg = wx.FileDialog(None, "Choose a file", '', "", "I*.*", wx.FD_MULTIPLE)
        dlg = wx.FileDialog(None,
                            message="Choose Files",
                            defaultDir=os.getcwd(),
                            defaultFile="",
                            wildcard="I*.*",
                            style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR)
        if dlg.ShowModal() == wx.ID_OK:
            short_filenames=dlg.GetFilenames()
            files_dir = dlg.GetDirectory()
        dlg.Destroy()
        for fn in short_filenames:
            filenames_out.append(files_dir + '/' + fn)
        return filenames_out


    class get_load_params_gui(wx.Dialog):
        def __init__(self, parent, id, title, caller):
        # initialize the outer namespace here:
            self.caller = caller
            wx.Dialog.__init__(self, parent, wx.ID_ANY, title, size=(600,300),
              style = wx.DEFAULT_DIALOG_STYLE | wx.NO_FULL_REPAINT_ON_RESIZE)

            #-- I think you need a panel if you're going to have more than one control in your frame.
            panel = wx.Panel(self, -1)


            #-- create the file selection activator button
            btn_fileselector = wx.Button(panel, -1, "Choose input &files: ")
            self.Bind(wx.EVT_BUTTON, self.ChooseInputFile, btn_fileselector)

            #-- Create the processing button, add it to the panel and wire it up to a function in the class
            btn_SaveExit = wx.Button(panel, -1, "&Save and exit")
            self.Bind(wx.EVT_BUTTON, self.saveExit, btn_SaveExit)

            #-- Create the close button, add it to the panel and wire it up to a function in the class
            btn_Close = wx.Button(panel, -1, "&Cancel")
            self.Bind(wx.EVT_BUTTON, self.onExit, btn_Close)

            #-- Now we have to create a grid to layout our controls
            sizer_main = wx.FlexGridSizer(rows=4,cols=1,hgap=1,vgap=5)

            filenames_string = ''
            for fn in self.caller.inFileNames:
                filenames_string += str(fn) + '\n'
            self.filenames_label = wx.TextCtrl(panel, -1, filenames_string,
              size=wx.Size(600,100), style=wx.TE_MULTILINE | wx.TE_DONTWRAP | wx.TE_READONLY )

            sizer_buttons = wx.FlexGridSizer(rows=1,cols=2, hgap=5,vgap=5)
            sizer_buttons.Add(btn_SaveExit)
            sizer_buttons.Add(btn_Close)

            #-- Add the grid to the panel and make it fit
            sizer_params = wx.FlexGridSizer(rows=3, cols=6, hgap=5, vgap=10)
            self.values = {}
            self.text_labels = {}

            params = {}
            params['twoth_zero_pixel'] = self.caller.twoth_zero_pixel
            params['twoth_offset'] = 0.0

            self.keys = [
            'twoth_zero_pixel',
            'twoth_offset',
            ]

            self.labels = [
              'Main beam center pixel at: ',
              'when A4 = '
              ]

            #self.q0_value = wx.TextCtrl(panel, 1, size=(100, -1))
            #self.q0_label = wx.StaticText(panel, -1, label)
            #sizer_params.Add(self.q0_label)
            #sizer_params.Add(self.q0_value)
            #self.q0_value.SetValue(
            #self.q0_value.V
            for key,label in zip(self.keys,self.labels):
                value = wx.TextCtrl(panel, 1, size=(100, -1))
                text_label = wx.StaticText(panel, -1, label)
                self.values[key] = value
                self.text_labels[key] = text_label
                sizer_params.Add(text_label)
                sizer_params.Add(value)
                value.SetValue(str(params[key]))

            #self.filenames = params['inFileNames']

            #-- Add the grid to the panel and make it fit

            sizer_main.Add(self.filenames_label)
            sizer_main.Add(btn_fileselector)
            sizer_main.Add(sizer_params)
            sizer_main.Add(sizer_buttons)
            panel.SetSizer(sizer_main)
            panel.Fit()
            #-- Show the window that we've just built

            self.ShowModal()

        def ChooseInputFile(self,e=None):
            #dlg = wx.FileDialog(self, "Choose a file", '', "", "I*.*", wx.FD_MULTIPLE)
            dlg = wx.FileDialog(None,
                                message="Choose Files",
                                defaultDir=os.getcwd(),
                                defaultFile="",
                                wildcard="I*.*",
                                style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR)
            if dlg.ShowModal() == wx.ID_OK:
                self.filenames=dlg.GetFilenames()
                self.dirname=dlg.GetDirectory()
                label = ''
                for fn in self.filenames:
                    label += str(self.dirname) + '/' + str(fn) + '\n'

                self.filenames_label.SetValue(label)

            dlg.Destroy()

        def saveExit(self, event):
            inFileNames = []
            for fn in self.filenames:
                inFileNames.append(self.dirname + '/' + fn)
            self.caller.inFileNames = inFileNames
            twoth_zero_pixel = float(self.values['twoth_zero_pixel'].GetValue())
            twoth_offset = float(self.values['twoth_offset'].GetValue())
            # now write value of q0_pixel including angle offset:
            self.caller.twoth_zero_pixel = twoth_zero_pixel - (self.caller.pixelsPerDegree * twoth_offset)
            self.Close(True)

        def onExit(self, event):
            self.Close(True)


    class get_params_gui(wx.Dialog):
        def __init__(self, parent, id, title, params):
            wx.Dialog.__init__(self, parent, wx.ID_ANY, title, size=(600,400), style = wx.DEFAULT_DIALOG_STYLE | wx.NO_FULL_REPAINT_ON_RESIZE)

            #-- I think you need a panel if you're going to have more than one control in your frame.
            panel = wx.Panel(self, -1)
            self.params = params

            ##-- create the file selection activator button
            #btn_fileselector = wx.Button(panel, -1, "Choose input &files: ")
            #self.Bind(wx.EVT_BUTTON, self.ChooseInputFile, btn_fileselector)

            #-- Create the processing button, add it to the panel and wire it up to a function in the class
            btn_SaveExit = wx.Button(panel, -1, "&Save and exit")
            self.Bind(wx.EVT_BUTTON, self.saveExit, btn_SaveExit)

            #-- Create the close button, add it to the panel and wire it up to a function in the class
            btn_Close = wx.Button(panel, -1, "&Cancel")
            self.Bind(wx.EVT_BUTTON, self.onExit, btn_Close)

            #-- Now we have to create a grid to layout our controls
            sizer_main = wx.FlexGridSizer(rows=4,cols=1,hgap=1,vgap=5)

            #filenames_string = ''
            #for fn in params['inFileNames']:
                #filenames_string += str(fn) + '\n'
            #self.filenames_label = wx.TextCtrl(panel, -1, filenames_string, size=wx.Size(600,100), style=wx.TE_MULTILINE | wx.TE_DONTWRAP | wx.TE_READONLY )

            sizer_buttons = wx.FlexGridSizer(rows=1,cols=2, hgap=5,vgap=5)
            sizer_buttons.Add(btn_SaveExit)
            sizer_buttons.Add(btn_Close)

            sizer_params = wx.FlexGridSizer(rows=3, cols=6, hgap=5, vgap=10)
            self.values = {}
            self.text_labels = {}

            self.keys = [
            'x_min',
            'x_max',
            'x_steps',
            'y_min',
            'y_max',
            'y_steps',
            'qZeroPixel',
            'qZeroAngle',
            ]

            self.labels = [
              'min Qx: ',
              'max Qx: ',
              'Qx steps: ',
              'min Qz: ',
              'max Qz: ',
              'Qz steps: ',
              'X pixel value for Q=0 : ',
              'when A4 = '
              ]
            for key,label in zip(self.keys,self.labels):
                value = wx.TextCtrl(panel, 1, size=(100, -1))
                text_label = wx.StaticText(panel, -1, label)
                self.values[key] = value
                self.text_labels[key] = text_label
                sizer_params.Add(text_label)
                sizer_params.Add(value)
                value.SetValue(str(params[key]))

            #self.filenames = params['inFileNames']

            #-- Add the grid to the panel and make it fit

            #sizer_main.Add(self.filenames_label)
            #sizer_main.Add(btn_fileselector)
            sizer_main.Add(sizer_params)
            sizer_main.Add(sizer_buttons)
            panel.SetSizer(sizer_main)
            panel.Fit()

            #-- Show the window that we've just built

            self.ShowModal()

        def saveExit(self, event):
            for key in self.keys:
                self.params[key] = float(self.values[key].GetValue())
            #self.params['inFileNames'] = self.filenames
            stored_name = '.stored_rebin_params'
            if os.access('./', os.W_OK):
                stored_file = open(stored_name, 'w')
                pickle.dump(self.params, stored_file)
                stored_file.close()
            self.Close(True)

        def onExit(self, event):
            self.Close(True)
