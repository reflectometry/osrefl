#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# check to see if all parameters are set on the command line
# if they are, then don't open GUI interface
from math import *
import sys
#from Tkinter import *
#import tkMessageBox
#import tkFileDialog
#from FileDialog import *
from osrefl.loaders.reduction import *
import osrefl.loaders.reduction as red
from numpy import *
from pylab import imshow,cm,colorbar,hot,show,xlabel,ylabel,connect, plot, figure, draw, axis, gcf
import matplotlib.colors as colors
from matplotlib.widgets import  RectangleSelector
from osrefl.viewers.colormap import change_colormap
from matplotlib.axis import XAxis, YAxis
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as Canvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg
import wx
import pickle
from matplotlib.image import FigureImage
from matplotlib.figure import Figure
from osrefl.viewers.zoom_colorbar import zoom_colorbar
from .reduction.cmapmenu import CMapMenu
from matplotlib.cm import get_cmap
import matplotlib.cbook as cbook
import matplotlib
from osrefl.viewers.plot_2d import plot_2d_data
from scipy import ndimage

def load_plottable_2d_data(filename=None):
    if filename == None:
        dlg = wx.FileDialog(None, "Load plottable_2d_data object from file:", '', "", "*.*", style=wx.FD_OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            fn = dlg.GetFilename()
            fd = dlg.GetDirectory()
        dlg.Destroy()
        filename = fd + '/' + fn
    initobj = pickle.load(filename)
    return initobj.instantiate()

def save_plottable_2d_data(obj, filename):
    obj.save_to_file(filename)

class p2d_initialization:
    """stores necessary data to recreate a plottable_2d_data object, for pickling"""
    def __init__(self, bin_data = None, params=None, creator=None, supervisor=None, plot_data=None, title=''):
        self.bin_data = bin_data
        self.params = params
        self.creator = creator
        self.supervisor = supervisor
        self.plot_data = plot_data
        self.title = title

    def instantiate(self):
        return plottable_2d_data(self.bin_data, self.params, self.creator, self.supervisor, self.plot_data, self.title)


class plottable_2d_data:
    """ container class, for holding an array of 2d data
    with normalization and pixel count data
    bin_data[:,:,0] = intensity
    bin_data[:,:,1] = pixels in bin
    bin_data[:,:,2] = monitor count
    bin_data[:,:,3] = average (intensity over total monitor count)

    also holds a dictionary of parameters with range boundaries and
    dimensions.
    Includes a plot method that shows a 2D log-intensity plot
    with correct axes
    Math methods rely on binning being the same for compared data
    """
    def __init__(self,
                 bin_data,
                 params_in={},
                 creator = None,
                 supervisor = None,
                 plot_data = True,
                 title = '',
                 base_data_obj = None,
                 **kwargs):
        #self.param_label = [ 'ymin', 'ymax', 'ysteps', 'xmin', 'xmax', 'xsteps' ]
        self.creator = creator # calling object
        self.params = params_in
        self.show_sliceplots = True #default = on
        self.bin_data = bin_data
        self.plot_data = plot_data
        # bin_data[:,:,0] = raw intensity (summed)
        # bin_data[:,:,1] = raw pixels in the bin
        # bin_data[:,:,2] = normalization (monitor or time, summed)
        # bin_data[:,:,3] = avg. intensity (raw intensity / normalization)
        self.data_array = []
        self.fig = None
        self.slice_y_data = None
        self.slice_x_data = None
        self.plot = None
        self.title = title
        self.supervisor = supervisor
        self.base_data_obj = base_data_obj

        if supervisor:
            self.register(supervisor)

        #self.xaxis_units = xaxis_units
        #self.yaxis_units = yaxis_units

    def __getstate__(self):
        return (self.bin_data, self.params, self.plot_data, self.title)

    def __setstate__(self, state):
        self.bin_data, self.params, self.plot_data, self.title = state


    def register(self, supervisor):
        supervisor.AddPlottable2dData(self, name = self.title, base_data_obj = self.base_data_obj)
        #self.number = supervisor.plottable_count
        self.supervisor = supervisor

    def __del__(self):
        if self.area_plot:
            self.area_plot.Close()

    def copy(self):
        new_data = plottable_2d_data(self.bin_data.copy(), self.params, supervisor = self.supervisor, title = self.title)
        return new_data

    def save_to_file(self, filename = None):
        if filename == None:
            dlg = wx.FileDialog(None, "Save plottable_2d_data object to file:", '', "", "", wx.FD_SAVE)
            if dlg.ShowModal() == wx.ID_OK:
                fn = dlg.GetFilename()
                fd = dlg.GetDirectory()
            dlg.Destroy()
            filename = fd + '/' + fn

        initobj = p2d_initialization(self.bin_data, self.params, self.creator, self.supervisor, self.plot_data, self.title)

        stored_file = open(filename, 'wb')
        pickle.dump(initobj, stored_file)

    def save2(self, filename = None):
        if filename == None:
            dlg = wx.FileDialog(None, "Save plottable_2d_data object to file:", '', "", "", wx.FD_SAVE)
            if dlg.ShowModal() == wx.ID_OK:
                fn = dlg.GetFilename()
                fd = dlg.GetDirectory()
            dlg.Destroy()
            filename = fd + '/' + fn

        stored_file = open(filename, 'wb')
        pickle.dump(self, stored_file)
        stored_file.close()

    def __sub__(self, other_data):
        if type(other_data) == type(self):
            return self.__subtract_otherdataset(other_data)
        elif isscalar(other_data):
            return self.__subtract_scalar(other_data)

    def __subtract_scalar(self, other_data):
        new_data = plottable_2d_data(zeros(self.bin_data.shape), self.params, supervisor = self.supervisor, title = '('+self.title+') - ('+str(other_data)+')')
        new_data.bin_data = self.bin_data.copy()
        new_data.bin_data[:,:,3] -= other_data
        return new_data

    def __subtract_otherdataset(self, other_data):
        if other_data.bin_data.shape == self.bin_data.shape:
            if self.check_compatible_binning(other_data):
                new_data = plottable_2d_data(zeros(self.bin_data.shape), self.params, supervisor = self.supervisor, title = '('+self.title+') - ('+other_data.title+')')
                self_nonzero = (self.bin_data[:,:,1] != 0)
                other_nonzero = (other_data.bin_data[:,:,1] != 0)
                mask = self_nonzero * other_nonzero
                print len(mask)
                new_data.bin_data[:,:,3][mask] = self.bin_data[:,:,3][mask] - other_data.bin_data[:,:,3][mask]
                new_data.bin_data[:,:,2][mask] = self.bin_data[:,:,2][mask] + other_data.bin_data[:,:,2][mask]
                new_data.bin_data[:,:,1][mask] = self.bin_data[:,:,1][mask] + other_data.bin_data[:,:,1][mask]
                new_data.bin_data[:,:,0][mask] = self.bin_data[:,:,0][mask] - other_data.bin_data[:,:,0][mask]
                return new_data
            else:
                print("error: datasets are not identically binned")
        elif other_data.bin_data.shape[1] == 1:
            new_data = plottable_2d_data(zeros(self.bin_data.shape), self.params, supervisor = self.supervisor, title = '('+self.title+') - ('+other_data.title+')')
            new_data.bin_data[:,:,0] = self.bin_data[:,:,0] - other_data.bin_data[:,:,0]
            new_data.bin_data[:,:,1] = self.bin_data[:,:,1] + other_data.bin_data[:,:,1]
            new_data.bin_data[:,:,2] = self.bin_data[:,:,2] + other_data.bin_data[:,:,2]
            new_data.bin_data[:,:,3] = self.bin_data[:,:,3] - other_data.bin_data[:,:,3]
            return new_data
        else:
            print("subtraction failed for some reason")

    def check_compatible_binning(self, other_data):
        compatible =  ( other_data.__class__.__name__ == self.__class__.__name__ )
        compatible &= ( other_data.params['y_steps'] == self.params['y_steps'] )
        compatible &= ( other_data.params['x_steps'] == self.params['x_steps'] )
        compatible &= ( other_data.params['x_min'] == self.params['x_min'] )
        compatible &= ( other_data.params['y_min'] == self.params['y_min'] )
        compatible &= ( other_data.params['x_max'] == self.params['x_max'] )
        compatible &= ( other_data.params['y_max'] == self.params['y_max'] )
        ### compatibility shouldn't depend on unit labels for axes, but uncomment if you want it to.
        #compatible &= ( other_data.params['x_units'] == self.params['x_units'] )
        #compatible &= ( other_data.params['y_units'] == self.params['y_units'] )

        return compatible

    def __add__(self, other_data):
        if type(other_data) == type(self):
            print 'adding two data sets \n'
            return self.__add_otherdataset(other_data)
        elif isscalar(other_data):
            print 'adding scalar to data set\n'
            return self.__add_scalar(other_data)

    def __add_scalar(self, other_data):
        new_data = self.copy()
        new_data.bin_data[:,:,3] += other_data
        return new_data

    def __add_otherdataset(self, other_data):
        if self.check_compatible_binning(other_data):
            new_data = plottable_2d_data(self.bin_data, self.params, supervisor = self.supervisor, title = '('+self.title+') + ('+other_data.title+')')
            new_data.bin_data = self.bin_data + other_data.bin_data
            #new_data.bin_data = self.bin_data.copy()
            #for j in range(self.params['y_steps']):
                #for i in range(self.params['x_steps']):
            #new_data.bin_data[i,j,0] = self.bin_data[i,j,0] + other_data.bin_data[i,j,0]
            #new_data.bin_data[i,j,1] = self.bin_data[i,j,1] + other_data.bin_data[i,j,1]
            #new_data.bin_data[i,j,2] = self.bin_data[i,j,2] + other_data.bin_data[i,j,2]
            #new_data.bin_data[i,j,3] = self.bin_data[i,j,3] + other_data.bin_data[i,j,3]
            return new_data
        else:
            print("error: datasets are not identically binned")
            return

    def __mul__(self, multiplier):
        multiplier = float(multiplier)
        new_data = plottable_2d_data(self.bin_data, self.params, supervisor = self.supervisor, title = '('+self.title+') * ' + str(multiplier))
        new_data.bin_data = self.bin_data.copy()
        for j in range(self.params['y_steps']):
            for i in range(self.params['x_steps']):
                new_data.bin_data[i,j,0] = self.bin_data[i,j,0] * multiplier
                new_data.bin_data[i,j,1] = self.bin_data[i,j,1]
                new_data.bin_data[i,j,2] = self.bin_data[i,j,2]
                new_pixelCount = new_data.bin_data[i,j,1]
                if ( new_pixelCount > 0 ):
                    new_monitorTotal = new_data.bin_data[i,j,2]
                    new_avg = new_data.bin_data[i,j,0] / double(new_monitorTotal)
                else:
                    new_avg = 0.0
                new_data.bin_data[i,j,3] = new_avg

        return new_data

    #def __div__(self, other_data):
        #"""divide one dataset by another - useful for
        #dividing by a background, for instance"""
        #if not self.check_compatible_binning(other_data):
            #print("error: datasets are not identically binned")
            #return self
        #if not other_data.bin_data[self.bin_data[:,:,0] != 0].min() > 0:
        ## checking for denominator positive for all data points
            #print("error: attempting to divide by zero")
            #return self
        #else:
            #new_data = plottable_2d_data(self.bin_data, self.params)
            #new_data.bin_data = self.bin_data.copy()
            #for j in range(self.params['y_steps']):
                #for i in range(self.params['x_steps']):

                    #new_data.bin_data[i,j,0] = self.bin_data[i,j,0] / other_data.bin_data[i,j,0]
                    #new_data.bin_data[i,j,1] = self.bin_data[i,j,1]
                    #new_data.bin_data[i,j,2] = self.bin_data[i,j,2]
                    #new_pixelCount = new_data.bin_data[i,j,1]
                    #if ( new_pixelCount > 0 ):
                #new_monitorTotal = new_data.bin_data[i,j,2]
                #new_avg = new_data.bin_data[i,j,0] / double(new_monitorTotal)
                    #else:
                #new_avg = 0.0
                    #new_data.bin_data[i,j,3] = new_avg

            #return new_data

    def __div__(self, other_data):
        """divide one dataset by another - useful for
        dividing by a background, for instance"""
        if not self.check_compatible_binning(other_data):
            print("error: datasets are not identically binned")
            return self
        if not other_data.bin_data[:,:,3][self.bin_data[:,:,1] != 0].min() > 0:
        # checking for denominator positive for all data points
            print("error: attempting to divide by zero")
            return self
        else:
            new_data = self.copy()
            new_data.bin_data *= 0
            self_nonzero = (self.bin_data[:,:,1] != 0)
            other_nonzero = (other_data.bin_data[:,:,1] != 0)
            mask = self_nonzero * other_nonzero
            new_data.bin_data[:,:,3][mask] = self.bin_data[:,:,3][mask] / other_data.bin_data[:,:,3][mask]
            new_data.bin_data[:,:,1][mask] = self.bin_data[:,:,1][mask] + other_data.bin_data[:,:,1][mask]
            new_data.bin_data[:,:,2][mask] = self.bin_data[:,:,2][mask] + other_data.bin_data[:,:,2][mask]
            new_data.bin_data[:,:,0][mask] = self.bin_data[:,:,0][mask] + other_data.bin_data[:,:,0][mask]
            return new_data

    def overlap(self, other_data):
        new_data = self.copy()
        mask = (other_data.bin_data[:,:,1] == 0)
        new_data.bin_data[:,:,0][mask] = 0
        new_data.bin_data[:,:,1][mask] = 0
        new_data.bin_data[:,:,2][mask] = 0
        new_data.bin_data[:,:,3][mask] = 0

        return new_data

    def smooth(self, smoothing_width = 3.0, axis = 1):
        working_copy = self.bin_data.copy()
        nz_mask = working_copy[:,:,1].nonzero()
        working_copy[:,:,0] = ndimage.gaussian_filter1d(working_copy[:,:,0], smoothing_width, axis=axis)
        working_copy[:,:,3][nz_mask] = working_copy[:,:,0][nz_mask] / working_copy[:,:,2][nz_mask]
        self.bin_data = working_copy
        if self.area_plot:
            self.area_plot.Destroy()
            self.area_plot = self.wxplot()

    def save(self, outFileName = None):
        if outFileName == None:
            dlg = wx.FileDialog(None, "Save 2d data as:", '', "", "", wx.FD_SAVE)
            if dlg.ShowModal() == wx.ID_OK:
                fn = dlg.GetFilename()
                fd = dlg.GetDirectory()
            dlg.Destroy()
            outFileName = fd + '/' + fn

        outFile = open(outFileName, 'w')
        # write header: parameter list, then column headers, then data.
        outFile.write("#" + str(self.params) + "\n")
        outFile.write("x\ty\tI/monitor\tI_0\n")
        y_max = self.params['y_max']
        y_min = self.params['y_min']
        y_steps = self.params['y_steps']
        y_units = self.params['y_units']
        x_max = self.params['x_max']
        x_min = self.params['x_min']
        x_steps = self.params['x_steps']
        x_units = self.params['x_units']
        xStepSize = ( y_max - y_min ) / y_steps
        yStepSize = ( x_max - x_min ) / x_steps
        yArray = arange(y_steps, dtype='float') / y_steps * (y_max - y_min) + y_min
        xArray = arange(x_steps, dtype='float') / x_steps * (x_max - x_min) + x_min
        for i in range(x_steps):
            for j in range(y_steps):
                x = xArray[i]
                y = yArray[j]
                I_0 = self.bin_data[i,j,0]
                avg = self.bin_data[i,j,3]
                outFile.write(str(x) + "\t" + str(y) + "\t")
                if ( avg > 0.0 ):
                    outFile.write(str(avg) + "\t" + str(I_0))
                else:
                    outFile.write('-0\t-0')
                outFile.write('\n')

        outFile.close()
        return

    def onselect(self, eclick, erelease):
        x_range = [eclick.xdata, erelease.xdata]
        y_range = [eclick.ydata, erelease.ydata]
        ax = eclick.inaxes
        self.sliceplot(x_range, y_range, ax=ax)
        print 'sliceplot([%f,%f],[%f,%f])' % (x_range[0],x_range[1],y_range[0],y_range[1])

    def sliceplot(self, x_range, y_range, ax = None):
        """sum along x and z within the box defined by qX- and qZrange.
        sum along qx is plotted to the right of the data,
        sum along qz is plotted below the data.
        Transparent white rectangle is overlaid on data to show summing region"""

        x, slice_y_data, y, slice_x_data = self.do_xy_slice(x_range, y_range)
        self.x = x
        self.slice_y_data = slice_y_data
        self.y = y
        self.slice_x_data = slice_x_data
        self.slice_xrange = x_range
        self.slice_yrange = y_range

        if self.area_plot:
            self.area_plot.show_slice_overlay(x_range, y_range, x, slice_y_data, y, slice_x_data)

    def do_xy_slice(self, x_range, y_range):
        """ slice up the data, once along x and once along z.
        returns 4 arrays:  a y-axis for the x data,
        an x-axis for the y data."""
        params = self.params
        print 'doing xy slice'
        data = self.bin_data[:,:,3].copy()
        pixels = self.bin_data[:,:,1]
        # zero out any pixels in the sum that have zero in the pixel count:
        data[pixels == 0] = 0

        normalization_matrix = ones(pixels.shape)
        normalization_matrix[pixels == 0] = 0
        x_min = min(x_range)
        x_max = max(x_range)
        y_min = min(y_range)
        y_max = max(y_range)

        x_size,y_size = data.shape
        global_x_range = (params['x_max'] - params['x_min'])
        global_y_range = (params['y_max'] - params['y_min'])

        x_pixel_min = round( (x_min - params['x_min']) / global_x_range * x_size )
        x_pixel_max = round( (x_max - params['x_min']) / global_x_range * x_size )
        y_pixel_min = round( (y_min - params['y_min']) / global_y_range * y_size )
        y_pixel_max = round( (y_max - params['y_min']) / global_y_range * y_size )

        #correct any sign switches:
        if (x_pixel_min > x_pixel_max):
            new_min = x_pixel_max
            x_pixel_max = x_pixel_min
            x_pixel_min = new_min

        if (y_pixel_min > y_pixel_max):
            new_min = y_pixel_max
            y_pixel_max = y_pixel_min
            y_pixel_min = new_min

        new_x_min = x_pixel_min / x_size * global_x_range + params['x_min']
        new_x_max = x_pixel_max / x_size * global_x_range + params['x_min']
        new_y_min = y_pixel_min / y_size * global_y_range + params['y_min']
        new_y_max = y_pixel_max / y_size * global_y_range + params['y_min']

        x_pixel_min = int(x_pixel_min)
        x_pixel_max = int(x_pixel_max)
        y_pixel_min = int(y_pixel_min)
        y_pixel_max = int(y_pixel_max)

        y_norm_factor = sum(normalization_matrix[x_pixel_min:x_pixel_max,y_pixel_min:y_pixel_max], axis=1)
        x_norm_factor = sum(normalization_matrix[x_pixel_min:x_pixel_max,y_pixel_min:y_pixel_max], axis=0)
        # make sure the normalization has a minimum value of 1 everywhere,
        # to avoid divide by zero errors:
        y_norm_factor[y_norm_factor == 0] = 1
        x_norm_factor[x_norm_factor == 0] = 1

        slice_y_data = sum(data[x_pixel_min:x_pixel_max,y_pixel_min:y_pixel_max], axis=1) / y_norm_factor
        slice_x_data = sum(data[x_pixel_min:x_pixel_max,y_pixel_min:y_pixel_max], axis=0) / x_norm_factor

        #slice_y_data = slice_y_data
        #slice_x_data = slice_x_data

        x_vals = arange(slice_y_data.shape[0], dtype = 'float') / slice_y_data.shape[0] * (new_x_max - new_x_min) + new_x_min
        y_vals = arange(slice_x_data.shape[0], dtype = 'float') / slice_x_data.shape[0] * (new_y_max - new_y_min) + new_y_min

        return x_vals, slice_y_data, y_vals, slice_x_data



    def log_lin_select(self,event):
        if not (isinstance(event.artist, XAxis) or isinstance(event.artist, YAxis)):
            return
        ax = event.artist.axes
        label = ax.get_label()
        if label == 'sz':
            scale = ax.get_yscale()
            if scale == 'log':
                ax.set_yscale('linear')
                ax.figure.canvas.draw()
            elif scale == 'linear':
                ax.set_yscale('log')
                ax.figure.canvas.draw()

        elif label == 'sx':
            scale = ax.get_xscale()
            if scale == 'log':
                ax.set_xscale('linear')
                ax.figure.canvas.draw()
            elif scale == 'linear':
                ax.set_xscale('log')
                ax.figure.canvas.draw()
        return

    def toggle_selector(self, event):
        print ' Key pressed.'
        if event.key in ['C', 'c']:
            print 'change colormap.'
            ax = event.inaxes
            change_colormap(ax.images[0])
        if event.key in ['Q', 'q'] and self.RS.active:
            print ' RectangleSelector deactivated.'
            self.RS.set_active(False)
        if event.key in ['A', 'a'] and not self.RS.active:
            print ' RectangleSelector activated.'
            self.RS.set_active(True)


    def save_slice(self, outFileName, header = ""):
        outFile = open(outFileName, 'w')
        outFile.write(header)
        if not (self.slice_qx_data == None):
            for i in range(self.slice_qx_data.shape[0]):
                x = self.qz_vals[i]
                y = self.slice_qx_data[i]
                outFile.write(str(x) + "\t" + str(y) + "\n")
        outFile.close()
        print('saved qx slice in %s' % (outFileName))
        return

    def save_x_slice(self, event=None, outFileName=None):
        if outFileName == None:
            dlg = wx.FileDialog(None, "Save 2d data as:", '', "", "", wx.FD_SAVE)
            if dlg.ShowModal() == wx.ID_OK:
                fn = dlg.GetFilename()
                fd = dlg.GetDirectory()
            dlg.Destroy()
            outFileName = fd + '/' + fn
        outFile = open(outFileName, 'w')
        outFile.write('#'+self.title+'\n')
        outFile.write('#xmin: ' + str(self.slice_xrange[0]) + '\n')
        outFile.write('#xmax: ' + str(self.slice_xrange[1]) + '\n')
        outFile.write('#ymin: ' + str(self.slice_yrange[0]) + '\n')
        outFile.write('#ymax: ' + str(self.slice_yrange[1]) + '\n')
        outFile.write("#y\tslice_x_data\n")
        if not (self.slice_x_data == None):
            for i in range(self.slice_x_data.shape[0]):
                x = self.y[i]
                y = self.slice_x_data[i]
                outFile.write(str(x) + "\t" + str(y) + "\n")
        outFile.close()
        print('saved x slice in %s' % (outFileName))
        return

    def save_y_slice(self, event=None, outFileName=None):
        if outFileName == None:
            dlg = wx.FileDialog(None, "Save 2d data as:", '', "", "", wx.FD_SAVE)
            if dlg.ShowModal() == wx.ID_OK:
                fn = dlg.GetFilename()
                fd = dlg.GetDirectory()
            dlg.Destroy()
            outFileName = fd + '/' + fn
        outFile = open(outFileName, 'w')
        outFile.write('#'+self.title+'\n')
        outFile.write('#xmin: ' + str(self.slice_xrange[0]) + '\n')
        outFile.write('#xmax: ' + str(self.slice_xrange[1]) + '\n')
        outFile.write('#ymin: ' + str(self.slice_yrange[0]) + '\n')
        outFile.write('#ymax: ' + str(self.slice_yrange[1]) + '\n')
        outFile.write("#x\tslice_y_data\n")
        if not (self.slice_y_data == None):
            for i in range(self.slice_y_data.shape[0]):
                x = self.x[i]
                y = self.slice_y_data[i]
                outFile.write(str(x) + "\t" + str(y) + "\n")
        outFile.close()
        print('saved y slice in %s' % (outFileName))
        return

    def plot_qx_slice(self, label = '', figNum = None):
        if figNum:
            fig = figure(figNum)
        else:
            fig = figure()
            figNum = fig.number
        plot(self.qz_vals, self.slice_qx_data, label = label)
        return

    def plot_qz_slice(self, label = '', figNum = None):
        if figNum:
            fig = figure(figNum)
        else:
            fig = figure()
            figNum = fig.number
        plot(self.qx_vals, self.slice_qz_data, label = label)
        return

    def logplot(self, show_sliceplots=True):
        from zoom_colorbar import zoom_colorbar
        x_min = self.params['x_min']
        x_max = self.params['x_max']
        y_min = self.params['y_min']
        y_max = self.params['y_max']
        self.show_data = zeros((self.params['x_steps'],self.params['y_steps']))

        self.minimum_intensity = inf
        for j in range(self.params['y_steps']):
            for i in range(self.params['x_steps']):
                avg = self.bin_data[i,j,3]
                if avg > 0.0:
                    self.show_data[i,j] = avg
                else:
                    self.show_data[i,j] = 0.0
                if (avg < self.minimum_intensity and avg > 0):
                    self.minimum_intensity = avg

        #self.show_data = transpose(log(self.show_data + self.minimum_intensity / 2.0))

        fig = figure()
        self.fig = fig
        connect('pick_event', self.log_lin_select)

        if show_sliceplots:
            ax = fig.add_subplot(221, label='qxqz_plot')
            fig.sx = fig.add_subplot(222, label='sx', picker=True)
            fig.sx.xaxis.set_picker(True)
            fig.sx.yaxis.set_picker(True)
            fig.sz = fig.add_subplot(223, label='sz', picker=True)
            fig.sz.xaxis.set_picker(True)
            fig.sz.yaxis.set_picker(True)
            self.RS = RectangleSelector(ax, self.onselect, drawtype='box', useblit=True)
            fig.slice_overlay = None
        else:
            ax = fig.add_subplot(111, label='qxqz_plot')
        fig.ax = ax
        ax.set_title(self.params['description'])
        connect('key_press_event', self.toggle_selector)
        transformed_show_data = transpose(log(self.show_data + self.minimum_intensity / 2.0))
        im = ax.imshow(transformed_show_data, interpolation='nearest', aspect='auto', origin='lower',cmap=cm.jet, extent=(x_min,x_max,y_min,y_max))
        fig.im = im
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(self.ylabel)
        zoom_colorbar(im)
        figure(fig.number)
        fig.canvas.draw()
        return im

    def plot(self, show_sliceplots=True):
        from zoom_colorbar import zoom_colorbar
        x_min = self.params['x_min']
        x_max = self.params['x_max']
        y_min = self.params['y_min']
        y_max = self.params['y_max']
        self.show_data = zeros((self.params['x_steps'],self.params['y_steps']))

        self.minimum_intensity = inf
        for j in range(self.params['y_steps']):
            for i in range(self.params['x_steps']):
                avg = self.bin_data[i,j,3]
                self.show_data[i,j] = avg
                if (avg < self.minimum_intensity and avg > 0):
                    self.minimum_intensity = avg

        fig = figure()
        self.fig = fig
        connect('pick_event', self.log_lin_select)

        if show_sliceplots:
            ax = fig.add_subplot(221, label='qxqz_plot')
            fig.sx = fig.add_subplot(222, label='sx', picker=True)
            fig.sx.xaxis.set_picker(True)
            fig.sx.yaxis.set_picker(True)
            fig.sz = fig.add_subplot(223, label='sz', picker=True)
            fig.sz.xaxis.set_picker(True)
            fig.sz.yaxis.set_picker(True)
            self.RS = RectangleSelector(ax, self.onselect, drawtype='box', useblit=True)
            fig.slice_overlay = None
        else:
            ax = fig.add_subplot(111, label='qxqz_plot')
        fig.ax = ax
        ax.set_title(self.params['description'])
        connect('key_press_event', self.toggle_selector)
        transformed_show_data = transpose(self.show_data)
        im = ax.imshow(transformed_show_data, interpolation='nearest', aspect='auto', origin='lower',cmap=cm.hot, extent=(x_min,x_max,y_min,y_max))
        fig.im = im
        ax.set_xlabel('Qx (inv Angstroms)')
        ax.set_ylabel('Qz (inv Angstroms)')
        #fig.colorbar(im, ax=ax)
        zoom_colorbar(im)
        #im.set_cmap(cm.hot)
        #fig.show()
        figure(fig.number)
        draw()
        return im

    def wxplot(self, destroy_older = True, scale = 'log'):
        #use the custom WX Frame defined below, with custom toolbar including slice button

        x_min = self.params['x_min']
        x_max = self.params['x_max']
        y_min = self.params['y_min']
        y_max = self.params['y_max']
        self.show_data = self.bin_data[:,:,3].copy()
        pixel_mask = self.bin_data[:,:,1].copy()

        extent=[x_min,x_max,y_min,y_max]
        #from plot_2d3 import plot_2d_data


        plot_title = self.params['description']
        x_label = self.params['x_units']
        y_label = self.params['y_units']
        frame = offspec_plot_2d_data(self.show_data, extent, self, scale = scale, pixel_mask = pixel_mask, window_title = self.title, plot_title = plot_title, x_label = x_label, y_label = y_label)
        frame.Show()
        self.area_plot = frame
        return frame



class offspec_plot_2d_data(plot_2d_data):
    """overriding the context menus to add interaction with other objects known to supervisor"""
    def get_all_plot_2d_instances(self):
        """get all other plots that are open (from supervisor?)"""
        if self.caller == None:
            return [], []

        supervisor = self.caller.supervisor
        instances = []
        instance_names = []
        #for dataset in supervisor.rebinned_data_objects:
            ##instances.append(dataset)
            #for subkey in dataset.__dict__.keys():
                    #if isinstance(dataset.__dict__[subkey], plottable_2d_data):
                ##print('plottable_2d_data yes')
                #instance_names.append(str(dataset.number) + ': ' + subkey + ': ' + dataset.description)
                #instances.append(dataset.__dict__[subkey])
        for plottable in supervisor.plottable_2d_data_objects:
            if hasattr(plottable, 'area_plot'):
                instances.append(plottable.area_plot)
                instance_names.append(plottable.title + ': ' + plottable.params['description'])

        return instances, instance_names

    def other_plots_menu(self):
        other_plots, other_plot_names = self.get_all_plot_2d_instances()
        other_menu = wx.Menu()
        for op in other_plot_names:
            item = other_menu.Append(wx.ID_ANY, op, op)
        return other_menu

    def other_plots_dialog(self):
        other_plots, other_plot_names = self.get_all_plot_2d_instances()
        #selection_num = wx.GetSingleChoiceIndex('Choose other plot', '', other_plot_names)
        dlg = wx.SingleChoiceDialog(None, 'Choose other plot', '', other_plot_names)
        dlg.SetSize(wx.Size(640,480))
        if dlg.ShowModal() == wx.ID_OK:
            selection_num=dlg.GetSelection()
        dlg.Destroy()
        return other_plots[selection_num]

    def dummy(self, evt):
        print 'the event is: ' + str(evt)

    def area_context(self, mpl_mouseevent, evt):
        area_popup = wx.Menu()
        item1 = area_popup.Append(wx.ID_ANY,'&Grid on/off', 'Toggle grid lines')
        wx.EVT_MENU(self, item1.GetId(), self.OnGridToggle)
        cmapmenu = CMapMenu(self, callback = self.OnColormap, mapper=self.mapper, canvas=self.canvas)
        item2 = area_popup.Append(wx.ID_ANY,'&Toggle log/lin', 'Toggle log/linear scale')
        wx.EVT_MENU(self, item2.GetId(), lambda evt: self.toggle_log_lin(mpl_mouseevent))
        item3 = area_popup.AppendMenu(wx.ID_ANY, "Colourmaps", cmapmenu)
        #other_plots, other_plot_names = self.get_all_plot_2d_instances()
        #if not (other_plot_names == []):
            #other_menu = wx.Menu()
            #for op in other_plot_names:
                #item = other_menu.Append(wx.ID_ANY, op, op)
        #other_menu = self.other_plots_menu()
        item4 = area_popup.Append(wx.ID_ANY, "copy intens. scale from", '')
        wx.EVT_MENU(self, item4.GetId(), lambda evt: self.copy_intensity_range_from(self.other_plots_dialog()) )
        item5 = area_popup.Append(wx.ID_ANY, "copy slice region from", '')
        wx.EVT_MENU(self, item5.GetId(), lambda evt: self.sliceplot(self.other_plots_dialog().slice_xy_range) )
        self.PopupMenu(area_popup, evt.GetPositionTuple())
