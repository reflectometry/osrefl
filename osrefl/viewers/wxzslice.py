# Copyright (C)
# All rights reserved.
# See LICENSE.txt for details.
# Author: Brian Maranville, Christopher Metting

#Starting Date:8/23/2010

from numpy import arctan2,indices,array,ma,pi,amin,amax,nan,degrees, isfinite

import matplotlib,wx
from matplotlib.widgets import RectangleSelector
from matplotlib.blocking_input import BlockingInput
from matplotlib.colors import LinearSegmentedColormap
from pylab import figure,show,imshow,draw, delaxes,cla

# Disable interactive mode so that plots are only updated on show() or draw().
# Note that the interactive function must be called before selecting a backend
# or importing pyplot, otherwise it will have no effect.
matplotlib.interactive(False)

# Specify the backend to use for plotting and import backend dependent classes.
# Note that this must be done before importing pyplot to have an effect.
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2Wx as Toolbar

# The Figure object is used to create backend-independent plot representations.
from matplotlib.figure import Figure

from matplotlib import pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from pylab import get_current_fig_manager as gcfm
# Wx-Pylab magic ...
from matplotlib import _pylab_helpers
from matplotlib.backend_bases import FigureManagerBase


def momentView(omf):

    app = wx.PySimpleApp()
    app.Frame = MomentFrame(-1,omf)
    app.Frame.Show(True)
    app.MainLoop()
    app.Destroy()
    return()


class MomentFrame(wx.Frame):
    '''

    '''
    def __init__(self,parent,omf):
        wx.Frame.__init__(self,None,-1,"Z Slice Viewer")
        self.omf = omf
        self.parent = parent
        self.panel = ZMomentplot(self,-1)
        self.Fit()


class ZMomentplot(wx.Panel):
    def __init__(self,frame,id):
        wx.Panel.__init__(self,frame,id,style = wx.BORDER_RAISED)

        z_layer = 0
        self.omf = frame.omf
        self.halfstep = (float(self.omf.parameters['zstepsize'])/2.0)*1.0e10

        self.figure = Figure(frameon = True)
        self.figure.set_facecolor('.82')

        self.canvas = FigureCanvas(self, -1, self.figure)

        fm = FigureManagerBase(self.canvas, 0)
        _pylab_helpers.Gcf.set_active(fm)

        self.wheel_cmap = LinearSegmentedColormap.from_list('wheel_rgby',
                               ['red', 'green', 'blue', 'yellow', 'red'])

        mpl_toolbar = Toolbar(self.canvas)
        mpl_toolbar.Realize()

        x,y = indices((100,100), dtype = float) - 50.
        wheel_angle = arctan2(x,y)

        self.ax1 = self.figure.add_axes([0.1,0.1,0.7,0.8])

        self.angle = ma.array(arctan2(-self.omf.my[:,:,z_layer],
                                 self.omf.mx[:,:,z_layer]),
                                  mask=(self.omf.M[:,:,z_layer] == 0.0))


        xmax = (float(self.omf.parameters['xnodes']) *
                float(self.omf.parameters['xstepsize']) * 1.0e10)

        ymax = (float(self.omf.parameters['ynodes']) *
                float(self.omf.parameters['ystepsize']) * 1.0e10)

        self.extent = [0., xmax, 0., ymax]

        self.im = self.ax1.imshow(self.angle.T, origin='lower',
                                  interpolation = 'nearest',
                                  extent = self.extent,
                                  cmap = self.wheel_cmap,
                                  vmin = -pi, vmax = pi)


        self.ax1.set_title('Z Slice = ' +
               str(z_layer*float(self.omf.parameters['zstepsize'])* 1.0e10
                   + self.halfstep) +' Ang' ,size = 'xx-large')

        self.ax1.set_xlabel('$x (\AA)$', size = 'x-large')
        self.ax1.set_ylabel('$y (\AA)$', size = 'x-large')


        self.ax1w = self.figure.add_axes([0.75,0.4,0.3,0.2], polar=True)

        self.ax1w.yaxis.set_visible(False)

        self.ax1w.imshow(wheel_angle, cmap=self.wheel_cmap,
                         extent=[0,2*pi,0,pi])

        self.ax1w.set_title('M direction\n(in-plane)')


        self.zselect = wx.Slider(self,-1,size = [300,40],minValue = int(0),
                                 maxValue = int(self.omf.dims[2]-1),
                                 style = wx.SL_AUTOTICKS)


        self.datavalue = wx.StatusBar(self,-1)
        self.datavalue.SetStatusText('Angle Value: ')
        self.label = wx.StaticText(self,-1,'Select Z layer: ')

        self.minVal = wx.StaticText(self,-1, '0.0')

        self.maxVal = wx.StaticText(self,-1, str(self.omf.dims[2]*
                         (float(self.omf.parameters['zstepsize'])* 1.0e10)
                         -self.halfstep))

        #Sizer Creation
        toolSize = wx.BoxSizer(wx.HORIZONTAL)
        BotSize = wx.BoxSizer(wx.HORIZONTAL)
        vertSize = wx.BoxSizer(wx.VERTICAL)

        BotSize.Add(self.label,0,wx.LEFT|wx.RIGHT,border = 5)
        BotSize.Add(self.minVal,0,wx.LEFT|wx.RIGHT,border = 5)
        BotSize.Add(self.zselect,0,wx.TOP|wx.LEFT|wx.RIGHT,border = 5)
        BotSize.Add(self.maxVal,0,wx.LEFT|wx.RIGHT,border = 5)

        toolSize.Add(mpl_toolbar,0,wx.RIGHT,border = 20)
        toolSize.Add(self.datavalue,0,wx.LEFT|wx.RIGHT|wx.TOP,border = 20)

        vertSize.Add(self.canvas,-1,wx.EXPAND|wx.LEFT|wx.RIGHT,border = 5)
        vertSize.Add(toolSize,0)
        vertSize.Add(BotSize,0,wx.ALL, border = 10)

        self.SetSizer(vertSize)

        self.Fit()


        self.zselect.Bind(wx.EVT_SCROLL,self.newMomentZ)
        self.canvas.mpl_connect('motion_notify_event',self.onMouseOver)

    def newMomentZ(self,event):

        z_layer = self.zselect.GetValue()
        Zvalue = str(float(self.omf.parameters['zstepsize'])*z_layer)


        self.ax1.set_title('Moment at Z = '+ str(Zvalue))

        self.angle = ma.array(arctan2(-self.omf.my[:,:,z_layer],
                                 self.omf.mx[:,:,z_layer]),
                                 mask=(self.omf.M[:,:,z_layer] == 0))

        self.im.set_data(self.angle.T)

        self.ax1.set_title('Z Slice = ' +
               str(z_layer*float(self.omf.parameters['zstepsize'])* 1.0e10
                   + self.halfstep) + ' Ang',size = 'xx-large')

        draw()

    def onMouseOver(self, event):
        """

        """

        if event.inaxes == self.ax1:
            if (event.xdata != None and event.ydata != None):

                xidx = (int(event.xdata/
                        (float(self.omf.parameters['xstepsize'])*1.0e10)))

                yidx = (int(event.ydata/
                        (float(self.omf.parameters['ystepsize'])*1.0e10)))

                value = self.angle[xidx,yidx]


                if ma.getmask(self.angle)[xidx,yidx]:
                    self.datavalue.SetStatusText('Angle Value: MASKED')

                else:
                    value = -degrees(value)
                    if (value < 0.0): value +=360
                    self.datavalue.SetStatusText('Angle Value: '
                                                 +str('%.2f'%value))
        else:

            self.datavalue.SetLabel('')

        return



def unitView(v,extent = None, step = None, n = None):
    '''
    **Overview:**

        This method is used to plot an array with a slicing option in the z
        direction. The plotter contains a slider bar which allows the user to
        scroll through the different layers and view the x-y slices. This method
        is more reliable and less expensive  than the mlab.contour3d method
        which is still included in the software functionality.

    **Parameters:**

        *v:* (float:3D array|angstroms^-2)
            This is the 3D array which will be plotted. It is generally used for
            viewing the unit cell SLD object, however may be extended to other
            uses.

    '''
    from numpy import shape,array
    dim = shape(array(v))
    print dim
    if extent == None:
        extent = array([[0,dim[0]],[0,dim[1]],[0,dim[2]]])
    if step == None:
        step = [1,1,1]
    if n == None:
        n = dim

    app = wx.PySimpleApp()
    app.Frame = unitFrame(-1,v,extent,step,n)
    app.Frame.Show(True)
    app.MainLoop()
    app.Destroy()
    return()


class unitFrame(wx.Frame):
    '''

    '''
    def __init__(self,parent,v,extent,step,n):
        wx.Frame.__init__(self,None,-1,"Z Slice Viewer")
        self.v = v
        self.extent = extent
        self.step = step
        self.n = n

        self.parent = parent
        self.panel = ZUnitPlot(self,-1)
        self.Fit()

class ZUnitPlot(wx.Panel):
    def __init__(self,frame,id):
        wx.Panel.__init__(self,frame,id,style = wx.BORDER_RAISED)

        self.v = frame.v
        self.extent = frame.extent

        self.step = frame.step
        self.n = frame.n

        z_layer = 0
        self.halfstep = self.step[2]/2.0

        self.figure = Figure(frameon = True)

        self.figure.set_facecolor('.82')
        self.canvas = FigureCanvas(self, -1, self.figure)

        fm = FigureManagerBase(self.canvas, 0)
        _pylab_helpers.Gcf.set_active(fm)
        self.ax1 = self.figure.add_axes([0.1,0.1,0.7,0.8])

        plotExtent = ([self.extent[0,0],self.extent[0,1],
                       self.extent[1,0],self.extent[1,1]])

        self.im = self.ax1.imshow(self.v[:,:,z_layer].T,
                                  origin='lower',
                                  interpolation = 'nearest',
                                  extent = plotExtent,
                                  vmin = amin(self.v),
                                  vmax = amax(self.v))

        self.ax1.set_title('Z Slice = ' +
                        str(z_layer *self.step[2] + self.halfstep) +
                        ' Ang' ,size = 'xx-large')

        self.figure.colorbar(self.im,format = '%.2e')

        mpl_toolbar = Toolbar(self.canvas)
        mpl_toolbar.Realize()

        self.zselect = wx.Slider(self,-1,size = [300,40],minValue = int(0),
                                 maxValue = int(self.n[2]-1),
                                 style = wx.SL_AUTOTICKS)

        print self.extent[2,1]
        print self.halfstep
        print 'TEST', str(self.extent[2,1] - self.halfstep)
        self.label = wx.StaticText(self,-1,'Select Z layer: ')

        self.minVal = wx.StaticText(self,-1, '0.0')

        self.maxVal = wx.StaticText(self,-1, str(self.extent[2,1]
                                                 - self.halfstep))


        BotSize = wx.BoxSizer(wx.HORIZONTAL)
        vertSize = wx.BoxSizer(wx.VERTICAL)

        BotSize.Add(self.label,0,wx.LEFT|wx.RIGHT,border = 5)
        BotSize.Add(self.minVal,0,wx.LEFT|wx.RIGHT,border = 5)
        BotSize.Add(self.zselect,0,wx.TOP|wx.LEFT|wx.RIGHT,border = 5)
        BotSize.Add(self.maxVal,0,wx.LEFT|wx.RIGHT,border = 5)

        vertSize.Add(self.canvas,-1,wx.EXPAND|wx.LEFT|wx.RIGHT)
        vertSize.Add(mpl_toolbar,0)
        vertSize.Add(BotSize,0,wx.ALL, border = 10)

        self.SetSizer(vertSize)

        self.Fit()
        self.zselect.Bind(wx.EVT_SCROLL,self.newUnitZ)

    def newUnitZ(self,event):

        z_layer = self.zselect.GetValue()
        Zvalue = str(float(self.step[2])*z_layer)


        self.ax1.set_title('Moment at Z = '+ str(Zvalue))

        newData = self.v[:,:,z_layer]
        self.im.set_data(newData.T)

        self.ax1.set_title('Z Slice = ' +
                        str(z_layer *self.step[2] + self.halfstep) +
                        ' Ang' ,size = 'xx-large')

        draw()

def _test():


    from omfLoader import *
    from numpy import asarray

    mag = Omf('/home/mettingc/Documents/BBM/c_500mTz_-Oxs_MinDriver-Magnetization-05-0005651.omf')
    #mag = Omf('/home/mettingc/Documents/test.omf')
    momentView(mag)


    import sample_prep
    Au = (sample_prep.Parallelapiped(SLD = 4.506842e-6,
                                     Ms = 8.6e5,dim=[5.0e4,5.0e4,2.0e4]))

    Cr = (sample_prep.Layer(SLD = 3.01e-6,Ms = 8.6e5,
                            thickness_value = 1000.0))

    #Au.on_top_of(Cr)
    scene = sample_prep.Scene([Au])

    GeoUnit = (sample_prep.GeomUnit(Dxyz = [10.0e4,10.0e4,2.2e4],
                                    n = [20,25,20],scene = scene))

    cell = GeoUnit.buildUnit()
    print cell.step
    unitView(cell.unit,asarray([[0,cell.Dxyz[0]],[0,cell.Dxyz[1]],[0,cell.Dxyz[2]]]),cell.step,cell.n)

if __name__=="__main__":_test()
