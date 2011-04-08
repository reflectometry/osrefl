#!/usr/bin/env python
from numpy import *

import matplotlib
from matplotlib.widgets import RectangleSelector
from matplotlib.blocking_input import BlockingInput
from ginput_rect import *

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

# Wx-Pylab magic ...
from matplotlib import _pylab_helpers
from matplotlib.backend_bases import FigureManagerBase

import wx

def fit_run(data, model,plot_extent):
    print 'Model', model[5,5]
    app = App(data=data, model=model,plot_extent=plot_extent)
    model = fliplr((app.frame.panel.scaled_data).T)
    
    app.MainLoop()
    print 'scaled data', app.frame.panel.scaled_data[5,5]
    return model

class App(wx.App):
    def __init__(self,data, model, plot_extent):
        wx.App.__init__(self)

        self.frame = NormFrame(self,data,model,plot_extent)
        self.frame.Show()
        self.SetTopWindow(self.frame)
        

    
class NormFrame(wx.Frame):
    def __init__(self,parent,data,model,plot_extent):
        wx.Frame.__init__(self,None,-1,"Normalize")
        self.data = data
        self.parent = parent
        self.model = model
        self.plot_extent = plot_extent
        self.panel = scalePanel(self,-1)
  
  
class scalePanel(wx.Panel):
        def __init__(self,parent,id):
            wx.Panel.__init__(self,parent,id,style=wx.BORDER_SUNKEN)
            
            self.scale_collection = array([])
            self.real_data = flipud(parent.data.T)
            self.parent = parent
            
            self.raw_data = flipud(parent.model.T)
            self.scaled_data = copy(self.raw_data)
            self.plot_extent = parent.plot_extent
            
            lower_lim = amin(self.scaled_data[nonzero(self.scaled_data.real)])
            
            finite_real = self.real_data[isfinite(self.real_data)]
            finite_real = finite_real[nonzero(finite_real)]
            
            #Hack
            #self.vmin = amin(log(finite_real.real))
            self.vmax = amax(log(finite_real.real))
            self.vmin = self.vmax - 12
            self.figure = Figure()
            self.axes = self.figure.add_subplot(211)
            self.canvas = FigureCanvas(self, -1, self.figure)
            
            fm = FigureManagerBase(self.canvas, 0)
            _pylab_helpers.Gcf.set_active(fm)

            # Instantiate the matplotlib navigation toolbar and explicitly show it.
            mpl_toolbar = Toolbar(self.canvas)
            mpl_toolbar.Realize()
            
            self.myImage = self.showImage(log(abs(self.scaled_data)),self.axes)
            self.modelColor = self.figure.colorbar(self.myImage)
            
            self.dataaxes = self.figure.add_subplot(212)
            self.datamyImage = self.showImage(log(self.real_data),self.dataaxes)
            self.dataColor = self.figure.colorbar(self.datamyImage)
            
            self.scale = wx.TextCtrl(self,-1)
            
            self.updateButton = wx.Button(self,-1,'UPDATE')
            self.resetButton = wx.Button(self,-1,'RESET')
            self.areaScaleButton = wx.Button(self,-1,'AUTO SCALE')

            
            BotSize = wx.BoxSizer(wx.HORIZONTAL)
            vertSize = wx.BoxSizer(wx.VERTICAL)
            
            BotSize.Add(self.scale,0,wx.LEFT|wx.RIGHT,border = 5)
            BotSize.Add(self.updateButton,0,wx.LEFT|wx.RIGHT,border = 5)
            BotSize.Add(self.resetButton,0,wx.LEFT|wx.RIGHT,border = 5)
            BotSize.Add(self.areaScaleButton,0,wx.LEFT|wx.RIGHT, border = 5)
            
            vertSize.Add(self.canvas,-1,wx.EXPAND)
            vertSize.Add(mpl_toolbar,0)
            vertSize.Add(BotSize,0,wx.ALL, border = 10)
            self.SetSizer(vertSize)
            self.Fit()
            
            self.updateButton.Bind(wx.EVT_BUTTON,self.update)
            self.resetButton.Bind(wx.EVT_BUTTON,self.reset)
            self.areaScaleButton.Bind(wx.EVT_BUTTON,self.autoScale)
            
    
        def plotPanel(self):
    
            a = array(arange(20))
            b = array(arange(30))
         
            a = a.reshape(20,1)
            
            b = b.reshape(1,30)
            
            self.real_data = a*b
            self.scaled_data = copy(self.real_data)
          
            return 
          
            
        def reset(self,event):
            self.scaled_data = copy(self.raw_data)
            lower_lim = amin(self.scaled_data[nonzero(self.scaled_data.real)])
            self.myImage = self.showImage(log(abs(self.scaled_data+lower_lim)),self.axes)
            
            self.canvas.draw()
            
            self.scale_collection = []
            
            
            
        def update(self,event):
            self.scale_collection = hstack((float(self.scale.GetValue()),self.scale_collection))
            print sum(self.scale_collection)
            
            self.scaled_data = self.raw_data * (self.scale_collection)
            print self.scaled_data[10]
            #lower_lim = amin(self.scaled_data[nonzero(self.scaled_data.real)])
           
            self.myImage = self.showImage(log(abs(self.scaled_data)),self.axes)
            self.canvas.draw()
        
        def autoScale(self,event):
            
            x1,x2,z1,z2 = ginput_rect(self.datamyImage)
            
        def showImage(self,data,axis):
            return axis.imshow(data, aspect = 'auto',extent = self.plot_extent,vmin = self.vmin, vmax = self.vmax)
        
    
if __name__ == '__main__':

    a = array(arange(20))
    b = array(arange(30))

    a = a.reshape(20,1)
    
    b = b.reshape(1,30)
    
    raw_data = a*b
    scaled_data = copy(raw_data)
    plot_extent = array([0,20,0,20])
    print 'start'
    c = fit_run(raw_data,scaled_data,plot_extent)
    print 'continued'