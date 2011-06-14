# Copyright (C) 2008 University of Maryland
# All rights reserved.
# See LICENSE.txt for details.
# Author: Christopher Metting

#Starting Date:8/23/2010


from numpy import size, shape, array, log, zeros, ones, argwhere, asarray, ma
from numpy import isneginf,log10,inf,nonzero
from numpy import amin, amax, sum, ceil,arange,linspace, exp
import matplotlib,wx
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.widgets import RectangleSelector
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from mpl_toolkits.axes_grid import ImageGrid,Grid
from matplotlib.widgets import RectangleSelector
from matplotlib.blocking_input import BlockingInput
from matplotlib.colors import LinearSegmentedColormap
from pylab import figure,show,imshow, draw, delaxes,cla,subplot,setp,plot
from pylab import figlegend,legend,xlim
from copy import copy
import sys
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
from pylab import get_current_fig_manager as gcfm

# Wx-Pylab magic ...
from matplotlib import _pylab_helpers
from matplotlib.backend_bases import FigureManagerBase


def MultiView(data, step, n, extent=None, titles = None,
              axisLabel = None, vlimit = None):

    app = wx.PySimpleApp()

    infoCol = PlotInfo(data,vlimit,titles,extent,step,n,axisLabel)
    controller = PlotCtrler(app,infoCol)
    app.MainLoop()
    app.Destroy()

    return


class MultiViewFrame(wx.Frame):
    '''
    This is the main Frame used to collect the Panels
    '''
    def __init__(self,parent, pos = (20,20), size = (1920,1280)):

        wx.Frame.__init__(self,parent,title = "Plot Comparison Application",pos=pos,size=size)

        return
    

class ParentVisPanel(wx.Panel):
    '''
    This is the parent panel which is owned by the frame and owns all of the sub
    panels.
    '''
    def __init__(self,parent):
        wx.Panel.__init__(self,parent,style = wx.BORDER_RAISED)
        #self.SetBackgroundColour('purple')

class plotNotebook(wx.Notebook):
    '''
    The GUI can hold multiple plots to view. This notebook has a page for each
    set of data plus the first page which is a thumbnail view of all of the
    plots.
    '''
    def __init__(self,parent):
        wx.Notebook.__init__(self, parent,style = wx.BORDER_RAISED)
        #self.SetBackgroundColour('blue')

class MultiViewPanel(wx.Panel):
    '''
    This is top panel which displays all of the 2D plots in the notebook it is
    somewhat rigid because it was written prior to the rest of the GUI.
    '''
    def __init__(self,parent,id, data, extent, titles, axLabel, vlimit,
                 selMask,scale = 'log'):

        wx.Panel.__init__(self,parent,id,style = wx.BORDER_SUNKEN)

        FIGPROPS = dict(figsize = (5.0,(5.0/1.618)),dpi = 64)
        #self.SetBackgroundColour('red')
        ADJUSTPROPS = dict(left = 0.1, bottom = 0.1, right = 0.97,
                           top = 0.94, wspace = 0.01, hspace = 0.01 )

        viewIdx = argwhere(selMask).flatten()
        viewCount = size(viewIdx)

        self.plotFigure = Figure(frameon = True,**FIGPROPS)
        self.plotFigure.set_facecolor('.82')
        self.plotCanvas = FigureCanvas(self, -1, self.plotFigure)
        self.plotFigure.subplots_adjust(**ADJUSTPROPS)

        fm = FigureManagerBase(self.plotCanvas, 0)
        _pylab_helpers.Gcf.set_active(fm)



        if (scale == 'linear'):
            plotData = data
            vmin = vlimit[0]
            vmax = vlimit[1]
        elif (scale == 'log'):
            plotData = log10(data)
            vmin = log10(vlimit[0])
            vmax = log10(vlimit[1])
            for i in range(viewCount):
                plotData[viewIdx[i]][data[viewIdx[i]] == 0.0] = vmin

        depth = ceil(float(viewCount/2.0))
        self.axes = [None]*int(viewCount)
        self.im = [None]*int(viewCount)

        if (viewCount == 1):
            col = 1
            row = 1
        elif (viewCount%3.0 == 0.0):
            col = 3
            row =(int(viewCount/3.0))
        elif (viewCount%2.0 == 0.0):
            col = 2
            row = (int(viewCount/2.0))
        else:
            print 'Plot Layout Could Not Be Built'
        for i in range(viewCount):

            if i == 0:
                self.axes[i] = self.plotFigure.add_subplot(row,col,i+1)

            else:
                self.axes[i] = self.plotFigure.add_subplot(row,col,i+1,
                                                       sharex = self.axes[0],
                                                       sharey = self.axes[0])

            self.im[i] = self.axes[i].imshow(plotData[viewIdx[i]].T,
                                             origin='lower',
                                     interpolation = 'nearest',
                                     extent = extent,aspect = 'auto',
                                     vmin = vmin, vmax = vmax)

            self.axes[i].set_gid(viewIdx[i])

            if viewCount > 1:
                self.axes[i].text(0.01,.92,titles[viewIdx[i]],fontsize=14,
                                  bbox = dict(facecolor = 'white',alpha = 0.5),
                                  transform = self.axes[i].transAxes)

                self.plotFigure.text(0.54,0.01,axLabel[0],
                                     ha = 'center',va = 'bottom',fontsize=14)

                self.plotFigure.text(0.01,.5,axLabel[1],
                                     ha = 'left',va = 'center',
                                     rotation = 'vertical',fontsize=14)
            else:
                self.axes[i].set_title(titles[viewIdx[i]])
                self.axes[i].set_xlabel(axLabel[0],fontsize=14)
                self.axes[i].set_ylabel(axLabel[1],fontsize=14)
                self.plotFigure.colorbar(self.im[i])

            if i+1 <= viewCount-col:
                setp(self.axes[i].get_xticklabels(), visible=False)


            if ((i+1)%(viewCount/row) == 0.0 and viewCount != 1):

                setp(self.axes[i].get_yticklabels(), visible=False)

        self.box = Toolbar(self.plotCanvas)
        self.box.Hide()

        multiViewSizer = wx.BoxSizer(wx.HORIZONTAL)
        multiViewVert = wx.BoxSizer(wx.VERTICAL)

        multiViewSizer.Add(self.plotCanvas,1,wx.EXPAND|wx.RIGHT|wx.LEFT,
                                                            border = 1)
        multiViewVert.Add(multiViewSizer,1,wx.EXPAND|wx.TOP|wx.BOTTOM,
                                                            border = 1)


        self.SetSizer(multiViewVert)


        return

class MouseReadPan(wx.Panel):
    '''
    A panel for reading out all three dimensions of a 2D plot. This is important
    for getting intensity values for scaleing. It works by first focusing on the
    panel for which the mouse is currently over and then reading the values off
    of the plot.
    '''
    def __init__(self,parent):
        wx.Panel.__init__(self,parent,style = wx.BORDER_RAISED)
        self.readoutLabel = wx.StaticText(self,-1)
        self.readoutLabel.SetLabel('    X              Y             Z ')
        self.readoutData = wx.StaticText(self,-1)
        newFont = self.readoutData.GetFont()
        newFont.SetPointSize(8)
        self.readoutData.SetFont(newFont)
        self.readoutLabel.SetFont(newFont)
        self.readoutData.SetLabel('  ')

        self.readSizer = wx.BoxSizer(wx.VERTICAL)
        self.readSizer.Add(self.readoutLabel)
        self.readSizer.Add(self.readoutData)
        self.SetSizer(self.readSizer)
        self.Fit()

class ScalePan(wx.Panel):
    '''
    This is the panel which allows the user to toggle the 1D plot between a log
    and linear scale.
    '''
    def __init__(self,parent):
        wx.Panel.__init__(self,parent,style = wx.BORDER_RAISED)

        choices = ['Linear','Semi-log']
        self.scaleBox = wx.RadioBox(self,1,'Y-Axis Scale',
                                      choices = choices)

        self.scaleSizer = wx.BoxSizer(wx.VERTICAL)
        self.scaleSizer.Add(self.scaleBox)
        self.Fit()

        return


class ZlimPan(wx.Panel):
    '''
    This panel allows the user to change the floor and ceiling limits on the
    data. The order of magnitude the theory calculates is much greater then what
    the instrument realistically has the capability to measure so, by adjusting
    the ceiling and floor to match the data, we can get a better handle on how
    well the theory matches.
    '''
    def __init__(self,parent):
        wx.Panel.__init__(self,parent,style = wx.BORDER_RAISED)
        wx.StaticText
        self.vminLab = wx.StaticText(self,-1,'   v min')
        self.vmaxLab = wx.StaticText(self,-1,'   v max')


        newFont = self.vminLab.GetFont()
        newFont.SetPointSize(8)

        self.vminLab.SetFont(newFont)
        self.vmaxLab.SetFont(newFont)

        self.vminBox = wx.TextCtrl(self, 1,size = (10, 20))
        self.vmaxBox = wx.TextCtrl(self, 2,size = (10, 20))
        self.submit = wx.Button(self,3, 'Submit',(10, 10))
        self.reset = wx.Button(self,4, 'Reset',(10, 10))
        
        self.vminBox.SetFont(newFont)
        self.vmaxBox.SetFont(newFont)
        
        self.butSize = wx.BoxSizer(wx.HORIZONTAL)
        self.panSize = wx.BoxSizer(wx.VERTICAL)
        self.labVal = wx.BoxSizer(wx.HORIZONTAL)
        self.enterVal = wx.BoxSizer(wx.HORIZONTAL)

        self.butSize.Add(self.submit)
        self.butSize.Add(self.reset)

        self.labVal.Add(self.vminLab,1,wx.EXPAND|wx.RIGHT|wx.LEFT|wx.TOP,
                                                                    border = 1)
        self.labVal.Add(self.vmaxLab,1,wx.EXPAND|wx.RIGHT|wx.LEFT,border = 1)

        self.enterVal.Add(self.vminBox,1,wx.RIGHT|wx.LEFT,border = 1)
        self.enterVal.Add(self.vmaxBox,1,wx.RIGHT|wx.LEFT,border = 1)

        self.panSize.Add(self.labVal,1,wx.EXPAND|wx.RIGHT|wx.LEFT,border = 1)
        self.panSize.Add(self.enterVal,1,wx.EXPAND|wx.RIGHT|wx.LEFT,border = 1)
        self.panSize.Add(self.butSize,1,wx.EXPAND|wx.ALL,border = 1)
        self.SetSizer(self.panSize)

        self.Fit()
        return

class HozVertPan(wx.Panel):
    '''
    Some of the other slicers out there use a 1D plot box on the right side for
    the vertical slice and a 1D plot on the bottom for the horizontal. In
    reality, the user is most likely only going to care about one slice at a
    timeand the vertical 1D plot is somewhat clunky so this button allows the
    user to toggle the plot axis for the single 1D plot at the bottom of the
    GUI.
    '''
    def __init__(self,parent):
        wx.Panel.__init__(self,parent,style = wx.BORDER_RAISED)

        choices = ['Vertical','Horizontal']
        self.HozVertBox = wx.RadioBox(self,1,'Slice Orientation',
                                      choices = choices)

        self.radioSizer = wx.BoxSizer(wx.VERTICAL)
        self.radioSizer.Add(self.HozVertBox)
        self.Fit()

        return

class slicePanel(wx.Panel):
    '''
    This is the 1D plot panel which displays the slices which the user selects
    '''

    def __init__(self,parent):
        wx.Panel.__init__(self,parent,style = wx.BORDER_SUNKEN)

        self.sliceFigure = Figure(frameon = True,figsize =(2,3))
        self.sliceFigure.set_facecolor('.82')
        self.sliceCanvas = FigureCanvas(self, -1, self.sliceFigure)
        self.sliceCanvas.SetAutoLayout(False)

        self.box = Toolbar(self.sliceCanvas)
        self.box.Hide()
        self.box.zoom()

        self.fm = FigureManagerBase(self.sliceCanvas, 1)
        _pylab_helpers.Gcf.set_active(self.fm)

        self.sliceAxis = self.sliceFigure.add_axes([0.1,0.2,0.8,0.7])

        sliceSizer = wx.BoxSizer(wx.VERTICAL)
        sliceSizer.Add(self.sliceCanvas,1,wx.EXPAND|wx.RIGHT|wx.LEFT,
                                                                    border = 1)
        self.SetSizer(sliceSizer)
        return

    def updatePlot(self,data,xaxis,legLab,xlabel,ylabel,title,scale,legTog):
        '''
        The user can select a new slice at any time to plot on the 1D panel.
        This method refreshes the plot so that the new data is displayed.
        '''
        self.sliceAxis.cla()

        self.linePlot = self.sliceAxis.plot(xaxis,(array(data).T))
        xlim((min(xaxis),max(xaxis)))

        self.sliceAxis.leg = legend(legLab)

        self.sliceAxis.leg.set_visible(legTog)

        self.sliceAxis.set_xlabel(xlabel)
        self.sliceAxis.set_ylabel(ylabel)
        self.sliceAxis.set_title(title)

        self.sliceAxis.set_yscale(str(scale))

        draw()

        return


class dataScale(wx.Panel):
    '''
    This panel allows the user to scale the theory to match it with the data.
    This is done through 'markers' which tell the panel whether or not the
    array being plotted is 'Theory' or 'Measured'.
    '''
    def __init__(self,parent):
        wx.Panel.__init__(self,parent,style = wx.BORDER_RAISED)
        self.submit = wx.Button(self,13, 'Submit',(10, 10))
        self.reset = wx.Button(self,12, 'Reset',(10, 10))
        self.scaleFactor = wx.TextCtrl(self, 11,size = (10, 10))
        self.scaleFactorLabel = wx.StaticText(self,1,'Theory Scaling Factor')
        newFont = self.scaleFactorLabel.GetFont()
        newFont.SetPointSize(8)

        self.scaleFactorLabel.SetFont(newFont)

        self.horScaleSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.horSize = wx.BoxSizer(wx.HORIZONTAL)
        self.vertSize = wx.BoxSizer(wx.VERTICAL)

        self.horScaleSizer.Add(self.scaleFactor,1,wx.EXPAND,border = 1)

        self.horSize.Add(self.submit,1,wx.EXPAND|wx.RIGHT|wx.LEFT|wx.TOP,
                                                                    border = 1)
        self.horSize.Add(self.reset,1,wx.EXPAND|wx.RIGHT|wx.LEFT,border = 1)

        self.vertSize.Add(self.scaleFactorLabel,1,wx.EXPAND|wx.RIGHT|wx.LEFT,
                                                                    border = 1)
        self.vertSize.Add(self.horScaleSizer,1,wx.EXPAND|wx.RIGHT|wx.LEFT,
                                                                    border = 1)
        self.vertSize.Add(self.horSize,1,wx.EXPAND|wx.RIGHT|wx.LEFT,border = 1)

        self.SetSizer(self.vertSize)
        return


class PlotCtrler(object):
    '''
    ***********CONTROLLER********************
    '''
    def __init__(self, app,plotInfo):
        self.plotInfo = plotInfo
        
        pos,size = self.window_placement(1920,1280)
        self.Frame = MultiViewFrame(None,pos,size)

        menuBar = wx.MenuBar()
        file = wx.Menu()
        options = wx.Menu()

        quit = wx.MenuItem(file,1,'&Quit\tCtrl+Q')
        save2D = wx.MenuItem(file,2,'&Save2D Img\tCtrl+S')
        export2D = wx.MenuItem(file,3,'Export 2D Data\ tCtrl+E')

        file.AppendItem(quit)
        file.AppendItem(save2D)
        file.AppendItem(export2D)

        logTog = wx.MenuItem(options, 11, '&Linear Scale\tCtrl+l',
                             'Changes Image Plot to Log Scale', wx.ITEM_CHECK)
        legTog = wx.MenuItem(options,12,'&Legend \tl',
                             'Toggles the Legend in the plot', wx.ITEM_CHECK)

        options.AppendItem(logTog)
        options.AppendItem(legTog)

        options.Check(12,True)

        menuBar.Append(file,'&File')
        menuBar.Append(options,'&Options')

        self.Frame.Bind(wx.EVT_MENU,self.imScaleChange,id = 11)
        self.Frame.Bind(wx.EVT_MENU,self.toggleLegend,id = 12)
        self.Frame.Bind(wx.EVT_MENU,self.onQuit,id = 1)
        self.Frame.Bind(wx.EVT_MENU,self.sliceSave,id = 2)
        self.Frame.Bind(wx.EVT_MENU, self.export2D, id = 3)
        self.Frame.SetMenuBar(menuBar)

        self.parPanel = ParentVisPanel(self.Frame)


        self.notePan = plotNotebook(self.parPanel)
        self.initNotebook()

        self.getCurrPage(None)

        self.slicePan = slicePanel(self.parPanel)
        self.orientation = HozVertPan(self.parPanel)
        self.scale = ScalePan(self.parPanel)
        self.readOut = MouseReadPan(self.parPanel)

        self.zlim = ZlimPan(self.parPanel)

        self.dataMatch = dataScale(self.parPanel)

        self.plotContrSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.toolSizer = wx.BoxSizer(wx.VERTICAL)

        self.toolSizer.Add(self.readOut,1,wx.EXPAND|wx.LEFT|wx.RIGHT,
                           border = 1)
        self.toolSizer.Add(self.orientation,1,wx.EXPAND|wx.LEFT|wx.RIGHT,
                           border = 1)
        self.toolSizer.Add(self.scale,1,wx.EXPAND|wx.LEFT|wx.RIGHT,
                           border = 1)
        self.toolSizer.Add(self.zlim,1,wx.EXPAND|wx.LEFT|wx.RIGHT,
                           border = 1)
        self.toolSizer.Add(self.dataMatch,1,wx.EXPAND|wx.LEFT|wx.RIGHT,
                           border = 1)


        self.plotContrSizer.Add(self.slicePan,1,wx.EXPAND|wx.LEFT|wx.RIGHT,
                                border = 1)
        self.plotContrSizer.Add(self.toolSizer,0,wx.EXPAND|wx.LEFT|wx.RIGHT,
                                border = 1)

        self.vertSizer = wx.BoxSizer(wx.VERTICAL)

        self.vertSizer.Add(self.notePan,2,wx.EXPAND|wx.LEFT|wx.RIGHT,border = 1)
        self.vertSizer.Add(self.plotContrSizer,1,wx.EXPAND|wx.LEFT|wx.RIGHT,
                           border = 1)

        self.parPanel.SetSizer(self.vertSizer)

        self.slicePan.sliceCanvas.Bind(wx.EVT_RIGHT_DOWN,self.sliceHome)
        self.parPanel.Fit()

        self.Frame.Show(True)

        self.notePan.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, self.OnPageChanged)

        self.orientation.Bind(wx.EVT_RADIOBOX,self.toggleOrientation)
        self.scale.Bind(wx.EVT_RADIOBOX,self.toggleScale)

        self.zlim.Bind(wx.EVT_BUTTON,self.setZ,self.zlim.submit)
        self.zlim.Bind(wx.EVT_BUTTON,self.resetZ,self.zlim.reset)


        self.dataMatch.Bind(wx.EVT_BUTTON, self.setDataScale,
                            self.dataMatch.submit)
        self.dataMatch.Bind(wx.EVT_BUTTON, self.resetDataScale,
                            self.dataMatch.reset)

        return
    
    def window_placement(self, desired_width, desired_height):
        """
        Determines the position and size of a window such that it fits on the
        user's screen without obstructing (or being obstructed by) the task bar.
        The returned size is bounded by the desired width and height passed in,
        but it may be smaller if the screen is too small.  Usually the returned
        position (upper left coordinates) will result in centering the window
        on the screen excluding the task bar area.  However, for very large
        monitors it will be placed on the left side of the screen.
        """

        # WORKAROUND: When running Linux and using an Xming (X11) server on a
        # PC with a dual monitor configuration, the reported display count may
        # be 1 (instead of 2) with a display size of both monitors combined.
        # (For example, on a target PC with an extended desktop consisting of
        # two 1280x1024 monitors, the reported monitor size was 2560x1045.)
        # To avoid displaying the window across both monitors, we check for
        # screen 'too big'.  If so, we assume a smaller width which means the
        # application will be placed towards the left hand side of the screen.
        
        x, y, w, h = wx.Display().GetClientArea() # size excludes task bar
        #print "*** x, y, w, h", x, y, w, h
        xpos, ypos = x, y
        h -= 20  # to make room for Mac window decorations
        if len(sys.argv) > 1 and '--platform' in sys.argv[1:]:
            j, k = wx.DisplaySize()  # size includes task bar area
            print "*** Reported screen size including taskbar is %d x %d"%(j, k)
            print "*** Reported screen size excluding taskbar is %d x %d"%(w, h)

        if w > 1920: w = 1280  # display on left side, not centered on screen
        if w > desired_width:  xpos = x + (w - desired_width)/2
        if h > desired_height: ypos = y + (h - desired_height)/2

        # Return the suggested position and size for the application frame.
        return (xpos, ypos), (min(w, desired_width), min(h, desired_height))


    def OnPageChanged(self, event):
        self.getCurrPage(event)
        event.Skip()

    def clearText(self):
        self.zlim.vminBox.Clear()
        self.zlim.vmaxBox.Clear()
        self.dataMatch.scaleFactor.Clear()

        return

    def resetDataScale(self,event):

        self.plotInfo.data = copy(self.plotInfo.preservedData)
        self.plotInfo.initVlim()
        curPage = self.notePan.GetSelection()
        self.imScaleChange(event)
        self.notePan.SetSelection(curPage)

        self.clearText()
        return

    def setDataScale(self,event):

        if self.dataMatch.scaleFactor.GetValue() != '':
            self.plotInfo.data[self.plotInfo.type == 'Theory'] = (
                  self.plotInfo.preservedData[self.plotInfo.type == 'Theory']*
                  float(self.dataMatch.scaleFactor.GetValue()))
            self.plotInfo.initVlim()
            curPage = self.notePan.GetSelection()
            self.imScaleChange(event)
            self.notePan.SetSelection(curPage)
        self.clearText()
        return


    def resetZ(self,event):
        self.plotInfo.initVlim()
        curPage = self.notePan.GetSelection()
        self.imScaleChange(event)
        self.notePan.SetSelection(curPage)
        self.clearText()
        return

    def setZ(self,event):

        vminHold = self.plotInfo.vlimit[0]
        vmaxHold = self.plotInfo.vlimit[1]

        if self.zlim.vminBox.GetValue()!= '':
            self.plotInfo.vlimit[0] = float(self.zlim.vminBox.GetValue())
            if self.plotInfo.imScale == 'log':
                self.plotInfo.vlimit[0] = 10**(self.plotInfo.vlimit[0])

        if self.zlim.vmaxBox.GetValue() !='':
            self.plotInfo.vlimit[1] = float(self.zlim.vmaxBox.GetValue())
            if self.plotInfo.imScale == 'log':
                self.plotInfo.vlimit[1] = 10**(self.plotInfo.vlimit[1])

        if self.plotInfo.vlimit[0] >= self.plotInfo.vlimit[1]:
            print 'ERROR: Limits incorrectly chosen'
            self.plotInfo.vlimit[0] = vminHold
            self.plotInfo.vlimit[1] = vmaxHold

        curPage = self.notePan.GetSelection()
        self.imScaleChange(event)
        self.notePan.SetSelection(curPage)
        self.clearText()
        return

    def imScaleChange(self,event):
        if event.GetSelection() == 0:
            self.plotInfo.imScale = 'log'
        else:
            self.plotInfo.imScale = 'linear'
        self.initNotebook()
        _pylab_helpers.Gcf.set_active(self.slicePan.fm)
        return

    def onQuit(self,event):
        self.Frame.Destroy()
        self.Frame.Close()

    def getCoord(self,event):

        if (event.xdata == None or event.xdata == None):
            self.readOut.readoutData.SetLabel('  ')

            return
        else:
            x = event.xdata
            y = event.ydata
            xcorrd = self._convertCoord(self.plotInfo.extent[0],
                                           self.plotInfo.step[0],[x])[0]
            ycorrd = self._convertCoord(self.plotInfo.extent[1],
                                           self.plotInfo.step[1],[y])[0]

            if (xcorrd >= (self.plotInfo.n[0]) or ycorrd >= (
                                                         self.plotInfo.n[1])):
                self.readOut.readoutData.SetLabel(' ')

            else:
                z = self.plotInfo.data[int(event.inaxes.get_gid())][xcorrd,ycorrd]
                self.readOut.readoutData.SetLabel('%.2e  '%x +'%.2e  '%y +'%.2e  '%z)

        return

    def getCurrPage(self,event):
        if event is not None:
            self.currPage = self.notePan.GetPage(event.GetSelection())
        else:
            self.currPage = self.notePan.GetCurrentPage()

        #for coordinates
        self.currPage.plotCanvas.mpl_connect('motion_notify_event',
                                             self.getCoord)

        #for plot
        self.coordBeg = self.currPage.plotCanvas.mpl_connect(
                                         'button_press_event',self.getSliceBeg)

        self.coordEnd = self.currPage.plotCanvas.mpl_connect(
                                     'button_release_event',self.getSliceEnd)

    def getPos(self,event):
        self.curX = event.xdata
        self.curY = event.ydata
        return

    def export2D(self,event):
        print 'NOT IMPLEMENTED'
        self.plotInfo.sliceData
        self.sliceX

    def sliceSave(self,event):
        self.slicePan.box.save(event)
        return
    def sliceHome(self,event):
        self.slicePan.box.home()
        return

    def getSliceBeg(self,event):

        if (event.xdata == None or event.xdata == None):
            self.plotInfo.selecting = False
            return
        else:
            self.plotInfo.selecting = True
            self.start = [event.x,event.y]
            self.valStart = [event.xdata,event.ydata]
            self.lastLocalx = event.x
            self.lastLocaly = event.y
            self.lastValuex = event.x
            self.lastValuey = event.y

            self.probeAx = event.inaxes

            self.moveEvt = self.currPage.plotCanvas.mpl_connect(
                                        'motion_notify_event',self.boxUpdate)
        return

    def _convertCoord(self,min,step,value):
        coord = zeros(size(value))

        for i in range(size(value)):
            coord[i] = round((value[i]+abs(min))/step)
        return sorted(coord)

    def getSliceEnd(self,event):
        self.currPage.box.draw_rubberband(self,0.0,0.0,0.0,0.0)
        #self.currPage.plotCanvas.mpl_disconnect(self.moveEvt)
        if hasattr(self, 'moveEvt'):
            self.currPage.plotCanvas.mpl_disconnect(self.moveEvt)

        if self.plotInfo.selecting == True:

            self.plotInfo.selecting = False

            self.currPage.plotCanvas.mpl_disconnect(self.coordBeg)
            self.currPage.plotCanvas.mpl_disconnect(self.coordEnd)

            self.currPage.plotCanvas.blit(self.probeAx.bbox)

            xSelect = sorted([self.valStart[0],self.lastValuex])
            ySelect = sorted([self.valStart[1],self.lastValuey])

            xArea = self._convertCoord(self.plotInfo.extent[0],
                                       self.plotInfo.step[0],xSelect)
            yArea = self._convertCoord(self.plotInfo.extent[1],
                                       self.plotInfo.step[1],ySelect)

            #A crude way to test if the area selected is to small
            if ((xArea[1] - xArea[0] > self.plotInfo.n[0]) or
                (yArea[1] - yArea[0] > self.plotInfo.n[1])):
                return
            else:
                self.plotInfo.axisSelect = [xSelect[0],ySelect[0],
                                            xSelect[1],ySelect[1]]

                self.plotInfo.areaSelect = [xArea[0],yArea[0],xArea[1],yArea[1]]

        else:
            return

        for i in range(self.plotInfo.dataCount):

            dataSubset = (self.plotInfo.data[i][self.plotInfo.areaSelect[0]:
                                               self.plotInfo.areaSelect[2],
                                               self.plotInfo.areaSelect[1]:
                                               self.plotInfo.areaSelect[3]])

            self.plotInfo.sliceData[i] = (sum(dataSubset,
                                 axis = (self.plotInfo.horizontal)).flatten())

        pstart = self.plotInfo.axisSelect[(abs(self.plotInfo.horizontal-1))]
        pstop = (self.plotInfo.axisSelect[(abs(self.plotInfo.horizontal-1))+2])
        pstep = size(self.plotInfo.sliceData[0])
        self.plotInfo.sliceX = linspace(pstart,pstop,pstep)

        avgStart = self.plotInfo.axisSelect[self.plotInfo.horizontal+2]
        avgEnd = self.plotInfo.axisSelect[self.plotInfo.horizontal]

        title = ('Averaged from ' + str('%.2e '%avgEnd) + 'to ' +
             str('%.2e '%avgStart))
        if self.plotInfo.horizontal == True: title += ' Vertically'
        else: title += ' Horizontally'

        floorSliceData = asarray(copy(self.plotInfo.sliceData))
        floorSliceData[floorSliceData < self.plotInfo.vlimit[0]] = (
                                                    self.plotInfo.vlimit[0])
        
        self.slicePan.updatePlot(floorSliceData,
                         self.plotInfo.sliceX,
                         self.plotInfo.titles,
                         self.plotInfo.axisLabel[1-self.plotInfo.horizontal],
                         'Intensity',title,
                         self.plotInfo.scale,self.plotInfo.legTog)

        self.getCurrPage(None)

        return

    def boxUpdate(self,event):

        self.getPos(event)
        self.currLocalx = self.lastLocalx
        self.currValuex = self.lastValuex
        self.currLocaly = self.lastLocaly
        self.currValuey = self.lastValuey

        if (event.inaxes != self.probeAx):
            self.currLocalx = self.lastLocalx
            self.currValuex = self.lastValuex
            self.currLocaly = self.lastLocaly
            self.currValuey = self.lastValuey
        else:
            self.currLocalx = event.x
            self.lastLocalx = event.x

            self.currValuex = event.xdata
            self.lastValuex = event.xdata

            self.currLocaly = event.y
            self.lastLocaly = event.y

            self.currValuey = event.ydata
            self.lastValuey = event.ydata

        self.currPage.box.draw_rubberband(self,self.start[0],self.start[1],
                                          self.currLocalx,self.currLocaly)

        return



    def imPanBuild(self,selMask = None):
        if selMask == None:
            selMask = ones(self.plotInfo.dataCount)

        viewPan = MultiViewPanel(self.notePan,-1,self.plotInfo.data,
                                 self.plotInfo.extent,self.plotInfo.titles,
                                 self.plotInfo.axisLabel,self.plotInfo.vlimit,
                                 selMask,self.plotInfo.imScale)

        return viewPan

    def initNotebook(self):
        self.notePan.DeleteAllPages()
        if self.plotInfo.dataCount == 1.0:
            tabs = [None]
        else:
            tabs = [None]*(self.plotInfo.dataCount+1)

        allTab = self.imPanBuild()
        allTab.Fit()

        if self.plotInfo.dataCount > 1.0:
            tabs = [None]*(self.plotInfo.dataCount+1)
            tabs[0] = allTab
            self.notePan.AddPage(tabs[0],'Thumbnail')

            for i in range(self.plotInfo.dataCount):
                selMask = zeros(self.plotInfo.dataCount)
                selMask[i] = 1
                tabs[i+1] = self.imPanBuild(selMask)

                self.notePan.AddPage(tabs[i+1], str(self.plotInfo.titles[i]))

        else:
            self.notePan.AddPage(allTab,'Data')

        self.noteHorz = wx.BoxSizer(wx.HORIZONTAL)
        self.noteVert = wx.BoxSizer(wx.VERTICAL)

        self.noteHorz.Add(self.notePan,1,wx.EXPAND|wx.ALL,border = 1)
        self.noteVert.Add(self.noteHorz,1,wx.EXPAND|wx.LEFT|wx.RIGHT,border = 1)
        self.notePan.SetSizer(self.noteVert)

        return

    def toggleOrientation(self,event):

        self.plotInfo.horizontal = event.GetSelection()
        return

    def toggleScale(self,event):
        choice = ['linear','log']
        idx =  event.GetSelection()
        self.plotInfo.scale = choice[idx]
        self.slicePan.sliceAxis.set_yscale(str(choice[idx]))
        draw()

    def toggleLegend(self,event):

        idx =  event.GetSelection()
        self.plotInfo.legTog = idx
        if idx == 0:
            self.slicePan.sliceAxis.leg.set_visible(False)

        else:
            self.slicePan.sliceAxis.leg.set_visible(True)
        draw()

class PlotInfo(object):
    def __init__(self,rawData,vlimit,titles,extent,step,n,axisLabel):

        self.data=[]
        self.type=[]
        
        for i in range(len(rawData)):self.data.append(rawData[i][0])
        for i in range(len(rawData)):self.type.append(rawData[i][1])

        self.data = asarray(self.data)
        self.type = asarray(self.type,dtype = object)

        self.preservedData = array(copy(self.data))

        if vlimit == None:
            self.initVlim()
        else:
            self.vlimit = vlimit

        self.titles = titles
        self.extent = extent
        self.step = step
        self.n = n
        self.axisLabel = axisLabel
        self.selecting = False
        self.horizontal = False
        self.linear = True
        self.scale = 'linear'
        self.imScale = 'log'
        self.legTog = True
        self.dataCount = (shape(self.data))[0]
        self.areaSelect = None
        self.axisSelect = None

        self.sliceData = [zeros(30)]*self.dataCount
        self.sliceX = arange(30.0)
        return


    def initVlim(self):

        self.vlimit= [amax(self.data[0])*1.0e-20,amax(self.data)]
        return


def _test():

    from omfLoader import *

    import sample_prep,scatter
    from numpy import flipud,fliplr,rot90
    Au = (sample_prep.Parallelapiped(SLD = 4.506842e-6,
                                     Ms = 8.6e5,dim=[5.0e4,5.0e4,2.0e4]))

    Cr = (sample_prep.Layer(SLD = 3.01e-6,Ms = 8.6e5,
                            thickness_value = 1000.0))

    #Au.on_top_of(Cr)
    scene = sample_prep.Scene([Au])

    GeoUnit = (sample_prep.GeomUnit(Dxyz = [10.0e4,10.0e4,2.2e4],
                                    n = [20,25,20],scene = scene))
    cell = GeoUnit.buildUnit()
    space = sample_prep.Q_space([-.0001,-0.001,0.00002],
                                [.0001,0.001,0.01],[200,5,200])

    lattice = sample_prep.Rectilinear([20,20,1],cell)
    beam = sample_prep.Beam(5.0,None,None,0.05,None)
    test = scatter.Calculator(lattice, beam, space, cell)
    extent = space.getExtent()
    test.BA()
    test.resolution_correction()
    baData = test.results
    resbaData = test.corrected_results
    a = arange(0,50)
    b = arange(0,100)
    c = a.reshape(1,50)
    d = b.reshape(100,1)

    data = asarray((c*d),dtype = 'float')

    test = [None]*4
    test[0] = data
    test[1] = (fliplr(rot90(data))).reshape(100,50)
    test[2] = flipud(data)
    test[3] = fliplr(data)

    label = ['Test 1','Test 2','Test 3','Test 4']
    label = ['BA','resBA']
    axLab = ['xaxis','yaxis']
    MultiView([[baData,'Theory'],[resbaData,'Theory']],
              [space.q_step[0],space.q_step[2]],[space.points[0],
             space.points[2]],titles=label,extent= extent,axisLabel = axLab)

if __name__=="__main__":_test()