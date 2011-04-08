# Copyright (C) 2008 University of Maryland
# All rights reserved.
# See LICENSE.txt for details.
# Author: Christopher Metting

#Starting Date:8/23/2010


from numpy import size, shape, array, log, zeros, ones, argwhere, asarray, ma, isneginf,log,inf,nonzero
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
from pylab import figure,show,imshow, draw, delaxes,cla,subplot,setp,plot,figlegend,legend,xlim


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
    
    '''
    def __init__(self,parent):
        wx.Frame.__init__(self,parent,title = "Plot Comparison Application")

        return

class ParentVisPanel(wx.Panel):
    def __init__(self,parent):
        wx.Panel.__init__(self,parent,style = wx.BORDER_RAISED)
        #self.SetBackgroundColour('purple')
        
class plotNotebook(wx.Notebook):
    def __init__(self,parent):
        wx.Notebook.__init__(self, parent,style = wx.BORDER_RAISED)
        #self.SetBackgroundColour('blue')

class MultiViewPanel(wx.Panel):
    def __init__(self,parent,id, data, extent, titles, axLabel, vlimit,
                 selMask,scale = 'log'):
        
        wx.Panel.__init__(self,parent,id,style = wx.BORDER_SUNKEN)
        
        FIGPROPS = dict(figsize = (7.0,(7.0/1.618)),dpi = 128)
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
            if vlimit[0] == None:
    
                vmin = amin(data[0])
            else:
                vmin = vlimit[0]
                
            if vlimit[1] == None:
                vmax = amax(data[0])
            else:
                vmax = vlimit[1]
                
        elif (scale == 'log'):
            plotData = log(data)
                
            if (vlimit[1] == None):
                vmax = log(amax(data[0][nonzero(data[0])]))
            else:
                vmax = log(vlimit[1])
            if (vlimit[0] == None):

                vmin = log(amax(data[0][nonzero(data[0])])) - 30.0

            else:
                vmin = log(vlimit[0])
                
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
        
        multiViewSizer.Add(self.plotCanvas,3,wx.EXPAND|wx.RIGHT|wx.LEFT,
                                                            border = 1)
        multiViewVert.Add(multiViewSizer,1,wx.EXPAND|wx.TOP|wx.BOTTOM,
                                                            border = 1) 
        

        self.SetSizer(multiViewVert)
        
        
        return

class MouseReadPan(wx.Panel):
    def __init__(self,parent):
        wx.Panel.__init__(self,parent,style = wx.BORDER_RAISED)
        self.readoutLabel = wx.StaticText(self,-1)
        self.readoutLabel.SetLabel('    X              Y             Z ')
        
        self.readoutData = wx.StaticText(self,-1)
        newFont = self.readoutData.GetFont()
        newFont.SetPointSize(8)
        self.readoutData.SetFont(newFont)
        self.readoutData.SetLabel('  ')
        
        self.readSizer = wx.BoxSizer(wx.VERTICAL)
        self.readSizer.Add(self.readoutLabel)
        self.readSizer.Add(self.readoutData)
        self.SetSizer(self.readSizer)
        self.Fit()
          
class ScalePan(wx.Panel):
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
    def __init__(self,parent):
        wx.Panel.__init__(self,parent,style = wx.BORDER_RAISED)
        wx.StaticText
        self.vminLab = wx.StaticText(self,-1,'   v min')
        self.vmaxLab = wx.StaticText(self,-1,'   v max')

        
        newFont = self.vminLab.GetFont()
        newFont.SetPointSize(12)
        
        self.vminLab.SetFont(newFont)
        self.vmaxLab.SetFont(newFont)
        
        self.vminBox = wx.TextCtrl(self, 1,size = (10, 20))
        self.vmaxBox = wx.TextCtrl(self, 2,size = (10, 20))
        self.submit = wx.Button(self,3, 'Submit',(10, 10))
        self.reset = wx.Button(self,4, 'Reset',(10, 10))
        
        self.butSize = wx.BoxSizer(wx.HORIZONTAL) 
        self.panSize = wx.BoxSizer(wx.VERTICAL)
        self.labVal = wx.BoxSizer(wx.HORIZONTAL)
        self.enterVal = wx.BoxSizer(wx.HORIZONTAL)
        
        self.butSize.Add(self.submit)
        self.butSize.Add(self.reset)
        
        self.labVal.Add(self.vminLab,1,wx.EXPAND|wx.RIGHT|wx.LEFT|wx.BOTTOM,border = 1)
        self.labVal.Add(self.vmaxLab,1,wx.EXPAND|wx.RIGHT|wx.LEFT|wx.BOTTOM,border = 1)
        
        self.enterVal.Add(self.vminBox,1,wx.RIGHT|wx.LEFT,border = 1)
        self.enterVal.Add(self.vmaxBox,1,wx.RIGHT|wx.LEFT,border = 1)
        
        self.panSize.Add(self.labVal,1,wx.EXPAND|wx.RIGHT|wx.LEFT,border = 1)
        self.panSize.Add(self.enterVal,1,wx.EXPAND|wx.RIGHT|wx.LEFT,border = 1)
        self.panSize.Add(self.butSize,1,wx.EXPAND|wx.ALL,border = 3)
        self.SetSizer(self.panSize)
        
        self.Fit()
        return
    
class HozVertPan(wx.Panel):
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
    def __init__(self,parent):
        wx.Panel.__init__(self,parent,style = wx.BORDER_SUNKEN)
        
        self.sliceFigure = Figure(frameon = True,figsize =(12,3))
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
        sliceSizer.Add(self.sliceCanvas,-1,wx.EXPAND|wx.RIGHT|wx.LEFT,border = 1)
        self.SetSizer(sliceSizer)
        return
    
    def updatePlot(self,data,xaxis,legLab,xlabel,ylabel,title,scale,legTog):

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
    
class PlotCtrler(object):
    '''
    ***********CONTROLLER********************
    '''
    def __init__(self, app,plotInfo):
        self.plotInfo = plotInfo
        
        self.Frame = MultiViewFrame(None)
        
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
        
        
        self.plotContrSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.toolSizer = wx.BoxSizer(wx.VERTICAL)
        
        self.toolSizer.Add(self.readOut,0,wx.EXPAND|wx.LEFT|wx.RIGHT,
                           border = 1)
        self.toolSizer.Add(self.orientation,0,wx.EXPAND|wx.LEFT|wx.RIGHT,
                           border = 1)
        self.toolSizer.Add(self.scale,0,wx.EXPAND|wx.LEFT|wx.RIGHT,
                           border = 1)
        self.toolSizer.Add(self.zlim,0,wx.EXPAND|wx.LEFT|wx.RIGHT,
                           border = 1)

        
        self.plotContrSizer.Add(self.slicePan,1,wx.EXPAND|wx.LEFT|wx.RIGHT,
                                border = 1)
        self.plotContrSizer.Add(self.toolSizer,0,wx.EXPAND|wx.LEFT|wx.RIGHT,
                                border = 1)
        
        self.vertSizer = wx.BoxSizer(wx.VERTICAL)
        
        self.vertSizer.Add(self.notePan,0,wx.EXPAND|wx.LEFT|wx.RIGHT,border = 1)
        self.vertSizer.Add(self.plotContrSizer,1,wx.EXPAND|wx.LEFT|wx.RIGHT,
                           border = 1)
        
        self.parPanel.SetSizer(self.vertSizer)
        
        self.slicePan.sliceCanvas.Bind(wx.EVT_RIGHT_DOWN,self.sliceHome)
        self.parPanel.Fit()
        self.Frame.SetClientSize(self.parPanel.GetSize())
        self.Frame.SetMinSize(self.parPanel.GetSize())
        self.Frame.Show(True)
        
        self.notePan.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED,self.getCurrPage)
        
        self.orientation.Bind(wx.EVT_RADIOBOX,self.toggleOrientation)
        self.scale.Bind(wx.EVT_RADIOBOX,self.toggleScale)
        
        self.zlim.Bind(wx.EVT_BUTTON,self.setZ,self.zlim.submit)
        self.zlim.Bind(wx.EVT_BUTTON,self.resetZ,self.zlim.reset)
        
        return
    
    def resetZ(self,event):
        self.plotInfo.vlimit = [None,None]
        curPage = self.notePan.GetSelection()
        self.imScaleChange(event)
        self.notePan.SetSelection(curPage)
        return
    
    def setZ(self,event):
        self.plotInfo.vlimit[0] = float(self.zlim.vminBox.GetValue())
        self.plotInfo.vlimit[1] = float(self.zlim.vmaxBox.GetValue())
        
        if self.plotInfo.imScale == 'log':
            self.plotInfo.vlimit[0] = exp(self.plotInfo.vlimit[0])
            self.plotInfo.vlimit[1] = exp(self.plotInfo.vlimit[1])
        curPage = self.notePan.GetSelection()
        self.imScaleChange(event)
        self.notePan.SetSelection(curPage)
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

            if (xcorrd >= (self.plotInfo.n[0]) or ycorrd >= (self.plotInfo.n[1])):
                self.readOut.readoutData.SetLabel(' ')
                
            else:
                z = self.plotInfo.data[int(event.inaxes.get_gid())][xcorrd,ycorrd]
                self.readOut.readoutData.SetLabel('%.2e  '%x +'%.2e  '%y +'%.2e  '%z)
 
        return
    
    
    def getCurrPage(self,event):
        self.currPage = self.notePan.GetCurrentPage()

        #for coordinates
        self.currPage.plotCanvas.mpl_connect('motion_notify_event',self.getCoord)
        
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
            coord[i] = int((value[i]+abs(min))/step)
        return sorted(coord)

    def getSliceEnd(self,event):
        self.currPage.box.draw_rubberband(self,0.0,0.0,0.0,0.0)
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

        self.slicePan.updatePlot(self.plotInfo.sliceData,
                                 self.plotInfo.sliceX,
                                 self.plotInfo.titles,
                                 self.plotInfo.axisLabel[0],
                                 self.plotInfo.axisLabel[1],title,
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
    def __init__(self,data,vlimit,titles,extent,step,n,axisLabel):
        self.data = data
        
        if vlimit == None:
            self.vlimit = [None,None]
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
    

       
def _test():

    from omfLoader import *
    from numpy import asarray,arange,reshape,shape
    
    '''
    mag = Omf('/home/mettingc/Documents/BBM/c_500mTz_-Oxs_MinDriver-Magnetization-05-0005651.omf')
    #mag = Omf('/home/mettingc/Documents/test.omf')
    momentView(mag
    '''
    
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
    space = sample_prep.Q_space([-.0001,-0.001,0.00002],[.0001,0.001,0.01],[200,5,200])
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
    MultiView([baData,resbaData],[space.q_step[0],space.q_step[2]],[space.points[0],space.points[2]],titles=label,extent= extent,axisLabel = axLab)

if __name__=="__main__":_test()
