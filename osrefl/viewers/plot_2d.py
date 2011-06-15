# -*- coding: utf-8 -*-
from pylab import imshow,cm,colorbar,hot,show,xlabel,ylabel,connect, plot, figure, draw, axis, gcf,legend
from numpy import ones, sum, arange, transpose, log
import matplotlib.colors as colors
from matplotlib.widgets import  RectangleSelector
from colormap import change_colormap
from matplotlib.axis import XAxis, YAxis
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as Canvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg
from matplotlib.font_manager import fontManager, FontProperties
import wx
from matplotlib.image import FigureImage
from matplotlib.figure import Figure
from matplotlib.pyplot import figure, get_fignums
from zoom_colorbar import zoom_colorbar
from osrefl.loaders.reduction.cmapmenu import CMapMenu
from matplotlib.cm import get_cmap
import matplotlib.cbook as cbook
import matplotlib
#from binned_data_class3 import plottable_2d_data
#from wxrebin3 import rebinned_data
#import __main__

class MyNavigationToolbar(NavigationToolbar2WxAgg):
  """
  Extend the default wx toolbar with your own event handlers
  """

  class Cursors:
    HAND, POINTER, SELECT_REGION, MOVE = range(4)

  SET_SLICEMODE = wx.NewId()


  def __init__(self, canvas, cankill, data_frame_instance):
    NavigationToolbar2WxAgg.__init__(self, canvas)
    self.cursors = self.Cursors()
    self.data_frame_instance = data_frame_instance
    # for simplicity I'm going to reuse a bitmap from wx, you'll
    # probably want to add your own.
    slice_xpm_data = [
      "16 16 2 1",
      "       c None",
      ".      c #91384E",
      "                ",
      "                ",
      "    .           ",
      "    . .         ",
      "    .           ",
      "... . .  ..  .. ",
      "..  . . .   ... ",
      "  . . . .   .   ",
      "... . .  ..  .. ",
      "                ",
      "                ",
      "                ",
      "                ",
      "                ",
      "                ",
      "                "]
    slice2_xpm_data = [
      "32 16 72 1",
      "       c None",
      ".      c None",
      "+      c #91384E",
      "@      c #DFC6CC",
      "#      c #A3596B",
      "$      c #933B51",
      "%      c #9E5063",
      "&      c #D3B0B9",
      "*      c #FCFBFB",
      "=      c #C1909C",
      "-      c #99475B",
      ";      c #974459",
      ">      c #CDA5AF",
      ",      c #FDFCFC",
      "'      c #C495A1",
      ")      c #9A485C",
      "!      c #943F54",
      "~      c #B67B8A",
      "{      c #F9F4F5",
      "]      c #DBBEC5",
      "^      c #FBF9FA",
      "/      c #F2E7EA",
      "(      c #BB8592",
      "_      c #C597A2",
      ":      c #AF6E7E",
      "<      c #F2E9EB",
      "[      c #F8F2F3",
      "}      c #C18F9B",
      "|      c #C99DA8",
      "1      c #AE6D7D",
      "2      c #BA8290",
      "3      c #BE8996",
      "4      c #9D4E62",
      "5      c #D0AAB4",
      "6      c #9D4D61",
      "7      c #EEE1E4",
      "8      c #EEE0E3",
      "9      c #F5EDEF",
      "0      c #9A495D",
      "a      c #E8D5DA",
      "b      c #B27483",
      "c      c #99465B",
      "d      c #9F5265",
      "e      c #DDC2C9",
      "f      c #923A50",
      "g      c #FDFBFC",
      "h      c #923B50",
      "i      c #91394F",
      "j      c #FEFEFE",
      "k      c #DFC5CB",
      "l      c #9B4A5E",
      "m      c #F3E9EC",
      "n      c #B67C8B",
      "o      c #ECDEE1",
      "p      c #FCFAFA",
      "q      c #DEC4CA",
      "r      c #C697A3",
      "s      c #CA9FAA",
      "t      c #B07181",
      "u      c #EFE3E6",
      "v      c #EBDBDF",
      "w      c #D9BAC1",
      "x      c #A15669",
      "y      c #E0C7CD",
      "z      c #C08E9B",
      "A      c #98465A",
      "B      c #984559",
      "C      c #CEA6B0",
      "D      c #FEFDFD",
      "E      c #CAA0AB",
      "F      c #A35A6C",
      "G      c #D9BBC2",
      "................................",
      "................................",
      "................................",
      ".........+......................",
      ".........+..+...................",
      ".........+......................",
      "..@#$%&..+..+..*=-;>..,')!~{....",
      "..-]^/(..+..+.._:<[}..|1<{23....",
      "..45<....+..+..67.....%8..90....",
      "..abcde..+..+..fg.....hffi++....",
      "....jkl..+..+..67.....4m........",
      "..nopq-..+..+..r:<[}..stugv~....",
      "..wxh#y..+..+..*zABC..DE%$FG....",
      "................................",
      "................................",
      "................................"]

    slice_bitmap = wx.BitmapFromXPMData(slice2_xpm_data)
    self._slicemode = False
    self.AddCheckTool(self.SET_SLICEMODE, slice_bitmap,
                      shortHelp='Click me', longHelp='Activate slice mode')
    self.Bind(wx.EVT_TOOL, self.set_slicemode, id=self.SET_SLICEMODE)

  def set_slicemode(self, *args):
    'activate slice to rectangle mode'
    if self._slicemode:
        self._slicemode = False
        self.data_frame_instance.sliceplots_off()
    else:
        self._slicemode = True
        self.data_frame_instance.sliceplots_on()




class plot_2d_data(wx.Frame):
  """Generic 2d plotting routine - inputs are:
  - data (2d array of values),
  - x and y extent of the data,
  - title of graph, and
  - pixel mask to be used during summation  - must have same dimensions as data
  (only data entries corresponding to nonzero values in pixel_mask will be summed)
  - plot_title, x_label and y_label are added to the 2d-plot as you might expect"""

  def __init__(self, data, extent, caller = None, scale = 'log', window_title = 'log plot', pixel_mask = None, plot_title = "data plot", x_label = "x", y_label = "y", parent=None):
      wx.Frame.__init__(self, parent=None, title=window_title, pos = wx.DefaultPosition, size=wx.Size(800,600))
      print parent
      self.extent = extent
      self.data = data
      self.caller = caller
      self.window_title = window_title
      x_range = extent[0:2]
      #x_range.sort()
      self.x_min, self.x_max = x_range
      y_range = extent[2:4]
      #y_range.sort()
      self.y_min, self.y_max = y_range
      self.plot_title = plot_title
      self.x_label = x_label
      self.y_label = y_label
      self.slice_xy_range = (x_range, y_range)
      self.ID_QUIT = wx.NewId()
      self.ID_LOGLIN = wx.NewId()
      self.ID_UPCOLLIM = wx.NewId()
      self.ID_LOWCOLLIM = wx.NewId()

      menubar = wx.MenuBar()
      filemenu = wx.Menu()
      quit = wx.MenuItem(filemenu, 1, '&Quit\tCtrl+Q')
      #quit.SetBitmap(wx.Bitmap('icons/exit.png'))
      filemenu.AppendItem(quit)

      plotmenu = wx.Menu()
      self.menu_log_lin_toggle = plotmenu.Append(self.ID_LOGLIN, 'Plot 2d data with log color scale', 'plot 2d on log scale', kind=wx.ITEM_CHECK)
      self.Bind(wx.EVT_MENU, self.toggle_2d_plot_scale, id=self.ID_LOGLIN)
      menu_upper_colormap_limit = plotmenu.Append(self.ID_UPCOLLIM, 'Set upper limit of color map', 'Set upper limit of color map')
      self.Bind(wx.EVT_MENU, self.set_new_upper_color_limit, id=self.ID_UPCOLLIM)
      menu_lower_colormap_limit = plotmenu.Append(self.ID_LOWCOLLIM, 'Set lower limit of color map', 'Set lower limit of color map')
      self.Bind(wx.EVT_MENU, self.set_new_lower_color_limit, id=self.ID_LOWCOLLIM)
      #live_on_off = wx.MenuItem(live_update, 1, '&Live Update\tCtrl+L')
      #quit.SetBitmap(wx.Bitmap('icons/exit.png'))
      #live_update.AppendItem(self.live_toggle)
      #self.menu_log_lin_toggle.Check(True)

      menubar.Append(filemenu, '&File')
      menubar.Append(plotmenu, '&Plot')
      self.SetMenuBar(menubar)
      self.Centre()

      if pixel_mask == None:
        pixel_mask = ones(data.shape)

      if pixel_mask.shape != data.shape:
        print "Warning: pixel mask shape incompatible with data"
        pixel_mask = ones(data.shape)

      self.pixel_mask = pixel_mask

      self.show_data = transpose(data.copy())
      #self.minimum_intensity = self.data[pixel_mask.nonzero()].min()
      # correct for floating-point weirdness:
      self.minimum_intensity = self.data[self.data > 1e-17].min()

      #if scale == 'log':
        #self.show_data = log ( self.data.copy().T + self.minimum_intensity/2.0 )
        #self._scale = 'log'
        #self.menu_log_lin_toggle.Check(True)

      #elif (scale =='lin' or scale == 'linear'):
        #self._scale = 'lin'
        #self.menu_log_lin_toggle.Check(True)


      #self.bin_data = caller.bin_data
      #self.params = caller.params
      #fig = figure()
      self.fig = Figure(dpi=80, figsize=(5,5))
      #self.fig = figure()
      fig = self.fig
      self.canvas = Canvas(self, -1, self.fig)
      self.show_sliceplots = False # by default, sliceplots on
      self.sizer = wx.BoxSizer(wx.VERTICAL)
      self.sizer.Add(self.canvas, 1, wx.TOP | wx.LEFT | wx.EXPAND)

      #self.toolbar = Toolbar(self.canvas)
      self.toolbar = MyNavigationToolbar(self.canvas, True, self)
      self.toolbar.Realize()
      if wx.Platform == '__WXMAC__':
        # Mac platform (OSX 10.3, MacPython) does not seem to cope with
        # having a toolbar in a sizer. This work-around gets the buttons
        # back, but at the expense of having the toolbar at the top
        self.SetToolBar(self.toolbar)
      else:
        # On Windows platform, default window size is incorrect, so set
        # toolbar width to figure width.
        tw, th = self.toolbar.GetSizeTuple()
        fw, fh = self.canvas.GetSizeTuple()
        # By adding toolbar in sizer, we are able to put it at the bottom
        # of the frame - so appearance is closer to GTK version.
        # As noted above, doesn't work for Mac.
        self.toolbar.SetSize(wx.Size(fw, th))
        self.sizer.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)

      self.statusbar = self.CreateStatusBar()
      self.statusbar.SetFieldsCount(2)
      self.statusbar.SetStatusWidths([-1, -2])
      self.statusbar.SetStatusText("Current Position:", 0)

      self.canvas.mpl_connect('motion_notify_event', self.onmousemove)
      #self.canvas.mpl_connect('button_press_event', self.right_click_handler)
      #self.axes = fig.add_subplot(111)
      #self.axes = self.fig.gca()
      #ax = self.axes
      self.mapper = FigureImage(self.fig)
      #im = self.axes.pcolor(x,y,V,shading='flat')
      #self.mapper.add_observer(im)



      #self.show_data = transpose(log(self.show_data + self.minimum_intensity / 2.0))

      #self.canvas.mpl_connect('pick_event', self.log_lin_select)

      ax = fig.add_subplot(221, label='2d_plot')
      fig.sx = fig.add_subplot(222, label='sx', picker=True)
      fig.sx.xaxis.set_picker(True)
      fig.sx.yaxis.set_picker(True)
      fig.sx.yaxis.set_ticks_position('right')
      fig.sx.set_zorder(1)
      fig.sz = fig.add_subplot(223, label='sz', picker=True)
      fig.sz.xaxis.set_picker(True)
      fig.sz.yaxis.set_picker(True)
      fig.sz.set_zorder(1)
      self.RS = RectangleSelector(ax, self.onselect, drawtype='box', useblit=True)
      fig.slice_overlay = None

      ax.set_position([0.125,0.1,0.7,0.8])
      fig.cb = fig.add_axes([0.85,0.1,0.05,0.8])
      fig.cb.set_zorder(2)

      fig.ax = ax
      fig.ax.set_zorder(2)
      self.axes = ax
      ax.set_title(plot_title)
      #connect('key_press_event', self.toggle_selector)
      if scale == 'log':
        self.show_data = log ( self.data.copy().T + self.minimum_intensity/2.0 )
        self.__scale = 'log'
        self.fig.cb.set_xlabel('$\log_{10}I$')
        self.menu_log_lin_toggle.Check(True)

      elif (scale =='lin' or scale == 'linear'):
        self.__scale = 'lin'
        self.fig.cb.set_xlabel('$I$')
        self.menu_log_lin_toggle.Check(False)

      im = self.axes.imshow(self.show_data, interpolation='nearest', aspect='auto', origin='lower',cmap=cm.jet, extent=extent)
      #im = ax.imshow(data, interpolation='nearest', aspect='auto', origin='lower',cmap=cm.jet, extent=extent)
      fig.im = im
      ax.set_xlabel(x_label, size='large')
      ax.set_ylabel(y_label, size='large')
      self.toolbar.update()
      #zoom_colorbar(im)

      #fig.colorbar(im, cax=fig.cb)
      zoom_colorbar(im=im, cax=fig.cb)
      #figure(fig.number)
      #fig.canvas.draw()
      #return


      self.SetSizer(self.sizer)
      self.Fit()

      self.canvas.Bind(wx.EVT_RIGHT_DOWN, self.OnContext)
      self.Bind(wx.EVT_CLOSE, self.onExit)
      self.sliceplots_off()
      self.SetSize(wx.Size(800,600))
      self.canvas.draw()
      return

  def onExit(self, event):
    self.Destroy()

  def exit(self, event):
    wx.GetApp().Exit()


  def set_new_upper_color_limit(self, evt = None):
    current_uplim = self.fig.im.get_clim()[1]
    current_lowlim = self.fig.im.get_clim()[0]
    dlg = wx.TextEntryDialog(None, "Change upper limit of color map (currently %f)" % current_uplim, defaultValue = "%f" % current_uplim)
    if dlg.ShowModal() == wx.ID_OK:
      new_val = dlg.GetValue()
      xlab = self.fig.cb.get_xlabel()
      ylab = self.fig.cb.get_ylabel()
      self.fig.im.set_clim((current_lowlim, float(new_val)))
      self.fig.cb.set_xlabel(xlab)
      self.fig.cb.set_ylabel(ylab)
      self.fig.canvas.draw()
    dlg.Destroy()

  def set_new_lower_color_limit(self, evt = None):
    current_uplim = self.fig.im.get_clim()[1]
    current_lowlim = self.fig.im.get_clim()[0]
    dlg = wx.TextEntryDialog(None, "Change lower limit of color map (currently %f)" % current_lowlim, defaultValue = "%f" % current_lowlim)
    if dlg.ShowModal() == wx.ID_OK:
      new_val = dlg.GetValue()
      xlab = self.fig.cb.get_xlabel()
      ylab = self.fig.cb.get_ylabel()
      self.fig.im.set_clim((float(new_val), current_uplim))
      self.fig.cb.set_xlabel(xlab)
      self.fig.cb.set_ylabel(ylab)
      self.fig.canvas.draw()
    dlg.Destroy()

  def OnContext(self, evt):
    print self.show_sliceplots
    mpl_x = evt.X
    mpl_y = self.fig.canvas.GetSize()[1] - evt.Y
    mpl_mouseevent = matplotlib.backend_bases.MouseEvent('button_press_event', self.canvas, mpl_x, mpl_y, button = 3)

    if (mpl_mouseevent.inaxes == self.fig.ax):
      self.area_context(mpl_mouseevent, evt)
    elif ((mpl_mouseevent.inaxes == self.fig.sx or mpl_mouseevent.inaxes == self.fig.sz) and (self.show_sliceplots == True)):
      self.lineplot_context(mpl_mouseevent, evt)

  def area_context(self, mpl_mouseevent, evt):
    area_popup = wx.Menu()
    item1 = area_popup.Append(wx.ID_ANY,'&Grid on/off', 'Toggle grid lines')
    wx.EVT_MENU(self, item1.GetId(), self.OnGridToggle)
    cmapmenu = CMapMenu(self, callback = self.OnColormap, mapper=self.mapper, canvas=self.canvas)
    item2 = area_popup.Append(wx.ID_ANY,'&Toggle log/lin', 'Toggle log/linear scale')
    wx.EVT_MENU(self, item2.GetId(), lambda evt: self.toggle_log_lin(mpl_mouseevent))
    item3 = area_popup.AppendMenu(wx.ID_ANY, "Colourmaps", cmapmenu)
    self.PopupMenu(area_popup, evt.GetPositionTuple())

  def figure_list_dialog(self):
    figure_list = get_fignums()
    figure_list_names = []
    for fig in figure_list:
      figure_list_names.append('Figure ' + str(fig))
    figure_list_names.insert(0, 'New Figure')
    figure_list.insert(0, None)
    #selection_num = wx.GetSingleChoiceIndex('Choose other plot', '', other_plot_names)
    dlg = wx.SingleChoiceDialog(None, 'Choose figure number', '', figure_list_names)
    dlg.SetSize(wx.Size(640,480))
    if dlg.ShowModal() == wx.ID_OK:
      selection_num=dlg.GetSelection()
    dlg.Destroy()
    print selection_num
    return figure_list[selection_num]

  def lineplot_context(self, mpl_mouseevent, evt):
    popup = wx.Menu()
    item1 = popup.Append(wx.ID_ANY,'&Toggle log/lin', 'Toggle log/linear scale of slices')
    wx.EVT_MENU(self, item1.GetId(), lambda evt: self.toggle_log_lin(mpl_mouseevent))
    if mpl_mouseevent.inaxes == self.fig.sx:
      item2 = popup.Append(wx.ID_ANY, "Save x slice", "save this slice")
      wx.EVT_MENU(self, item2.GetId(), self.save_x_slice)
      item3 = popup.Append(wx.ID_ANY, '&Popout plot', 'Open this data in a figure window')
      wx.EVT_MENU(self, item3.GetId(), lambda evt: self.popout_x_slice())
    elif mpl_mouseevent.inaxes == self.fig.sz:
      item2 = popup.Append(wx.ID_ANY, "Save y slice", "save this slice")
      wx.EVT_MENU(self, item2.GetId(), self.save_y_slice)
      item3 = popup.Append(wx.ID_ANY, '&Popout plot', 'Open this data in a new plot window')
      wx.EVT_MENU(self, item3.GetId(), lambda evt: self.popout_y_slice())
    self.PopupMenu(popup, evt.GetPositionTuple())


  def popout_y_slice(self, event=None, figure_num = None, label = None):
    if figure_num == None:
        figure_num = self.figure_list_dialog()
    fig = figure(figure_num) # if this is None, matplotlib automatically increments figure number to highest + 1
    ax = self.fig.sz
    slice_desc = '\nsliceplot([%f,%f],[%f,%f])' % (self.slice_xy_range[0][0],self.slice_xy_range[0][1],self.slice_xy_range[1][0],self.slice_xy_range[1][1])
    if figure_num == None:
        default_title = self.plot_title + slice_desc
        dlg = wx.TextEntryDialog(None, 'Enter title for plot', defaultValue = default_title)
        if dlg.ShowModal() == wx.ID_OK:
            title = dlg.GetValue()
        else:
            title = default_title
        dlg.Destroy()
        new_ax = fig.add_subplot(111)
        new_ax.set_title(title, size='large')
        new_ax.set_xlabel(self.x_label, size='x-large')
        new_ax.set_ylabel('$I_{summed}$', size='x-large')
    else:
        new_ax = fig.axes[0]
    if label == None:
        default_label = self.window_title + ': ' + self.plot_title + slice_desc
        dlg = wx.TextEntryDialog(None, 'Enter data label (for plot legend)', defaultValue = default_label)
        if dlg.ShowModal() == wx.ID_OK:
            label = dlg.GetValue()
        else:
            label = default_label
        dlg.Destroy()
    xy = ax.lines[0].get_data()
    x = xy[0]
    y = xy[1]
    new_ax.plot(x,y, label = label)
    font = FontProperties(size='small')
    lg = legend(prop=font)
    drag_lg = DraggableLegend(lg)
    drag_lg.connect()
    fig.canvas.draw()
    fig.show()

  def popout_x_slice(self, event=None, figure_num = None, label = None):
    if figure_num == None:
        figure_num = self.figure_list_dialog()
    fig = figure(figure_num)
    ax = self.fig.sx
    slice_desc = '\nsliceplot([%f,%f],[%f,%f])' % (self.slice_xy_range[0][0],self.slice_xy_range[0][1],self.slice_xy_range[1][0],self.slice_xy_range[1][1])
    if figure_num == None:
        default_title = self.plot_title + slice_desc
        dlg = wx.TextEntryDialog(None, 'Enter title for plot', defaultValue = default_title)
        if dlg.ShowModal() == wx.ID_OK:
            title = dlg.GetValue()
        else:
            title = default_title
        dlg.Destroy()
        new_ax = fig.add_subplot(111)
        new_ax.set_title(title, size='large')
        new_ax.set_xlabel(self.y_label, size='x-large')
        new_ax.set_ylabel('$I_{summed}$', size='x-large')
    else:
        new_ax = fig.axes[0]
    if label == None:
        default_label = self.window_title + ': ' + self.plot_title + slice_desc
        dlg = wx.TextEntryDialog(None, 'Enter data label (for plot legend)', defaultValue = default_label)
        if dlg.ShowModal() == wx.ID_OK:
            label = dlg.GetValue()
        else:
            label = default_label
        dlg.Destroy()
    xy = ax.lines[0].get_data()
    x = xy[1]
    y = xy[0]
    new_ax.plot(x,y, label = label)
    font = FontProperties(size='small')
    lg = legend(prop=font)
    drag_lg = DraggableLegend(lg)
    drag_lg.connect()
    fig.canvas.draw()
    fig.show()

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
    outFile.write('#xmin: ' + str(self.slice_xy_range[0][0]) + '\n')
    outFile.write('#xmax: ' + str(self.slice_xy_range[0][1]) + '\n')
    outFile.write('#ymin: ' + str(self.slice_xy_range[1][0]) + '\n')
    outFile.write('#ymax: ' + str(self.slice_xy_range[1][1]) + '\n')
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


  def OnGridToggle(self, event):
      self.fig.ax.grid()
      self.fig.canvas.draw_idle()

  def OnColormap(self, name):
      print "Selected colormap",name
      self.fig.im.set_cmap(get_cmap(name))
      self.fig.canvas.draw()

  def toggle_2d_plot_scale(self, event=None):
    if self.__scale == 'log':
      self.show_data = self.data.T
      self.fig.im.set_array(self.show_data)
      self.fig.im.autoscale()
      self.fig.cb.set_xlabel('$I$')
      self.__scale = 'lin'
      self.menu_log_lin_toggle.Check(False)
      self.statusbar.SetStatusText("%s scale" % self.__scale, 0)
      self.fig.canvas.draw_idle()
    elif self.__scale == 'lin':
      self.show_data = log ( self.data.copy().T + self.minimum_intensity/2.0 )
      self.fig.im.set_array(self.show_data)
      self.fig.im.autoscale()
      self.fig.cb.set_xlabel('$\log_{10}I$')
      self.__scale = 'log'
      self.menu_log_lin_toggle.Check(True)
      self.statusbar.SetStatusText("%s scale" % self.__scale, 0)
      self.fig.canvas.draw_idle()


  def toggle_log_lin(self,event):

    ax = event.inaxes
    label = ax.get_label()

    if label == '2d_plot':
      self.toggle_2d_plot_scale()

    if label == 'sz':
      scale = ax.get_yscale()
      if scale == 'log':
        ax.set_yscale('linear')
        ax.figure.canvas.draw_idle()
      elif scale == 'linear':
        ax.set_yscale('log')
        ax.figure.canvas.draw_idle()

    elif label == 'sx':
      scale = ax.get_xscale()
      if scale == 'log':
        ax.set_xscale('linear')
        ax.figure.canvas.draw_idle()
      elif scale == 'linear':
        ax.set_xscale('log')
        ax.figure.canvas.draw_idle()


  def onmousemove(self,event):
    # the cursor position is given in the wx status bar
    #self.fig.gca()
    if event.inaxes:
        x, y = event.xdata, event.ydata
        self.statusbar.SetStatusText("%s scale x = %.3g, y = %.3g" % (self.__scale,x,y), 1)
        #self.statusbar.SetStatusText("y = %.3g" %y, 2)


  def onselect(self, eclick, erelease):
    x_range = [eclick.xdata, erelease.xdata]
    y_range = [eclick.ydata, erelease.ydata]
    ax = eclick.inaxes
    self.sliceplot((x_range, y_range), ax)
    print 'sliceplot(([%f,%f],[%f,%f]))' % (x_range[0],x_range[1],y_range[0],y_range[1])

  def sliceplots_off(self):
    self.fig.ax.set_position([0.125,0.1,0.7,0.8])
    self.fig.cb.set_position([0.85,0.1,0.05,0.8])
    #self.fig.cb.set_visible(True)
    self.fig.sx.set_visible(False)
    self.fig.sz.set_visible(False)
    if self.fig.slice_overlay:
      self.fig.slice_overlay[0].set_visible(False)
    self.RS.set_active(False)
    self.show_sliceplots = False
    self.fig.canvas.draw()

  def sliceplots_on(self):
    self.fig.ax.set_position([0.125,0.53636364, 0.35227273,0.36363636])
    self.fig.cb.set_position([0.49,0.53636364, 0.02, 0.36363636])
    self.fig.sx.set_position([0.58,0.53636364, 0.35227273,0.36363636])
    self.fig.sx.set_visible(True)
    self.fig.sz.set_visible(True)
    #self.fig.cb.set_visible(False)
    if self.fig.slice_overlay:
      self.fig.slice_overlay[0].set_visible(True)
    self.RS.set_active(True)
    self.show_sliceplots = True
    self.fig.canvas.draw()

  def toggle_sliceplots(self):
    """switch between views with and without slice plots"""
    if self.show_sliceplots == True:
      self.sliceplots_off()
    else: # self.show_sliceplots == False
      self.sliceplots_on()

  def show_slice_overlay(self, x_range, y_range, x, slice_y_data, y, slice_x_data):
    """sum along x and z within the box defined by qX- and qZrange.
    sum along qx is plotted to the right of the data,
    sum along qz is plotted below the data.
    Transparent white rectangle is overlaid on data to show summing region"""
    from matplotlib.ticker import FormatStrFormatter, ScalarFormatter

    if self.fig == None:
      print('No figure for this dataset is available')
      return

    fig = self.fig
    ax = fig.ax
    extent = fig.im.get_extent()

    if fig.slice_overlay == None:
      fig.slice_overlay = ax.fill([x_range[0],x_range[1],x_range[1],x_range[0]],[y_range[0],y_range[0],y_range[1],y_range[1]],fc='white', alpha=0.3)
      fig.ax.set_ylim(extent[2],extent[3])
    else:
      fig.slice_overlay[0].xy = [(x_range[0],y_range[0]), (x_range[1],y_range[0]), (x_range[1],y_range[1]), (x_range[0],y_range[1])]
    fig.sz.clear()
    default_fmt = ScalarFormatter(useMathText=True)
    default_fmt.set_powerlimits((-2,4))
    fig.sz.xaxis.set_major_formatter(default_fmt)
    fig.sz.yaxis.set_major_formatter(default_fmt)
    fig.sz.xaxis.set_major_formatter(FormatStrFormatter('%.2g'))
    fig.sz.set_xlim(x[0], x[-1])
    fig.sz.plot(x, slice_y_data)
    fig.sx.clear()
    fig.sx.yaxis.set_major_formatter(default_fmt)
    fig.sx.xaxis.set_major_formatter(default_fmt)
    fig.sx.yaxis.set_ticks_position('right')
    fig.sx.yaxis.set_major_formatter(FormatStrFormatter('%.2g'))
    fig.sx.set_ylim(y[0], y[-1])
    fig.sx.plot(slice_x_data, y)

    fig.im.set_extent(extent)
    fig.canvas.draw()

  def copy_intensity_range_from(self, other_plot):
    if isinstance(other_plot, type(self)):
      xlab = self.fig.cb.get_xlabel()
      ylab = self.fig.cb.get_ylabel()

      self.fig.im.set_clim(other_plot.fig.im.get_clim())
      self.fig.cb.set_xlabel(xlab)
      self.fig.cb.set_ylabel(ylab)
      self.fig.canvas.draw()

  def sliceplot(self, xy_range, ax = None):
    """sum along x and z within the box defined by qX- and qZrange.
    sum along qx is plotted to the right of the data,
    sum along qz is plotted below the data.
    Transparent white rectangle is overlaid on data to show summing region"""
    self.sliceplots_on()
    x_range, y_range = xy_range
    x, slice_y_data, y, slice_x_data = self.do_xy_slice(x_range, y_range)
    self.x = x
    self.slice_y_data = slice_y_data
    self.y = y
    self.slice_x_data = slice_x_data
    self.slice_xy_range = xy_range

    self.show_slice_overlay(x_range, y_range, x, slice_y_data, y, slice_x_data)

  def do_xy_slice(self, x_range, y_range):
    """ slice up the data, once along x and once along z.
    returns 4 arrays:  a y-axis for the x data,
    an x-axis for the y data."""
    #params = self.params
    print 'doing xy slice'
    data = self.data
    pixels = self.pixel_mask
    # zero out any pixels in the sum that have zero in the pixel count:
    data[pixels == 0] = 0

    normalization_matrix = ones(data.shape)
    normalization_matrix[pixels == 0] = 0
    x_min = min(x_range)
    x_max = max(x_range)
    y_min = min(y_range)
    y_max = max(y_range)

    x_size,y_size = data.shape
    global_x_range = (self.x_max - self.x_min)
    global_y_range = (self.y_max - self.y_min)

    x_pixel_min = round( (x_min - self.x_min) / global_x_range * x_size )
    x_pixel_max = round( (x_max - self.x_min) / global_x_range * x_size )
    y_pixel_min = round( (y_min - self.y_min) / global_y_range * y_size )
    y_pixel_max = round( (y_max - self.y_min) / global_y_range * y_size )

    #correct any sign switches:
    if (x_pixel_min > x_pixel_max):
      new_min = x_pixel_max
      x_pixel_max = x_pixel_min
      x_pixel_min = new_min

    if (y_pixel_min > y_pixel_max):
      new_min = y_pixel_max
      y_pixel_max = y_pixel_min
      y_pixel_min = new_min

    new_x_min = x_pixel_min / x_size * global_x_range + self.x_min
    new_x_max = x_pixel_max / x_size * global_x_range + self.x_min
    new_y_min = y_pixel_min / y_size * global_y_range + self.y_min
    new_y_max = y_pixel_max / y_size * global_y_range + self.y_min

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

class DraggableLegend:
    lock = None  # only one can be animated at a time
    def __init__(self, leg):
        self.leg = leg
        self.frame = leg.get_frame()
        self.figure = self.leg.figure
        self.axes = self.leg.parent
        self.press = None
        self.background = None

    def connect(self):
        'connect to all the events we need'
        self.cidpress = self.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)

    def on_press(self, event):
        'on button press we will see if the mouse is over us and store some data'
        if event.inaxes != self.axes: return
        if DraggableLegend.lock is not None: return
        contains, attrd = self.frame.contains(event)
        if not contains: return

        #bbox = self.leg.get_bbox_to_anchor().transformed(self.axes.transData.inverted())
        bbox = self.leg.get_bbox_to_anchor()
        #print 'event contains', self.leg.get_window_extent()
        self.press = bbox, event.x, event.y
        DraggableLegend.lock = self

        self.event = event
        # draw everything but the selected legend and store the pixel buffer
        canvas = self.leg.figure.canvas
        axes = self.axes
        self.leg.set_animated(True)
        canvas.draw()
        self.background = canvas.copy_from_bbox(self.axes.bbox)

        # now redraw just the rectangle
        axes.draw_artist(self.leg)

        # and blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_motion(self, event):
        'on motion we will move the legend if the mouse is over us'
        if DraggableLegend.lock is not self:
            return
        if event.inaxes != self.axes: return
        bbox, xpress, ypress = self.press
        dx = event.x - xpress
        dy = event.y - ypress

        new_bbox = (bbox.translated(dx, dy)).transformed(self.axes.transAxes.inverted())

        self.leg.set_bbox_to_anchor(new_bbox, transform = self.axes.transAxes)
        #new_bbox = bbox.translated(dx, dy)
        #self.leg.set_bbox_to_anchor(new_bbox, transform = self.transform)

        canvas = self.figure.canvas
        axes = self.axes
        # restore the background region
        canvas.restore_region(self.background)

        # redraw just the current rectangle
        axes.draw_artist(self.leg)

        # blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_release(self, event):
        'on release we reset the press data'
        if DraggableLegend.lock is not self:
            return

        self.press = None
        DraggableLegend.lock = None

        # turn off the rect animation property and reset the background
        self.leg.set_animated(False)
        self.background = None

        # redraw the full figure
        self.figure.canvas.draw()

    def disconnect(self):
        'disconnect all the stored connection ids'
        self.figure.canvas.mpl_disconnect(self.cidpress)
        self.figure.canvas.mpl_disconnect(self.cidrelease)
        self.figure.canvas.mpl_disconnect(self.cidmotion)


