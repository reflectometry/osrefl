class zoom_colorbar:
    """ a zoomable colorbar based on matplotlib colorbar
    (changes the color mapping of the image it is attached to)
    NOTE: also enables mouse zoom on image it is attached to.
    """

    def __init__(self, im=None):
        """if figure.image is not passed as parameter,
        use current image"""
        from pylab import colorbar, gci
        if im == None:
            im = gci()

        ax = im.axes
        cb = colorbar(im, ax=ax)
        self.colorbar_instance = cb
        canvas = ax.figure.canvas
        if not hasattr(canvas,'_colorbars'):
            canvas._colorbars = {}
        canvas._colorbars[cb.ax] = cb
        canvas.mpl_connect('scroll_event', self.onWheel)
        return

    def onWheel(self,event):
        """
        Process mouse wheel as zoom events
        """
        from matplotlib import colors

        ax = event.inaxes
        try:
            step = event.step
        except:
            if event.button == 'up':
                step = 1
            else:
                step = -1
        #print "zoom",step

        # Ick! Can't tell if the axes contains a colorbar.
        if hasattr(event.canvas,'_colorbars') and ax in event.canvas._colorbars:
            mappable = event.canvas._colorbars[ax].mappable
            # rescale colormap: the axes are already scaled to 0..1,
            # so use bal instead of pt for centering
            lo,hi = mappable.get_clim()
            if isinstance(mappable.norm, colors.LogNorm):
                vscale = 'log'
            else:
                vscale = 'linear'
            lo,hi = self._rescale(lo,hi,step,bal=event.ydata,scale=vscale)
            mappable.set_clim(lo,hi)

        elif ax != None:
            # Event occurred inside a plotting area
            lo,hi = ax.get_xlim()
            lo,hi = self._rescale(lo,hi,step,pt=event.xdata,scale=ax.get_xscale())
            ax.set_xlim((lo,hi))

            lo,hi = ax.get_ylim()
            lo,hi = self._rescale(lo,hi,step,pt=event.ydata,scale=ax.get_yscale())
            ax.set_ylim((lo,hi))
        else:
            # Check if zoom happens in the axes
            xdata,ydata = None,None
            x,y = event.x,event.y
            for ax in event.canvas.figure.get_axes():
                if ax.xaxis.contains(event):
                    xdata,_ = ax.transAxes.inverted().transform_point((x,y))
                    #print "xaxis",x,"->",xdata
                if ax.yaxis.contains(event):
                    _,ydata = ax.transAxes.inverted().transform_point((x,y))
                    #print "yaxis",y,"->",ydata
            if xdata is not None:
                lo,hi = ax.get_xlim()
                lo,hi = self._rescale(lo,hi,step,bal=xdata,scale=ax.get_xscale())
                ax.set_xlim((lo,hi))
            if ydata is not None:
                lo,hi = ax.get_ylim()
                lo,hi = self._rescale(lo,hi,step,bal=ydata,scale=ax.get_yscale())
                ax.set_ylim((lo,hi))

        event.canvas.draw_idle()

    def _rescale(self,lo,hi,step,pt=None,bal=None,scale='linear'):
        """
        Rescale (lo,hi) by step, returning the new (lo,hi)
        The scaling is centered on pt, with positive values of step
        driving lo/hi away from pt and negative values pulling them in.
        If bal is given instead of point, it is already in [0,1] coordinates.

        This is a helper function for step-based zooming.
        """
        # Convert values into the correct scale for a linear transformation
        # TODO: use proper scale transformers
        if scale=='log':
            lo,hi = math.log10(lo),math.log10(hi)
            if pt is not None: pt = math.log10(pt)

        # Compute delta from axis range * %, or 1-% if percent is negative
        if step > 0:
            delta = float(hi-lo)*step/100
        else:
            delta = float(hi-lo)*step/(100-step)

        # Add scale factor proportionally to the lo and hi values, preserving the
        # point under the mouse
        if bal is None:
            bal = float(pt-lo)/(hi-lo)
        lo = lo - bal*delta
        hi = hi + (1-bal)*delta

        # Convert transformed values back to the original scale
        if scale=='log':
            lo,hi = math.pow(10.,lo),math.pow(10.,hi)

        return (lo,hi)
