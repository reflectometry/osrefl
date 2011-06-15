"""
Select a region of the graph.

ginput_rect(ax) selects a region from the axes and continues with the plot.

s = BlockingRectangleSelector(ax) adds a rectangle selector to the axes and
lets the script call s.select() whenever a new region is needed.

demo() shows the selector in action.
"""
from matplotlib.pyplot import gca
from matplotlib.widgets import RectangleSelector
from matplotlib.blocking_input import BlockingInput
import pylab
class BlockingRectangleSelector:
    """
    Blocking rectangle selector selects once then continues with script.
    """
    def __init__(self, ax=None):
        if ax is None: ax=gca()
        self.ax = ax

        # drawtype is 'box' or 'line' or 'none'
        self.selector = RectangleSelector(self.ax, self._callback,
                               drawtype='box',useblit=True,
                               minspanx=5,minspany=5,spancoords='pixels')
        self.selector.set_active(False)
        self.block = BlockingInput(self.ax.figure)


    def _callback(self, event1, event2):
        """
        Selection callback.  event1 and event2 are the press and release events
        """
        x1, y1 = event1.xdata, event1.ydata
        x2, y2 = event2.xdata, event2.ydata
        if x1>x2: x1,x2 = x2,x1
        if y1>y2: y1,y2 = y2,y1
        self.x1,self.x2,self.y1,self.y2 = x1,x2,y1,y2
        print 'stopping event loop'
        self.ax.figure.canvas.stop_event_loop_default()


    def select(self):
        """
        Wait for box to be selected on the axes.
        """

        # Wait for selector to complete
        self.selector.set_active(True)
        self.ax.figure.canvas.draw_idle()
        self.block()
        self.selector.set_active(False)

        # Make sure the graph is redrawn next time the event loop is shown
        self.ax.figure.canvas.draw_idle()

    def remove(self):
        """
        Remove the selector from the axes.

        Note: this currently does nothing since matplotlib doesn't allow
        widgets to be removed from axes.
        """
        pylab.close('all')



def ginput_rect(ax=None,autoRemove = True):
    """
    Wait for user to select a region on the axes.

    Returns x1,x2,y1,y2
    """
    s = BlockingRectangleSelector(ax=ax)
    s.select()
    if autoRemove == True:
        s.remove()
    return s.x1,s.x2,s.y1,s.y2


def demo():
    from numpy import arange, sin, cos, pi
    from pylab import subplot, plot, axis, show
    from ginput_rect import ginput_rect

    subplot(111)
    # If N is large one can see improvement by blitting!
    N=100000
    x=10.0*arange(N)/(N-1)

    plot(x,sin(.2*pi*x),lw=3,c='b',alpha=.7)
    plot(x,cos(.2*pi*x),lw=3.5,c='r',alpha=.5)
    plot(x,-sin(.2*pi*x),lw=3.5,c='g',alpha=.3)
    print "\n      click  -->  release"

    x1,x2,y1,y2 = ginput_rect(autoRemove = False)
    print "(%3.2f, %3.2f) --> (%3.2f, %3.2f)"%(x1,y1,x2,y2)

    axis([x1, x2, y1, y2])
    show()


if __name__ == "__main__":
    demo()
