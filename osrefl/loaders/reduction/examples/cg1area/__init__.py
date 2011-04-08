"""
Sample data from NCNR-AND/R.

This is area detector data with all frames preserved.

TODO: Document contents of this data file.
"""

import numpy, os
import osrefl.loaders.reduction as reflred

PATH = os.path.dirname(os.path.realpath(__file__))

data = reflred.load(os.path.join(PATH,'psdca022.cg1.gz'))
#data = reflred.load(os.path.join(PATH,'Ipsdca022.cg1'))
#data = reflred.load(os.path.join(PATH,'small.cg1'))
if __name__ == "__main__":
    import pylab
    print "shape=",data.detector.counts.shape
    pylab.imshow(numpy.log(data.detector.counts[1]+1).T)
    #pylab.pcolor(numpy.log(data.detector.counts+1))
    pylab.show()
    pass
