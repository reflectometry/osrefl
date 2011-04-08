# Copyright (C) 2008 University of Maryland
# All rights reserved.
# See LICENSE.txt for details.
# Author: Christopher Metting

#Starting Date:5/6/2010


from numpy import *
from pylab import *

def imageload(filename = None):
    
    '''
    Overview:
        Loads a .png image file into a numpy array. In the case that a file name
    is not provided, a system browser will appear and ask the user to select a
    file to load.
 
    
    Parameters:
    
    filename(str) = The path and name of the .png file to be loaded.
    '''
    from pylab import imread,imshow,show
    from numpy import array
    
    if filename == None:
        filename = menuOpen()

    imdata = imread(str(filename))
    return imdata
    
def menuOpen():
    """
    Overview:
        Generic file opening code.
        
    Note:
    -This should be generalized elsewhere so that it may be called by any piece
    of the software.
    """
    import wx, os
    app = wx.PySimpleApp()
    
    dirname = ''
    fileType = "Portable Network Graphics file (.png)|*.png"
    dlg = wx.FileDialog(None, "Choose a file", dirname, "", fileType, wx.OPEN)
    if dlg.ShowModal() == wx.ID_OK:
        filename = dlg.GetFilename()
        dirname = dlg.GetDirectory()
        full_path = os.path.join(dirname, filename)
    dlg.Destroy()
    return full_path


def coursenRes(img,newres):
    '''
    Overview:
        This module is called to coarsen the resolution of a .png image. This 
    is used in the case that the array loaded in from the imread call is to
    large for the needs of the application. It uses a rebinning technique to 
    ensure minimal loss of information.
    
    
    Parameters:
    
    img:(array,[]) = The image array that is being rebinned.
    
    newres:(array,[2]|count) = The new size of the [0] and [1] axis of the 
    image.
    
    '''
    newres = asarray(newres,'float')
    temp = axResRed(img,newres,ax = 0)
    newim =  axResRed(temp,newres,ax = 1)
    
    return newim


def axResRed(img,newres,ax):
    '''
    Overview:
        This is the rebinning algorithm used to downgrade the resolution of a
    .png file. It keeps track of the fractional values as it reads across the
    axis and scales each pixel accordingly.
    
    
    Parameters:
    
    img:(array,[]) = The image array that is being rebinned.
        
    newres:(array,[2]|count) = The new size of the [0] and [1] axis of the 
    image.
    
    ax:(int|count) = The number of the axis which is to be rebinned.
    
    '''
    if ax == 1:
        img =rot90(img)
        newres =[newres[1],newres[0]]
    
    dim = shape(img)

    tempImg = zeros([newres[0],dim[1]])

    xalt = (dim[0]/newres[0])

    
    leftDec = 0.0
    rightDec = xalt-int(xalt)
    lastDecTrac = int(xalt-1)

    tempImg[0,:] = ((sum(img[0:int(xalt-1),:],axis = 0)+(img[xalt,:] *
                                                         rightDec))/xalt)


    ran = arange(1,newres[0]-1)

    for ii in (ran):

        leftDec = 1.0 - rightDec
        
        count = int(xalt-leftDec)
        
        rightDec = (xalt-leftDec) - count

        lastDecTrac = int(lastDecTrac+count+1)

        tempImg[ii,:] = ((img[lastDecTrac,:]*leftDec + 
                     sum(img[lastDecTrac+1:lastDecTrac+1+count,:],axis = 0) +
                    img[lastDecTrac+count+1,:]*rightDec )/xalt)
    
    if rightDec ==0.0:   
        tempImg[-1,:] = ((img[lastDecTrac,:]*leftDec + 
                     sum(img[lastDecTrac+1:lastDecTrac+1+count,:],axis = 0))/
                     xalt)
    else:
        tempImg[-1,:] = ((img[lastDecTrac,:]*leftDec + 
                     sum(img[lastDecTrac+1:lastDecTrac+1+count,:],axis = 0) +
                    img[lastDecTrac+count+1,:]*rightDec )/xalt)
    
    if ax == 1:
        tempImg = rot90(tempImg,-1)
    
    return tempImg


def test():
    import pylab
    import numpy

    test = pylab.imread('/home/mettingc/Documents/sample1_sld.png')[:,:,0]
    new = coursenRes(test,numpy.array([200.0,500.0]))


    figure(0)
    imshow(test)
    colorbar()
    figure(1)
    imshow(new)
    colorbar()
    show()

if __name__=="__main__":test()