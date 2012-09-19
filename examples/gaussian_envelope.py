from numpy import *
#from pylab import *

def normgauss(x, sigma, x0=0):
    return 1.0/sqrt(2*pi*sigma**2) * exp(-(x-x0)**2 / (2.0*sigma**2))

def FWHM_to_sigma(FWHM):
    # from FWHM to sigma parameter:
    return FWHM / sqrt(8.0 * log(2.0))

def rotxz(xp, zp, theta_deg):
    theta = theta_deg * pi / 180.0
    x = xp*cos(theta) - zp*sin(theta)
    z = xp*sin(theta) + zp*cos(theta)
    return x, z
    
def gengauss(sigma, x0=0, threshold=0.001, npts = 200):
    """ when gaussian drops below 0.001 of max amplitude, cut off x """
    cutoff = sqrt( -2.0 * sigma**2 * log( threshold ) )
    x = linspace( x0 - cutoff, x0 + cutoff, npts )
    y = normgauss( x, sigma, x0 )
    return x, y

def plot_intersection():
    pass

def gengaussFWHM(FWHM, x0=0, threshold=0.001, npts = 200):
    sigma = FWHM_to_sigma(FWHM)
    return gengauss(sigma, x0, threshold, npts)


        
#y_offset + amplitude * exp( - ( x - center )**2 * 4 * log(2) / FWHM**2 
