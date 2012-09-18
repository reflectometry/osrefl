from numpy import arange, linspace, float64, indices, zeros_like, ones_like, pi, sin, complex128, array, exp, newaxis, cumsum, sum, log10, complex64
from GISANS_problem import Shape, GISANS_problem
from OFFSPEC_problem import OFFSPEC_problem

def rectangle(x0, y0, dx, dy, sld=0.0, sldi=0.0):
    #generate points for a rectangle
    rect = Shape('rectangle')
    rect.points = [[x0,y0], [x0+dx, y0], [x0+dx, y0+dy], [x0, y0+dy]]
    rect.sld = sld
    rect.sldi = sldi
    return rect

# rectangles for inplane stripes: have width = 25 nm with
# alternating SLD
 
wavelength = 1.24 # x-ray wavelength, Angstroms
 
delta = [1.0e-6, 3.0e-6] * 6
beta = [1.0e-7, 3.0e-7] * 6
width = [300, 200] * 6 # Angstroms

y0 = cumsum(array(width))

sldn = pi/(wavelength**2) * 2.0 * array(delta)
avg_sldn = sum(sldn*array(width))/sum(array(width))
sldi = pi/(wavelength**2) * 2.0 * array(beta)
avg_sldi = sum(sldi*array(width))/sum(array(width))

rects = [rectangle(0, y0[i], 3000, width[i], sldn[i], sldi[i]) for i in range(12)]

thickness = 500.0 # Angstrom, thickness of layer
front_sld = 0.0 # air
back_sld = pi/(wavelength**2) * 2.0 * 5.0e-6 # substrate
back_sldi = pi/(wavelength**2) * 2.0 * 7.0e-8 # absorption in substrate
dz = thickness

qz = linspace(0.01, 0.11, 501)
qy = linspace(-0.1, 0.1, 500)
qx = array([1e-10], dtype=complex128)
#qx = ones_like(qy, dtype=complex128) * 1e-10

sublayers = [[rects, avg_sldn, avg_sldi, thickness] ]
matrix = rectangle(0,0, 3000, 3000, 0.0, 0.0) # empty matrix

g_problem = GISANS_problem(sublayers, matrix, front_sld, 0.0, back_sld, back_sldi, wavelength, qx, qy, qz)

oqz = linspace(0.03, 0.21, 501)
oqy = array([1e-10], dtype=complex128)
oqx = linspace(-2e-3+1e-10, 2e-3, 501)

orects = [rectangle(y0[i]*100.0, 0, width[i]*100.0, 300000, sldn[i], sldi[i]) for i in range(12)] # along x!
omatrix = rectangle(0,0, 300000, 300000, 0.0, 0.0) # empty matrix
osublayers = [[orects, avg_sldn, avg_sldi, thickness] ]

o_problem = OFFSPEC_problem(osublayers, omatrix, front_sld, 0.0, back_sld, 0.0, wavelength, oqx, oqy, oqz)
