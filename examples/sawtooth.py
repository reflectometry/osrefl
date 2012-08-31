from greens_thm_form import greens_form_line, greens_form_shape
from numpy import arange, linspace, float64, indices, zeros_like, ones_like, pi, sin, complex128, array, exp, newaxis, cumsum, sum, log10
import itertools

class shape:
    def __init__(self, name):
        self.name = name
        self.points = []
        self.sld = 0.0
        self.sldi = 0.0

def rectangle(x0, y0, dx, dy, sld=0.0, sldi=0.0):
    #generate points for a rectangle
    rect = shape('rectangle')
    rect.points = [[x0,y0], [x0+dx, y0], [x0+dx, y0+dy], [x0, y0+dy]]
    rect.sld = sld
    rect.sldi = sldi
    return rect

def sawtooth(z, n=6, x_length=3000.0, base_width=500.0, height=300.0,  sld=0.0, sldi=0.0, sld_front=0.0, sldi_front=0.0):
    if z>height:
        return [], sld_front
    
    width = (z / height) * base_width
    front_width = base_width - width
    rects = [rectangle(0, base_width*(i+0.5) - width/2.0, x_length, width, sld, sldi) for i in range(n)]
    # now rectangles for the gaps between the sawtooths...
    if (sld_front !=0.0 and sldi_front != 0.0):
        front_rects = [rectangle(0, 0, x_length, front_width/2.0, sld_front, sldi_front)]
        front_rects.extend([rectangle(0, base_width*(i+0.5)+width/2.0, x_length, front_width, sld_front, sldi_front) for i in range(1,n-1)])
        front_rects.append(rectangle(0, base_width*(n-0.5)+width/2.0, x_length, front_width/2.0, sld_front, sldi_front))
        rects.extend(front_rects)
    
    # now calculate the average SLD (nuclear) for the layer
    avg_sld = (width * sld + front_width * sld_front) / base_width
    avg_sldi = (width * sldi + front_width * sldi_front) / base_width
    return rects, avg_sld, avg_sldi


def conj(sld):
    conjugate_sld = sld.copy()
    conjugate_sld[:,2] *= -1
    return conjugate_sld
    
# rectangles for inplane stripes: have width = 25 nm with
# alternating SLD
 
wavelength = 1.24 # x-ray wavelength, Angstroms

front_sld = 0.0 # air
back_sld = pi/(wavelength**2) * 2.0 * 5.0e-6 # substrate
back_sldi = pi/(wavelength**2) * 2.0 * 7.0e-8 # absorption in substrate


delta = [1.0e-6, 3.0e-6] * 6
beta = [1.0e-7, 3.0e-7] * 6
width = [300, 200] * 6 # Angstroms

y0 = cumsum(array(width))
#arange(12.0) * 250

#Qx, Qy = (indices((501, 501), dtype=float64) -249.99)* 1.0 / 2500.0 # range is 0->1 nm^{-1} for each of Qx, Qy
qz = linspace(0.01, 0.41, 501)
qy = linspace(-0.1, 0.1, 500)
qx = ones_like(qy, dtype=complex128) * 1e-8

zs = linspace(0.0, 300.0, 31)
dz = zs[1] - zs[0]

sublayers = [sawtooth(z, sld=back_sld, sldi=back_sldi) for z in zs]
matrix = rectangle(0,0, 3000, 3000, front_sld, 0.0)
FTs = []
SLDArray =  [ [0,0,0], ]# air]
for sl in sublayers:
    FT = zeros_like(qx, dtype=complex128)
    for rect in sl[0]:
        FT += greens_form_shape(rect.points, qx, qy) * (rect.sld)
    FT += greens_form_shape(matrix.points, qx, qy) * (-sl[1]) # subtract FT of average SLD
    FTs.append(FT)
    SLDArray.append([sl[1], dz, sl[2]])

SLDArray.append([back_sld, 0, back_sldi])
SLDArray = array(SLDArray)

#alpha_in = 0.18 # incoming beam angle



thickness = 500.0 # Angstrom, thickness of layer

from osrefl.theory.DWBAGISANS import dwbaWavefunction

#qy, qz = indices((501, 501), dtype=float64)
#qy = (qy -249.99) / 250.0
#qz = qz / 500.0 + 0.1
def calc_gisans(alpha_in, show_plot=True):
    kz_in_0 = 2*pi/wavelength * sin(alpha_in * pi/180.0)
    kz_out_0 = kz_in_0 - qz

    wf_in = dwbaWavefunction(kz_in_0, SLDArray)
    wf_out = dwbaWavefunction(-kz_out_0, conj(SLDArray))
    
    kz_in_l = wf_in.kz_l # inside the layers
    kz_out_l = -wf_out.kz_l # inside the layers

    zs = cumsum(SLDArray[1:-1,1])
    dz = zs[1] - zs[0]
    #dz = SLDArray[1:-1,1][:,newaxis]
    z_array = array(zs)[:,newaxis]

    qrt_inside =  kz_in_l[1:-1] - kz_out_l[1:-1]
    qtt_inside =  kz_in_l[1:-1] + kz_out_l[1:-1]
    qtr_inside = -kz_in_l[1:-1] + kz_out_l[1:-1]
    qrr_inside = -kz_in_l[1:-1] - kz_out_l[1:-1]
    
    
    # the overlap is the forward-moving amplitude c in psi_in multiplied by 
    # the forward-moving amplitude in the time-reversed psi_out, which
    # ends up being the backward-moving amplitude d in the non-time-reversed psi_out
    # (which is calculated by the wavefunction calculator)
    # ... and vice-verso for d and c in psi_in and psi_out
    overlap  = wf_out.d[1:-1] * wf_in.c[1:-1] / (1j * qtt_inside) * (exp(1j * qtt_inside * dz) - 1.0)*exp(1j*qtt_inside*z_array)
    overlap += wf_out.c[1:-1] * wf_in.d[1:-1] / (1j * qrr_inside) * (exp(1j * qrr_inside * dz) - 1.0)*exp(1j*qrr_inside*z_array)
    overlap += wf_out.d[1:-1] * wf_in.d[1:-1] / (1j * qtr_inside) * (exp(1j * qtr_inside * dz) - 1.0)*exp(1j*qtr_inside*z_array)
    overlap += wf_out.c[1:-1] * wf_in.c[1:-1] / (1j * qrt_inside) * (exp(1j * qrt_inside * dz) - 1.0)*exp(1j*qrt_inside*z_array)

    overlap_BA  = 1.0 / (1j * qz) * (exp(1j * qz * dz) - 1.0) * exp(1j*qz*z_array)
    #overlap_BA += 1.0 / (-1j * qz) * (exp(-1j * qz * dz) - 1.0) * exp(-1j*qz*z_array)
    #######################################################################3

#    kz_in = 2*pi/wavelength * sin(alpha_in * pi/180.0)
#    kz_out = kz_in - qz

#    wf_in = dwbaWavefunction(kz_in, SLDArray)
#    wf_out = dwbaWavefunction(-kz_out, conj(SLDArray))

    #z_array = array(zs)[:,newaxis]
    #qz_inside = wf_in.kz_l[1:-1] - wf_out.kz_l[1:-1]
#    overlap  = array(wf_out.d)[1:-1] * array(wf_in.c)[1:-1] / (1j * qz_inside) * \
#        (exp(1j * qz_inside * dz) - exp(1j * qz_inside * 0.0)) * exp(1j*qz_inside*z_array)
#    overlap += array(wf_out.c)[1:-1] * array(wf_in.d)[1:-1] / (-1j * qz_inside) * \
#        (exp(-1j * qz_inside * dz) - exp(-1j * qz_inside * 0.0)) * exp(-1j*qz_inside*z_array)

#    overlap_BA  = 1.0 / (1j * qz) * \
#        (exp(1j * qz * dz) - exp(1j * qz * 0.0)) * exp(1j*qz*z_array)
#    overlap_BA += 1.0 / (-1j * qz) * \
#        (exp(-1j * qz * dz) - exp(-1j * qz * 0.0)) * exp(-1j*qz*z_array)
     
    #FTx0 = zeros_like(qy, dtype=complex128)
    #for rect in rects:
    #    FTx0 += greens_form_shape(rect['points'], qx, qy) * (rect['sld'] - avg_sldn)

       
    gisans = sum(overlap[:,:,newaxis] * array(FTs)[:,newaxis,:], axis=0)
    gisans_BA = sum(overlap_BA[:,:,newaxis] * array(FTs)[:,newaxis,:], axis=0)
    extent = [qy.min(), qy.max(), qz.min(), qz.max()]
    
    if show_plot == True:
        from pylab import imshow, figure, colorbar
        figure()
        imshow(log10(abs(gisans)**2), origin='lower', extent=extent, aspect='auto')
        colorbar()
        
        figure()
        imshow(log10(abs(gisans_BA)**2), origin='lower', extent=extent, aspect='auto')
        colorbar()
    return gisans, gisans_BA
    

