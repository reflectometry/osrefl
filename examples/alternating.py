from greens_thm_form import greens_form_line, greens_form_shape
from numpy import arange, linspace, float64, indices, zeros_like, ones_like, pi, sin, complex128, array, exp, newaxis, cumsum, sum, log10
from osrefl.theory.DWBAGISANS import dwbaWavefunction

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

def conj(sld):
    conjugate_sld = sld.copy()
    conjugate_sld[:,2] *= -1
    return conjugate_sld

# rectangles for inplane stripes: have width = 25 nm with
# alternating SLD
 
wavelength = 1.24 # x-ray wavelength, Angstroms
 
delta = [1.0e-6, 3.0e-6] * 6
beta = [1.0e-7, 3.0e-7] * 6
width = [300, 200] * 6 # Angstroms

y0 = cumsum(array(width))
#arange(12.0) * 250

sldn = pi/(wavelength**2) * 2.0 * array(delta)
avg_sldn = sum(sldn)/sldn.shape[0]
sldi = pi/(wavelength**2) * 2.0 * array(beta)

rects = [rectangle(0, y0[i], 3000, width[i], sldn[i], sldi[i]) for i in range(12)]

def inplane_FT():
    # just for testing.
    Qx, Qy = (indices((501, 501), dtype=float64) -249.99)* 1.0 / 2500.0 # range is 0->1 nm^{-1} for each of Qx, Qy
    FT = zeros_like(Qx, dtype=complex128)
    
    for rect in rects:
        FT += greens_form_shape(rect.points, Qx, Qy) * (rect.sld - avg_sldn)
    
    return FT


thickness = 500.0 # Angstrom, thickness of layer
front_sld = 0.0 # air
back_sld = pi/(wavelength**2) * 2.0 * 5.0e-6 # substrate
back_sldi = pi/(wavelength**2) * 2.0 * 7.0e-8 # absorption in substrate

SLDArray =  [ [0,0,0], # air
             [sum(sldn)/len(sldn), thickness, sum(sldi)/len(sldi)], # sample
             [back_sld, 0, back_sldi] ]
SLDArray = array(SLDArray)

matrix = rectangle(0,0, 3000, 3000, front_sld, 0.0)

qz = linspace(0.01, 0.11, 501)
qy = linspace(-0.1, 0.1, 500)
qx = ones_like(qy, dtype=complex128) * 1e-10

FTx0 = zeros_like(qy, dtype=complex128)
for rect in rects:
    FTx0 += greens_form_shape(rect.points, qx, qy) * (rect.sld)
# subtracting the avg sld from a box surrounding the whole layer.
FTx0 += greens_form_shape(matrix.points, qx, qy) * (-avg_sldn)


def calc_gisans(alpha_in, show_plot=True):
    #alpha_in = 0.18 # incoming beam angle
    
    kz_in_0 = 2*pi/wavelength * sin(alpha_in * pi/180.0)
    kz_out_0 = kz_in_0 - qz

    wf_in = dwbaWavefunction(kz_in_0, SLDArray)
    wf_out = dwbaWavefunction(-kz_out_0, conj(SLDArray))
    
    kz_in_l = wf_in.kz_l # inside the layers
    kz_out_l = -wf_out.kz_l # inside the layers

    zs = cumsum(SLDArray[1:-1,1])
    dz = SLDArray[1:-1,1][:,newaxis]
    z_array = array(zs)[:,newaxis]

    qrt_inside =  kz_in_l[1] - kz_out_l[1]
    qtt_inside =  kz_in_l[1] + kz_out_l[1]
    qtr_inside = -kz_in_l[1] + kz_out_l[1]
    qrr_inside = -kz_in_l[1] - kz_out_l[1]
    
    
    # the overlap is the forward-moving amplitude c in psi_in multiplied by 
    # the forward-moving amplitude in the time-reversed psi_out, which
    # ends up being the backward-moving amplitude d in the non-time-reversed psi_out
    # (which is calculated by the wavefunction calculator)
    # ... and vice-verso for d and c in psi_in and psi_out
    overlap  = wf_out.d[1] * wf_in.c[1] / (1j * qtt_inside) * (exp(1j * qtt_inside * thickness) - 1.0)
    overlap += wf_out.c[1] * wf_in.d[1] / (1j * qrr_inside) * (exp(1j * qrr_inside * thickness) - 1.0)
    overlap += wf_out.d[1] * wf_in.d[1] / (1j * qtr_inside) * (exp(1j * qtr_inside * thickness) - 1.0)
    overlap += wf_out.c[1] * wf_in.c[1] / (1j * qrt_inside) * (exp(1j * qrt_inside * thickness) - 1.0)

    overlap_BA  = 1.0 / (1j * qz) * (exp(1j * qz * thickness) - 1.0) 
    #overlap_BA += 1.0 / (-1j * qz) * (exp(-1j * qz * thickness) - 1.0)
    
#    kz_in = 2*pi/wavelength * sin(alpha_in * pi/180.0)
#    kz_out = kz_in - qz

#    wf_in = dwbaWavefunction(kz_in, SLDArray)
#    wf_out = dwbaWavefunction(-kz_out, conj(SLDArray))


#    qz_inside = wf_in.kz_l[1] - wf_out.kz_l[1]
#    overlap  = wf_out.d[1] * wf_in.c[1] / (1j * qz_inside) * \
#        (exp(1j * qz_inside * thickness) - exp(1j * qz_inside * 0.0)) 
#    overlap += wf_out.c[1] * wf_in.d[1] / (-1j * qz_inside) * \
#        (exp(-1j * qz_inside * thickness) - exp(-1j * qz_inside * 0.0)) 

#    overlap_BA  = 1.0 / (1j * qz) * \
#        (exp(1j * qz * thickness) - exp(1j * qz * 0.0)) 
#    overlap_BA += 1.0 / (-1j * qz_inside) * \
#        (exp(-1j * qz * thickness) - exp(-1j * qz * 0.0)) 
       
    gisans = overlap[:,newaxis] * FTx0[newaxis, :]
    gisans_BA = overlap_BA[:,newaxis] * FTx0[newaxis, :]
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

