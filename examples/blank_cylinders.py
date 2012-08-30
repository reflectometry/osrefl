from greens_thm_form import greens_form_line, greens_form_shape
from numpy import arange, linspace, float64, indices, zeros_like, ones_like, pi, sin, complex128, array, exp, newaxis, cumsum, sum, cos, sin, log, log10
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
    rect.area = dx * dy
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

def arc(r, theta_start, theta_end, x_center, y_center, theta_step=1.0, close=True, sld=0.0, sldi=0.0, ):
    a = shape('arc')
    a.theta_start = theta_start
    a.theta_end = theta_end
    a.area = pi * r**2 * abs(theta_end - theta_start)/360.0
    if close == True:
        a.points.append([x_center, y_center]) # center point
    numpoints = (theta_end - theta_start) / theta_step + 1
    thetas = linspace(theta_start, theta_end, numpoints) * pi/180 # to radians
    for th in thetas:
        a.points.append([r*cos(th) + x_center, r*sin(th) + y_center])
    a.sld = sld
    a.sldi = sldi
    return a


def limit_cyl(arc, xmin=0.0, xmax=0.0, ymin=0.0, ymax=0.0):
    new_arc = shape('arc')
    new_arc.sld = arc.sld
    new_arc.sldi = arc.sldi
    new_arc.theta_start = arc.theta_start
    new_arc.theta_end = arc.theta_end
    #new_arc.area = arc.area
    
    for point in arc.points:
        if (point[0] >= xmin) and (point[0] <= xmax) and (point[1] >=ymin) and (point[1] <= ymax):
            new_arc.points.append(point)
    
    if len(new_arc.points) < 3:
        new_arc.area = 0.0
    else:
        new_arc.area = (len(new_arc.points) - 2) / 360.0 * arc.area
    
    return new_arc
    
def conj(sld):
    conjugate_sld = sld.copy()
    conjugate_sld[:,2] *= -1
    return conjugate_sld
     
# alternating SLD
wavelength = 1.24 # x-ray wavelength, Angstroms

spacing = 600.0 # distance between cylinder centers
radius = 200.0 # Angstroms, radius of cylinders
thickness = 300.0 # Angstrom, thickness of cylinder layer
sublayer_thickness = 200.0 # Angstrom, full layer of matrix below cylinders
matrix_sld = pi/(wavelength**2) * 2.0 * 1.0e-6 # substrate
matrix_sldi = pi/(wavelength**2) * 2.0 * 1.0e-7 # absorption in substrate
cyl_sld = 0.0
cyl_sldi = 0.0 # cylinders are holes in matrix

unit_dx = 2.0 * spacing
unit_dy = 1.0 * spacing

matrix = rectangle(0,0, 3000, 3000, matrix_sld, matrix_sldi)

cylinders = []
centers = []

for i in range(3):
    for j in range(6):
        x0 = i * 2.0 * spacing
        y0 = j * spacing
        x1 = x0 + spacing # basis 
        y1 = y0 + spacing/2.0
        cylinders.append(arc(radius, 0.0, 360.0, x0, y0, sld=cyl_sld, sldi=cyl_sldi))
        cylinders.append(arc(radius, 0.0, 360.0, x1, y1, sld=cyl_sld, sldi=cyl_sldi))

cyl_area = 0.0
for cyl in cylinders:
    cyl_area += cyl.area

clipped_cylinders = [limit_cyl(cyl, xmin=0.0, xmax=3000.0, ymin=0.0, ymax=3000.0) for cyl in cylinders]
clipped_cyl_area = 0.0
for cyl in clipped_cylinders:
    clipped_cyl_area += cyl.area

print "clipped_cyl_area / matrix.area = ", clipped_cyl_area / matrix.area
print "ratio should be 0.3491 for FCT planar array with a/b = 2 and r = a/6"

avg_sld = (matrix.area * matrix_sld + clipped_cyl_area * cyl_sld) / matrix.area
avg_sldi = (matrix.area * matrix_sldi + clipped_cyl_area * cyl_sldi) / matrix.area

front_sld = 0.0 # air
back_sld = pi/(wavelength**2) * 2.0 * 5.0e-6 # substrate
back_sldi = pi/(wavelength**2) * 2.0 * 7.0e-8 # absorption in substrate

qz = linspace(0.01, 0.21, 501)
qy = linspace(-0.1, 0.1, 500)
qx = ones_like(qy, dtype=complex128) * 1e-8

SLDArray =  [ [0,0,0], # air
              [avg_sld, thickness, avg_sldi], # sample
              [matrix_sld, sublayer_thickness, matrix_sldi], # full matrix layer under cylinders
              [back_sld, 0, back_sldi] ]

FT = zeros_like(qx, dtype=complex128)
for cyl in clipped_cylinders:
    FT += greens_form_shape(cyl.points, qx, qy) * (cyl.sld)

FT += greens_form_shape(matrix.points, qx, qy) * (matrix.sld)
FT += greens_form_shape(matrix.points, qx, qy) * (-avg_sld)

SLDArray = array(SLDArray)

def calc_gisans(alpha_in, show_plot=True):

    #alpha_in = 0.25 # incoming beam angle

    kz_in = 2*pi/wavelength * sin(alpha_in * pi/180.0)
    kz_out = kz_in - qz

    wf_in = dwbaWavefunction(kz_in, SLDArray)
    wf_out = dwbaWavefunction(-kz_out, conj(SLDArray))

    zs = cumsum(SLDArray[1:-1,1])
    dz = SLDArray[1:-1,1][:,newaxis]
    z_array = array(zs)[:,newaxis]

    qz_inside = wf_in.kz_l[1] - wf_out.kz_l[1]
    
    # the overlap is the forward-moving amplitude c in psi_in multiplied by 
    # the forward-moving amplitude in the time-reversed psi_out, which
    # ends up being the backward-moving amplitude d in the non-time-reversed psi_out
    # (which is calculated by the wavefunction calculator)
    # ... and vice-verso for d and c in psi_in and psi_out
    overlap  = wf_out.d[1] * wf_in.c[1] / (1j * qz_inside) * \
        (exp(1j * qz_inside * thickness) - exp(1j * qz_inside * 0.0)) 
    overlap += wf_out.c[1] * wf_in.d[1] / (-1j * qz_inside) * \
        (exp(-1j * qz_inside * thickness) - exp(-1j * qz_inside * 0.0)) 

    overlap_BA  = 1.0 / (1j * qz) * \
        (exp(1j * qz * thickness) - exp(1j * qz * 0.0)) 
    overlap_BA += 1.0 / (-1j * qz_inside) * \
        (exp(-1j * qz * thickness) - exp(-1j * qz * 0.0)) 


    gisans = overlap[:,newaxis] * FT[newaxis, :]
    gisans_BA = overlap_BA[:,newaxis] * FT[newaxis, :]
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

    

