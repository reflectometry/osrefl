from numpy import arange, linspace, float64, indices, zeros_like, ones_like, pi, sin, complex128, array, exp, newaxis, cumsum, sum, cos, sin, log, log10
from GISANS_problem import Shape, GISANS_problem

def rectangle(x0, y0, dx, dy, sld=0.0, sldi=0.0):
    #generate points for a rectangle
    rect = Shape('rectangle')
    rect.points = [[x0,y0], [x0+dx, y0], [x0+dx, y0+dy], [x0, y0+dy]]
    rect.sld = sld
    rect.sldi = sldi
    rect.area = dx * dy
    return rect

def arc(r, theta_start, theta_end, x_center, y_center, theta_step=1.0, close=True, sld=0.0, sldi=0.0, ):
    a = Shape('arc')
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
    new_arc = Shape('arc')
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
    
     
# alternating SLD
wavelength = 1.24 # x-ray wavelength, Angstroms

spacing = 600.0 # distance between cylinder centers
radius = 200.0 # Angstroms, radius of cylinders
thickness = 300.0 # Angstrom, thickness of cylinder layer
sublayer_thickness = 200.0 # Angstrom, full layer of matrix below cylinders
matrix_sld = pi/(wavelength**2) * 2.0 * 1.0e-6 # substrate
matrix_sldi = pi/(wavelength**2) * 2.0 * 1.0e-7 # absorption in substrate
cyl_sld = -matrix_sld
cyl_sldi = -matrix_sldi # cylinders are holes in matrix

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

def draw_cylinders():
    from pylab import figure, plot
    figure()
    for c in cylinders.clipped_cylinders:
        cp = array(c.points)
        if cp.shape[0] > 0: plot(cp[:,0], cp[:,1])


front_sld = 0.0 # air
back_sld = pi/(wavelength**2) * 2.0 * 5.0e-6 # substrate
back_sldi = pi/(wavelength**2) * 2.0 * 7.0e-8 # absorption in substrate

qz = linspace(0.01, 0.11, 501)
qy = linspace(-0.1, 0.1, 500)
qx = array([1e-8], dtype=complex128)
#qx = ones_like(qy, dtype=complex128) * 1e-8

shapes = clipped_cylinders
dz = thickness
sublayers = [[clipped_cylinders, avg_sld, avg_sldi, thickness], \
             [[], matrix_sld, matrix_sldi, sublayer_thickness]]

problem = GISANS_problem(sublayers, matrix, front_sld, 0.0, back_sld, back_sldi, wavelength, qx, qy, qz)
