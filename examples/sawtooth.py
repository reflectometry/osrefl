from numpy import arange, linspace, float64, indices, zeros_like, ones_like, pi, sin, complex128, array, exp, newaxis, cumsum, sum, log10
from GISANS_problem import Shape, rectangle, GISANS_problem

#def rectangle(x0, y0, dx, dy, sld=0.0, sldi=0.0):
#    #generate points for a rectangle
#    rect = Shape('rectangle')
#    rect.points = [[x0,y0], [x0+dx, y0], [x0+dx, y0+dy], [x0, y0+dy]]
#    rect.sld = sld
#    rect.sldi = sldi
#    return rect

def sawtooth(z, dz, n=6, x_length=3000.0, base_width=500.0, height=300.0,  sld=0.0, sldi=0.0, sld_front=0.0, sldi_front=0.0):
    if z>height:
        return [], sld_front
    
    width = (z / height) * base_width
    front_width = base_width - width
    if width == 0:
        rects = []
    else: 
        rects = [rectangle(0, base_width*(i+0.5) - width/2.0, x_length, width, sld, sldi) for i in range(n)]
        
    #### below is now taken care of with "matrix" rectangle that surrounds every layer.
    # now rectangles for the gaps between the sawtooths...
#    if (sld_front !=0.0 and sldi_front != 0.0):
#        front_rects = [rectangle(0, 0, x_length, front_width/2.0, sld_front, sldi_front)]
#        front_rects.extend([rectangle(0, base_width*(i+0.5)+width/2.0, x_length, front_width, sld_front, sldi_front) for i in range(1,n-1)])
#        front_rects.append(rectangle(0, base_width*(n-0.5)+width/2.0, x_length, front_width/2.0, sld_front, sldi_front))
#        rects.extend(front_rects)
    
    # now calculate the average SLD (nuclear) for the layer
    avg_sld = (width * sld + front_width * sld_front) / base_width
    avg_sldi = (width * sldi + front_width * sldi_front) / base_width
    return rects, avg_sld, avg_sldi, dz

# rectangles for inplane stripes: have width = 25 nm with
# alternating SLD

def draw_planview(shapes, xview = (0,5000), yview=(0,5000)):
    from pylab import plot, figure, draw, Polygon
    slds = [shape.sld for shape in shapes]
    max_sld = max(slds)
    min_sld = min(slds + [0,])
    fig = figure()
    ax = fig.add_subplot(111)
    ax.set_xlim(xview)
    ax.set_ylim(yview)
    draw()
    ps = [Polygon(array(shape.points)) for shape in shapes]
    for p in ps:
        ax.add_patch(p)
    draw()
    
def draw_sideview(sublayers, yview=(0,5000), zview=(-50,400)):
    from pylab import plot, figure, draw, Polygon
    dz = [sl[3] for sl in sublayers]
    thickness = sum(array(dz))
    
    fig = figure()
    ax = fig.add_subplot(111)
    ax.set_xlim(yview)
    ax.set_ylim(zview)
    draw()
    z = 0
    for sl in sublayers:
        for shape in sl[0]:
            sp = array(shape.points)
            ymax = sp[:,1].max()
            ymin = sp[:,1].min()
            sideview = array([[ymin, z],[ymax, z],[ymax, z+sl[3]], [ymin, z+sl[3]]])
            p = Polygon(sideview)
            ax.add_patch(p)
        z += sl[3] # dz
    draw()
 
wavelength = 1.24 # x-ray wavelength, Angstroms
Lx = 3000.
Ly = 3000.

front_sld = 0.0 # air
back_sld = pi/(wavelength**2) * 2.0 * 5.0e-6 # substrate
back_sldi = pi/(wavelength**2) * 2.0 * 7.0e-8 # absorption in substrate


delta = [1.0e-6, 3.0e-6] * 6
beta = [1.0e-7, 3.0e-7] * 6
width = [300, 200] * 6 # Angstroms

y0 = cumsum(array(width))
#arange(12.0) * 250

qz = linspace(0.01, 0.41, 501)[newaxis,newaxis,:]
qy = linspace(-0.1, 0.1, 501)[newaxis,:,newaxis]
qx = array([1e-10])[:,newaxis,newaxis]
#qx = ones_like(qy, dtype=complex128) * 1e-8

zs = linspace(0.0, 300.0, 31)
dz = zs[1] - zs[0]

sublayers = [sawtooth(z, dz, sld=back_sld, sldi=back_sldi) for z in zs]
matrix = rectangle(0,0, Lx, Ly, front_sld, 0.0)

g_problem = GISANS_problem(sublayers, matrix, front_sld, 0.0, back_sld, back_sldi, wavelength, qx, qy, qz, Lx, Ly)

