from numpy import zeros_like, complex128, exp, array, empty, sum, newaxis

def greens_form_line(x0, y0, x1, y1, qx, qy):
    if y0 == y1:
        result = zeros_like(qx, dtype=complex128)
    elif x0 == x1:
        # m is infinity...
        result = -1.0 / (qx*qy) * exp(1j * qx * x0)
        result *= (exp(1j * qy * y1) - exp(1j * qy * y0))        
    else:
        m = (y1 - y0)/(x1 - x0)
        result = -1.0 / (qx*(qy + qx/m)) * exp(1j * qx * (x0 - y0/m))
        result *= (exp(1j * (qx/m + qy) * y1) - exp(1j * (qx/m + qy) * y0))
    return result

def greens_form_shape(points, qx, qy):
    result = zeros_like(qx, dtype=complex128)
    numpoints = len(points)
    for i in range(numpoints):
        x0,y0 = points[i]
        x1,y1 = points[(i+1) % numpoints] # loops back to zero for last point.
        result += greens_form_line(x0, y0, x1, y1, qx, qy)
    return result
    
def div_form_line(x0, y0, x1, y1, qx, qy):
    dx = x1 - x0
    dy = y1 - y0
    qx = qx[:,newaxis]
    qy = qy[:,newaxis]
    qlensq = qx**2 + qy**2
    result = 1.0/qlensq * (-qx*dy + qy*dx) / (qx*dx + qy*dy)
    result *= exp(1j*(qx*x1 + qy*y1)) - exp(1j*(qx*x0 + qy*y0))
    return result

def div_form_shape(points, qx, qy):
    result = zeros_like(qx, dtype=complex128)
    if len(points)>0:
        arr_points = array(points + [points[0],])
        x0 = arr_points[:-1,0]
        y0 = arr_points[:-1,1]
        x1 = arr_points[1:,0]
        y1 = arr_points[1:,1]
        subresult = div_form_line(x0, y0, x1, y1, qx, qy)
        result = sum(subresult, axis=1)
    return result
