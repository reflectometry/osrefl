from numpy import zeros_like, complex128, exp

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
