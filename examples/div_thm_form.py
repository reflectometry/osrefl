from numpy import zeros_like, complex128, exp, array, empty, sum, newaxis, zeros, logical_and, logical_not, size

def greens_form_line(x0, y0, x1, y1, qx, qy):
    if y0 == y1:
        result = 0.0
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
    qxl = qx[..., newaxis] 
    qyl = qy[..., newaxis]
    #qx_z = qx[zero_qs][..., newaxis] 
    #qy_z = qy[zero_qs][..., newaxis] 
    dx = x1 - x0
    dy = y1 - y0
    
    
    qlensq = qx**2 + qy**2
    result = zeros(qlensq.shape + x0.shape, dtype=complex128)
    zero_qs = (qlensq == 0.0)
    qlensq.shape += (1,) # add axis for spatial component
    nonzero_qs = logical_not(zero_qs)
    print qlensq[nonzero_qs].shape
    result[nonzero_qs] = complex128(1.0)/(qlensq[nonzero_qs]) * (-qxl*dy + qyl*dx)[nonzero_qs] / ((qxl*dx + qyl*dy)[nonzero_qs])
    result[nonzero_qs] *= (exp(1j*(qxl*x1 + qyl*y1)) - exp(1j*(qxl*x0 + qyl*y0)))[nonzero_qs]
    result[zero_qs, :] = (x0 * dy - y0 * dx)
    return result

def div_form_line_old(x0, y0, x1, y1, qx, qy):
    qxl = qx[..., newaxis] 
    qyl = qy[..., newaxis] 
#    qxl = qx[newaxis,:,newaxis] # put qx 2nd-last axis
#    qyl = qy[newaxis,newaxis,:] # put qy as last axis
#    x0l = x0[:,newaxis,newaxis] # spatial info on first axis
#    x1l = x1[:,newaxis,newaxis]
#    y0l = y0[:,newaxis,newaxis]
#    y1l = y1[:,newaxis,newaxis]
#    dxl = x1l - x0l
#    dyl = y1l - y0l
    dx = x1 - x0
    dy = y1 - y0


    qlensq = qxl**2 + qyl**2
    result = complex128(1.0)/qlensq * (-qxl*dy + qyl*dx) / (qxl*dx + qyl*dy)
    result *= (exp(1j*(qxl*x1 + qyl*y1)) - exp(1j*(qxl*x0 + qyl*y0)))
    return result

def div_form_shape(points, qx, qy):
    result = 0.0
    if len(points)>0:
        arr_points = array(points + [points[0],])
        x0 = arr_points[:-1,0]
        y0 = arr_points[:-1,1]
        x1 = arr_points[1:,0]
        y1 = arr_points[1:,1]
        subresult = div_form_line(x0, y0, x1, y1, qx, qy)
        result = sum(subresult, axis=(subresult.ndim-1))
    return result

def greens_form_shape_array(points, qx, qy):
    #qxl = qx[:,newaxis] # put qx 2nd-last axis
    #qyl = qy[newaxis,:] # put qy as last axis
    result = zeros((qx.shape[0], qy.shape[0]), dtype=complex128)
    numpoints = len(points)
    for i in range(numpoints):
        x0,y0 = points[i]
        x1,y1 = points[(i+1) % numpoints] # loops back to zero for last point.
        result += greens_form_line(x0, y0, x1, y1, qx, qy)
    return result
    
    
