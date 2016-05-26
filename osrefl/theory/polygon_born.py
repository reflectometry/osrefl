from pylab import *
from numpy import *
from .DWBAGISANS import *
import copy

EPSILON = 1e-7
k0 = 2*pi/5.0 # inverse Angstroms - initial k for MAGIK
rho_Si = 2.07e-6
rho_Ni = 9.408e-6
rho_Ti = -1.925e-6
rho_Cr = 3.027e-6
rho_Au = 4.5e-6

def s(qx, qz, x0, x1, z0, z1):
    o = 1/(qx**2 + qz**2)
    dlx = x1-x0
    dlz = z1-z0
    o *= (qx * dlz + qz * dlx) / (qx * dlx + qz * dlz)
    o*= (exp(1j * (qx * x1 + qz * z1)) - exp(1j * (qx * x0 + qz * z0)))
    return o    


def f(K, qx, qz, rl, rlp1, rm, rmp1, Dx, sigmax, sigmaz, epsilon=EPSILON):
    """ K is the index of the qx Bragg peak, 
    rl is the position (x,z) of the lth vertex
    rlp1 is the position of the l+1 vertex
    rm is the position of the mth vertex
    rmp1 is the position of the m+1 vertex """
    
    qKx = 2*pi*K/Dx
    qK = sqrt(qKx**2 + qz**2)
    #kz_in = k0 * sin(theta_in)
    #kz_out = k0 * sin(theta_out)
    dl = rlp1 - rl
    dm = rmp1 - rm
    tanl = (qKx * dl[1] + qz * dl[0]) / ( qKx * dl[0] + qz * dl[1] - epsilon )
    tanm = (qKx * dm[1] + qz * dm[0]) / ( qKx * dm[0] + qz * dm[1] + epsilon )
    x_prefac = exp(-(qKx - qx)**2 * sigmax**2)
    cos_sum = 0.0
    cos_sum +=  exp(-(rmp1[1] - rlp1[1])**2 / (4 * sigmaz**2)) * cos(qKx*(rmp1[0]-rlp1[0]) + qz*(rmp1[1]-rlp1[1]) + epsilon)
    cos_sum += -exp(-(rm[1] - rlp1[1])**2 / (4 * sigmaz**2)) * cos(qKx*(rm[0]-rlp1[0]) + qz*(rm[1]-rlp1[1]))
    cos_sum +=  exp(-(rm[1] - rl[1])**2 / (4 * sigmaz**2)) * cos(qKx*(rm[0]-rl[0]) + qz*(rm[1]-rl[1]) - epsilon)
    cos_sum += -exp(-(rmp1[1] - rl[1])**2 / (4 * sigmaz**2)) * cos(qKx*(rmp1[0]-rl[0]) + qz*(rmp1[1]-rl[1]))
    result = x_prefac * tanl * tanm * cos_sum / qK**4
    return result
    

K = 1
qx = linspace(-0.00033, 0.00033, 300).reshape(300,1)
qz = linspace(0.001, 0.25, 151).reshape(1, 151)

inplane_angle = 0;
stripe_repeat = 10e4 # 10 um
Dx = stripe_repeat / cos(inplane_angle) # effective spacing for probe
sigmax = 1000e4 # 1 mm
sigmaz = 10000.0



def projected_coh(theta, a, b):
    # in the limit a>>b (a is the coherence along the beam direction, b is the transverse
    #  coherence), then the projected length c is  b * sqrt((tan(theta)^2 + 1)/(tan(theta)^2 + b^2 / a^2))
    # For small values of tan(theta) which we are encountering in reflectivity, this becomes 
    # c approx. b / sqrt(tan(theta)^2 + b^2 / a^2)
    return a * b * sqrt( (((tan(theta))**2 + 1 ) / (a**2 * (tan(theta))**2 + b**2)) )
    


def fd(rs, rhos, K=0, epsilon=EPSILON, Dx=Dx, sigmax=sigmax, sigmaz=sigmaz, qx=qx, qz=qz):
    # rs is a list of lists: each represents a shape
    # rhos is a list of lists that matches the shape of r: the rho for each shape
    # the r is counterclockwise, so for a substrate it is right-to-left, then for a stripe
    # on top it is left-to-right at the same interface.
    result = 0
    qsqr = qx**2 + qz**2
    thetaQ = arcsin(sqrt(qsqr) / (2*k0))
    #tilt = arctan2(qKx, qz)
    tilt = arctan2(qx, qz)
    theta_in = thetaQ - tilt
    theta_out = thetaQ + tilt
    coherence_x = projected_coh(theta_in, sigmax, sigmaz)
    #coherence_x = sigmax
    segments = []
    rho = []
    for r, rh in zip(rs, rhos):
        for i in range(len(r)-1):
            segments.append([r[i], r[i+1]])
            rho.append(rh[i])
    
    print len(segments), len(rho)
    for i in range(len(segments) - 1):
        s = segments[i]
        result += rho[i]**2 * 0.5 * f(K, qx, qz, s[0], s[1], s[0], s[1], Dx, coherence_x, sigmaz, epsilon=epsilon)
        #print "(%d,%d)" % (i,i)
        for j in range(i+1, len(segments)-1):
            #print "(%d,%d)" % (i,j)
            other = segments[j]
            result += rho[i] * rho[j] * f(K, qx, qz, s[0], s[1], other[0], other[1], Dx, coherence_x, sigmaz, epsilon=epsilon)              
    return result * 4 * pi**2 / Dx**2 / (k0**2 * abs(sin(theta_out)) * abs(sin(theta_in)) - 4.0 * pi * rho_Si)

extent = [qx.min(), qx.max(), qz.min(), qz.max()]
layer_thickness = 100 # angstroms

def box(min_x, max_x, thickness, z0):
    """ counterclockwise closed box """
    return array([[max_x, z0], [max_x, z0+thickness], [min_x, z0+thickness], [min_x, z0], [max_x, z0]])

def stack(rhos, thicknesses, min_x, max_x, z0):
    """ lists start at the substrate side """
    rho_stack = []
    r_stack = []
    z = z0
    for r,t in zip(rhos, thicknesses):
        r_stack.append(box(min_x, max_x, t, z))
        rho_stack.append(array([r,]*4))
        z += t
    return rho_stack, r_stack 
