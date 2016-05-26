from pylab import *
from numpy import *
import copy
from osrefl.theory.polygon_born import *

EPSILON = 1e-7
k0 = 2*pi/5.0 # inverse Angstroms - initial k for MAGIK
rho_Si = 2.07e-6
rho_Ni = 9.408e-6
rho_Ti = -1.925e-6
rho_Cr = 3.027e-6
rho_Au = 4.5e-6

K = 1
qx = linspace(-0.00033, 0.00033, 300).reshape(300,1)
qz = linspace(0.001, 0.25, 151).reshape(1, 151)

inplane_angle = 0;
stripe_repeat = 10e4 # 10 um
Dx = stripe_repeat / cos(inplane_angle) # effective spacing for probe
sigmax = 1000e4 # 1 mm
sigmaz = 10000.0

extent = [qx.min(), qx.max(), qz.min(), qz.max()]
layer_thickness = 100 # angstroms

def canceling_multilayer_R(Dx, sticklayer_thickness, offset_thickness, layer_thickness, caplayer_thickness, sticklayer_rho, offset_rho, rho1, rho2, substrate_rho, caplayer_rho, FractionSample1=0.54, N_bilayers=5):
    # assumes layer1 and layer2 in bilayer are same thickness (same for sticking layers 1 and 2)
    r_substrate = array([[Dx, 0], [0,0]])

    min_x1 = 0
    max_x1 = Dx*FractionSample1
    rhos1 = [sticklayer_rho, offset_rho, sticklayer_rho]
    rhos1.extend([rho1, rho2] * N_bilayers)
    rhos1.extend([caplayer_rho,])
    thicknesses1 = [sticklayer_thickness, offset_thickness, sticklayer_thickness]
    thicknesses1 += [layer_thickness,] * N_bilayers * 2
    thicknesses1 += [caplayer_thickness,]
    rho_stack1, r_stack1 = stack(rhos1, thicknesses1, min_x1, max_x1, 0)
    rho_stack1 += [array([substrate_rho,]),]
    r_stack1 += [r_substrate,]
    
    min_x2 = Dx*FractionSample1
    max_x2 = Dx
    rhos2 = [sticklayer_rho,]
    rhos2.extend([rho1, rho2] * N_bilayers)
    rhos2.extend([caplayer_rho,])
    thicknesses2 = [sticklayer_thickness,]
    thicknesses2.extend([layer_thickness,] * N_bilayers * 2)
    thicknesses2.extend([caplayer_thickness,])
    rho_stack2, r_stack2 = stack(rhos2, thicknesses2, min_x2, max_x2, 0)
    rho_stack2.append(array([substrate_rho,]))
    r_stack2.append(r_substrate)

    R_cancel = sum(array([fd(r_stack1 + r_stack2, rho_stack1 + rho_stack2, K, Dx=Dx) for K in range(-5,5,1)]), axis=0)
    return R_cancel
    
figure()
inplane_angle = 0
Dx = stripe_repeat / cos(inplane_angle) # effective spacing for probe
R_cancel_5bilayers_0deg = canceling_multilayer_R(Dx, 56.0, 60.0, 147.0, 199.8, rho_Cr, rho_Au, rho_Au, rho_Ni, rho_Si, rho_Au, FractionSample1 = 0.54, N_bilayers=5)
plot(qz[0], abs(sum(R_cancel_5bilayers_0deg[145:155, :], axis=0)), label="5 bilayers + inverse Ni/Cr at %.0f deg (parallel)" % (degrees(inplane_angle),))

#inplane_angle = pi/4.0
#Dx = stripe_repeat / cos(inplane_angle) # effective spacing for probe
#R_cancel_5bilayers_45deg = canceling_multilayer_R(Dx, layer_thickness, rho_Ni, rho_Cr, rho_Si, N_bilayers=5)
#plot(qz[0], abs(sum(R_cancel_5bilayers_45deg[145:155, :], axis=0)), label="5 bilayers + inverse Ni/Cr at %.0f deg (parallel)" % (degrees(inplane_angle),))


# NOTE: the intensity at Qx=0 for the "forbidden" peaks seems to come from the fuzzing of the Qx!=0 peaks...
# this makes physical sense in some ways, but is confusing in realspace.
figure()
imshow(log10(abs(R_cancel_5bilayers_0deg)+1e-7).T, origin='lower', aspect='auto', extent=extent)
#figure()
#imshow(log10(abs(R_cancel_5bilayers_45deg)+1e-7).T, origin='lower', aspect='auto', extent=extent)

show()

#exact solution:


