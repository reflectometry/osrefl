# Copyright (C) 2008 University of Maryland
# All rights reserved.
# See LICENSE.txt for details.
# Author: Christopher Metting

#Starting Date:6/12/2009
'''
File Overview:

    This file holds the approximations used to calculate the form factor.
Although not all calculation components are held here, This is the first file
that the scatter modules go to for their approximations.

'''
from pylab import imshow,colorbar,show,pcolormesh
from numpy import *
import wavefunction_kernel
import osrefl.viewers.view
import czt


def partial_magnetic_BA(struc_cell,mag_cell,Q,lattice,beam):

    '''
    **Overview:**

        This calculation does the magnetic born approximation but assumes that
        the contribution to the magnetic SLD from the qx and qy components of
        the magnetic moment are negligible and the whole system can be estimated
        as only containing magnetic contribution in the qz direction.

        .. note::
            Unfortunately, this method is not accurate for magnetic moments
            aligned in the q directions.



    **Parameters:**

        *struc_cell:* (float:3D array|angstroms^2)
            The structural scattering potential of the feature being scattered
            off of.

        *mag_cell:* (float:3D array|angstroms^2)
            The magnetic scattering potential of the feature being scattered off
            of.

        *Q:* (q_space)
            A Q_space object that holds all of the information about the desired
            q space output.

        *lattice:* (Lattice)
            A lattice object that holds all of the information needed to solve
            the structure factor of the scattering.

        *beam:* (Beam)
            Holds all of the information about the experimental beam needed to
            apply beam dependent corrections to the data.

        '''

    rho_NSF, rho_SF = partial_rho(struc_cell.mag_unit,struc_cell.mag_theta,
                                  struc_cell.mag_phi)

    form_factor = [None,None,None,None]
    intensity = [None,None,None,None]

    effective_rho = struc_cell.unit + rho_NSF

    form_factor[0] = BA_FT(effective_rho,struc_cell.step,Q)

    effective_rho = struc_cell.unit - rho_NSF

    form_factor[1] = BA_FT(effective_rho,struc_cell.step,Q)

    form_factor[2] = BA_FT(rho_SF,struc_cell.step,Q)
    form_factor[3] = BA_FT(rho_SF,struc_cell.step,Q)

    structure_factor = lattice.struc_calc(Q)

    for i in range(4):

        form_factor[i] = complete_formula(form_factor[i],struc_cell.step,
                                          Q.q_list)

        intensity[i] = abs(structure_factor)**2 * abs(form_factor[i])**2

        intensity[i] = (normalize(intensity[i], struc_cell, Q, lattice)).real

    return intensity

def partial_magnetic_BA_long(struc_cell,mag_cell,Q,lattice,beam):
    '''
    This approximation employs the same math as the partial_magnetic_BA method
    however instead of using the czt to calculate the form, it uses the long
    form solution of the exponential.
    '''
    form_factor = [zeros(Q.points,dtype='complex'),
                   zeros(Q.points,dtype='complex'),
                   zeros(Q.points,dtype='complex'),
                   zeros(Q.points,dtype='complex')]

    intensity = [None,None,None,None]
    rho_NSF, rho_SF = partial_rho(struc_cell.mag_unit,struc_cell.mag_theta,
                                  struc_cell.mag_phi)

    effective_rho = [None,None,None,None]
    effective_rho[0] = struc_cell.unit + rho_NSF
    effective_rho[1] = struc_cell.unit - rho_NSF
    effective_rho[2] = rho_SF
    effective_rho[3] = rho_SF

    for i in range (size(Q.q_list[2])):
        print i
        for ii in range (size(Q.q_list[1])):
            for iii in range (size(Q.q_list[0])):
                arg = [Q.q_list[2][i], Q, effective_rho,struc_cell]
                form_factor[0][iii,ii,i],form_factor[1][iii,ii,i],
                form_factor[2][iii,ii,i], form_factor[3][iii,ii,i] = (
                          part_mag_func(Q.q_list[0][iii],Q.q_list[1][ii],arg))

    structure_factor = lattice.struc_calc(Q)

    for i in range(4):
        form_factor[i] = complete_formula(form_factor[i],struc_cell.step,
                                          Q.q_list)


        intensity[i] = abs(structure_factor)**2 * abs(form_factor[i])**2
        intensity[i] = (normalize(intensity[i], struc_cell, Q, lattice)).real

    return intensity

def cudaMagBA(struc_cell, Q, lattice, beam, omf, precision = 'float32',
              refract = False):

    from threadSMBA import magCudaBA_form
    intensity = [None]*4

    form_factor = magCudaBA_form(struc_cell,Q,lattice,beam,omf,
                              precision = precision,refract = refract)

    #structure_factor = lattice.struc_calc(Q)
    structure_factor = lattice.gauss_struc_calc(Q)

    for i in range(4):
        #intensity[i] = sum(abs(form_factor[i])**2,axis = 1)
        intensity[i] = (structure_factor)**2 * abs(form_factor[i])**2
        intensity[i] = normalize(intensity[i],struc_cell, Q, lattice).real

    return intensity

def magneticBA(struc_cell,Q,lattice,beam, omf):
    '''
    **Overview:**

        This calculation solves the Born Approximation for a magnetic sample
        using an :class:`~omfLoader.Omf` object which hold the information about
        the magnetic moments in a sample.

    **Parameters:**

        *struc_cell:* (float:3D array|angstroms^2)
            The structural scattering potential of the feature being scattered
            off of.

        *Q:* (q_space)
            A Q_space object that holds all of the information about the desired
            q space output.

        *lattice:* (Lattice)
            A lattice object that holds all of the information needed to solve
            the structure factor of the scattering.

        *beam:* (Beam)
            Holds all of the information about the experimental beam needed to
            apply beam dependent corrections to the data.

        *omf:* (Omf)
            This is an object which holds the magnetic moment information about
            the sample. It contains three arrays of the the same size as the
            unit cell which hold each of the x, y, and z components of the
            magnetic moment.

    '''


    form_factor = [zeros(Q.points,dtype='complex'),
                   zeros(Q.points,dtype='complex'),
                   zeros(Q.points,dtype='complex'),
                   zeros(Q.points,dtype='complex')]

    intensity = [None]*4

    qxn,qyn,qzn = Q.normalize()
    rho_m = omf.ConvertRho()

    for i in range (size(Q.q_list[0])):
        print i
        for ii in range (size(Q.q_list[1])):
            for iii in range (size(Q.q_list[2])):

                q = [Q.q_list[0][i],Q.q_list[1][ii],Q.q_list[2][iii]]
                qn = [qxn[i,0,0],qyn[0,ii,0],qzn[0,0,iii]]

                (form_factor[0][i,ii,iii],form_factor[1][i,ii,iii],
                 form_factor[2][i,ii,iii], form_factor[3][i,ii,iii]) = mag_func(
                                        struc_cell, omf, q, qn, rho_m)

    structure_factor = lattice.gauss_struc_calc(Q)

    for i in range(4):
        form_factor[i] = complete_formula(form_factor[i],struc_cell.step,
                                         Q.q_list)
        #intensity[i] = sum(abs(form_factor[i])**2,axis = 1)
        intensity[i] = (structure_factor)**2 * abs(form_factor[i])**2
        intensity[i] = normalize(intensity[i],struc_cell, Q, lattice).real

    return intensity


def mag_func(strucUnit,omf, q, qn, magCell):

        recip = [None,None,None,None]
        form_factor = [None,None,None,None]

        QdotM = (qn[0] * omf.mx) + (qn[1] * omf.my) + (qn[2] * omf.mz)

        qcx = omf.mx - asarray(qn[0]*QdotM)
        qcy = omf.my - asarray(qn[1]*QdotM)
        qcz = omf.mz - asarray(qn[2]*QdotM)

        e = exp_part(q[0], q[1], q[2], strucUnit.value_list)


        ruu = (strucUnit.unit) + (qcx*magCell)

        rdd = (strucUnit.unit) - (qcx*magCell)


        rud = (qcy + (1j*qcz)) * magCell

        rdu = (qcy - (1j*qcz)) * magCell

        recip[0] = ruu * e
        recip[1] = rdd * e
        recip[2] = rud * e
        recip[3] = rdu * e


        for i in range(4):
            form_factor[i] = sum(recip[i])

        return form_factor[0],form_factor[1],form_factor[2],form_factor[3]


def cudaBA(cell,Q,lattice,beam,precision = 'float32',refract = False):
    '''
    Overview:
        This calculation uses the Cuda version of the scattering calculator to
    solve the Born Approximation. Although this is slower than the CZT version,
    it is still faster than a direct Python solution.


    Parameters:

    cell(float:3D array|angstroms^2) = The structural scattering potential
    of the feature being scattered off of.

    Q:(q_space) = A Q_space object that holds all of the information about the
    desired q space output.

    lattice:(Lattice) = A lattice object that holds all of the information
    needed to solve the structure factor of the scattering.

    beam:(Beam) = Holds all of the information about the experimental beam
    needed to apply beam dependent corrections to the data.
    '''
    from threadSMBA import cudaBA_form

    form_factor = cudaBA_form(cell,Q,lattice,beam,
                              precision = precision,refract = refract)

    #structure_factor = lattice.struc_calc(Q)
    structure_factor = lattice.gauss_struc_calc(Q)

    intensity = (structure_factor)**2 * abs(form_factor)**2
    #intensity = abs(form_factor)**2

    return normalize(intensity, cell, Q, lattice).real


def cudaSMBA(cell,Q,lattice,beam,precision = 'float32', refract = True):
    '''
    Overview:
        This calculation uses the Cuda version of the scattering calculator to
    solve the Substrate Modified Born Approximation. It perturb the scattering
    by the solution to the wave function for the substrate and incident media.


    Parameters:

    cell(float:3D array|angstroms^2) = The structural scattering potential
    of the feature being scattered off of.

    Q:(q_space) = A Q_space object that holds all of the information about the
    desired q space output.

    lattice:(Lattice) = A lattice object that holds all of the information
    needed to solve the structure factor of the scattering.

    beam:(Beam) = Holds all of the information about the experimental beam
    needed to apply beam dependent corrections to the data.
    '''

    from smba_driver import cudaSMBA_form

    form_factor = cudaSMBA_form(cell,Q,lattice,beam, precision = precision)

    #structure_factor = lattice.struc_calc(Q)
    structure_factor = lattice.gauss_struc_calc(Q)

    intensity = (structure_factor) * abs(form_factor)**2

    return normalize(intensity, cell, Q, lattice).real


def SMBA(cell,Q,lattice,beam):
    '''
    Overview:
        This calculation uses the Python version of the scattering calculator to
    solve the Substrate Modified Born Approximation. It perturb the scattering
    by the solution to the wave function for the substrate and incident media.


    Parameters:

    cell(float:3D array|angstroms^2) = The structural scattering potential
    of the feature being scattered off of.

    Q:(q_space) = A Q_space object that holds all of the information about the
    desired q space output.

    lattice:(Lattice) = A lattice object that holds all of the information
    needed to solve the structure factor of the scattering.

    beam:(Beam) = Holds all of the information about the experimental beam
    needed to apply beam dependent corrections to the data.
    '''

    form_factor = zeros(Q.points,dtype='complex')
    #structure_factor = calc_struc(cell,Q,lattice)

    for i in range (size(Q.q_list[2])):
        print i
        for ii in range (size(Q.q_list[1])):
            for iii in range (size(Q.q_list[0])):
                arg = [Q.q_list[2][i], Q, cell,lattice,beam]
                form_factor[iii,ii,i] = SMBA_form(Q.q_list[0][iii],
                                                  Q.q_list[1][ii],arg)

    #structure_factor = lattice.struc_calc(Q)
    structure_factor = lattice.gauss_struc_calc(Q)
    intensity = (structure_factor) * abs(form_factor)**2

    intensity = (normalize(intensity,cell, Q, lattice)).real

    return intensity


def SMBA_form(qx,qy,arg):
    qz = arg[0]
    Q = arg[1]
    cell = arg[2]
    lattice = arg[3]
    beam = arg[4]
    psi_in_terms = []
    psi_out_terms = []
    ki_z,kf_z = QxQyQz_to_k(qx,0.0,qz,beam.wavelength)
    form_factor = 0.0


    if ki_z < 0.0:
        wf1 = wavefunction_BBM.neutron_wavefunction(-ki_z,flipud(cell.inc_sub))
        r_in = wf1.r
        t_in = wf1.t

        kz_in_sample = -wf1.kz_transmitted
        qx -= beam.wavelength*cell.inc_sub[-1,0]

        psi_in_terms.append(t_in*exp(1j*kz_in_sample*cell.step[2]))

    else:
        wf1 = wavefunction_BBM.neutron_wavefunction(ki_z,cell.inc_sub)
        r_in = wf1.r

        kz_in_sample = ki_z

        psi_in_terms.append(1.0*exp(1j*kz_in_sample*cell.step[2]))
        psi_in_terms.append(r_in*exp(-1j*kz_in_sample*cell.step[2]))


    if kf_z < 0.0:
        wf2 = wavefunction_BBM.neutron_wavefunction(-kf_z , cell.inc_sub)
        r_out = wf2.r

        kz_out_sample = kf_z

        psi_out_terms.append(1.0*exp(-1j*kz_out_sample*cell.step[2]))
        psi_out_terms.append(r_out*exp(1j*kz_out_sample*cell.step[2]))

    else:
        wf2 = wavefunction_BBM.neutron_wavefunction(kf_z ,cell.inc_sub)
        r_out = wf2.r
        t_out = wf2.t

        kz_out_sample = wf2.kz_transmitted
        qx += beam.wavelength*cell.inc_sub[-1,0]

        psi_out_terms.append(t_out*exp(-1j*kz_out_sample*cell.step[2]))

    e = exp_part(qx,qy,qz,cell.value_list)

    for t_i in psi_in_terms:
        for t_o in psi_out_terms:
            misc = complete_formula(1,cell.step,[qx,qy,qz])
            form_factor += sum((cell.unit * e*misc*t_i *t_o ))

    structure = lattice.gauss_struc_calc(Q)
    intensity = (structure)**2 * abs(form_factor)**2
    intensity = (complete_formula(form_factor,cell.step,[qx,qy,qz])[0,0,0])

    return intensity


def SMBA_wavecalc(qx,stepqz,stack,wavelength, ki_z,kf_z):


    if ki_z < 0.0:

        wf1 = wavefunction_BBM.neutron_wavefunction(-ki_z,flipud(stack))
        r_in = wf1.r
        t_in = wf1.t
        kz_in_sample = -wf1.kz_transmitted
        qx -= wavelength * stack[-1,0]

        pio = (t_in*exp(1j*kz_in_sample*stepqz))
        pit = 0.0

    else:

        wf1 = wavefunction_BBM.neutron_wavefunction(ki_z,stack)
        r_in = wf1.r
        t_in = wf1.t
        kz_in_sample = complex(ki_z)
        pio = (1.0*exp(1j*kz_in_sample*stepqz))
        pit = (r_in*exp(-1j*kz_in_sample*stepqz))


    if kf_z < 0.0: #leaving
        wf2 = wavefunction_BBM.neutron_wavefunction(-kf_z , stack)
        r_out = wf2.r

        kz_out_sample = complex(kf_z)

        poo = (1.0*exp(-1j*kz_out_sample*stepqz))
        pot = (r_out*exp(1j*kz_out_sample*stepqz))

    else:
        wf2 = wavefunction_BBM.neutron_wavefunction(kf_z ,stack)
        r_out = wf2.r
        t_out = wf2.t

        kz_out_sample = wf2.kz_transmitted
        qx += wavelength * stack[-1,0]

        poo = (t_out*exp(-1j*kz_out_sample*stepqz))
        pot = 0.0
    qx_refract = qx

    return pio,pit,poo,pot,qx_refract



def longBA(cell,Q,lattice,beam):

    flo64 = zeros(Q.points,dtype='complex128')
    flo32 = zeros(Q.points,dtype='complex64')

    new_q = [None]*3
    new_cell = [None]*3

    for i in range(size(Q.q_list[2])):
        print i
        for ii in range (size(Q.q_list[1])):
            for iii in range (size(Q.q_list[0])):

                flo64[iii,ii,i] = longBAForm(Q.q_list[0][iii],
                         Q.q_list[1][ii], Q.q_list[2][i], cell)

    flo64 = complete_formula(flo64,cell.step,Q.q_list)
    form_factor = flo64

    #structure_factor = lattice.struc_calc(Q)
    structure_factor = lattice.gauss_struc_calc(Q)

    intensity = (structure_factor)**2 * abs(form_factor)**2

    return normalize(intensity, cell, Q, lattice).real


def longBAForm(qx,qy,qz,cell):

    e = exp_part(qx,qy,qz,cell.value_list)
    recip = cell.unit * e
    form_factor = sum(recip)

    return form_factor


def exp_part(qx,qy,qz,d):

        I = complex64(1j)
        if (qx.dtype == float64):
            ex = exp(1j*qx*d[0]/2.0)
            ey = exp(1j*qy*d[1]/2.0)
            ez = exp(1j*qz*d[2]/2.0)

        elif (qx.dtype == float32):
            ex = exp(I*qx*d[0])
            ey = exp(I*qy*d[1])
            ez = exp(I*qz*d[2])

        ex = ex.reshape((size(d[0]),1,1))
        ey = ey.reshape((1,size(d[1]),1))
        ez = ez.reshape((1,1,size(d[2])))

        e = ex * ey * ez

        return e

def BA(cell,Q,lattice,beam):
    '''
    Overview:
        Solves the Born Approximation for the off specular scattering using the
    chirp-z transform which allows for the direct selection of areas and
    spacings in reciprocal space.

    Parameters:

    cell(float:3D array|angstroms^2) = The structural scattering potential
    of the feature being scattered off of.

    Q:(q_space) = A Q_space object that holds all of the information about the
    desired q space output.

    lattice:(Lattice) = A lattice object that holds all of the information
    needed to solve the structure factor of the scattering.

    beam:(Beam) = Holds all of the information about the experimental beam
    needed to apply beam dependent corrections to the data.
    '''

    #form_factor = describes the scattering off of the unit cell
    form_factor = BA_FT(cell.unit,cell.step,Q)
    form_factor = complete_formula(form_factor,cell.step,Q.q_list)

    #structure_factor = describes the scattering off of the lattice
    #structure_factor = lattice.struc_calc(Q)
    structure_factor = lattice.gauss_struc_calc(Q)

    intensity = (structure_factor)**2 * abs(form_factor)**2

    return normalize(intensity, cell, Q, lattice).real


def thetaBA(cell,theta,lattice,beam):
    formfactor = zeros(theta.points)
    q = theta.q_calc(beam.wavelength)

    for i in range(size(theta.theta_list[0])):
        print 'at theta_i number ',i
        for ii in range(size(theta.theta_list[1])):

            formfactor[i,ii] = longBAForm(q[0][i,ii],asarray(0.0),q[1][i,ii],cell)

    delta_x = cell.step[0]
    delta_y = cell.step[1]
    delta_z = cell.step[2]
    qx = q[0]
    qz = q[1]
    x_comp = (( 1. - exp(1j*qx*delta_x))*(-1j / qx))
    y_comp = delta_y
    z_comp = (( 1. - exp(1j*qz*delta_z))*(-1j / qz))

    x_comp[qx == 0.0] = delta_x
    z_comp[qz == 0.0] = delta_z

    formfactor *= x_comp*y_comp*z_comp
    structure_factor = lattice.theta_struc_calc(theta)
    return abs(formfactor)**2 * structure_factor**2



def DWBA(cell,Q,lattice,beam):
    '''
    This function solves the Distorted Wave Born Approximation.

    '''

    #form_factor = DWBA_form(cell,Q,lattice,beam)
    structure_factor = calc_struc(cell,Q,lattice)
    intensity = (structure_factor**2)# * abs(form_factor)**2


    return normalize(intensity, cell, Q, lattice)


def part_mag_func(qx,qy,arg):
    recip = [None,None,None,None]
    qz = arg[0]
    Q = arg[1]
    effective_rho = arg[2]
    struc_cell = arg[3]
    form_factor = [zeros((Q.points[0],Q.points[1]),dtype = 'complex'),
                   zeros((Q.points[0],Q.points[1]),dtype = 'complex'),
                   zeros((Q.points[0],Q.points[1]),dtype = 'complex'),
                   zeros((Q.points[0],Q.points[1]),dtype = 'complex')]

    e = exp_part(qx,qy,qz,struc_cell.value_list)

    recip[0] = effective_rho[0] * e
    recip[1] = effective_rho[1] * e
    recip[2] = effective_rho[2] * e
    recip[3] = effective_rho[3] * e

    for i in range(4):
        form_factor[i] = sum(sum(sum(recip[i],axis=2),axis=1))

    return form_factor[0],form_factor[1],form_factor[2],form_factor[3]


def partial_rho(mag_unit,mag_theta,mag_phi):

    rho_NSF = mag_unit * cos(radians(mag_theta)) * sin(radians(mag_phi))
    rho_SF = mag_unit * sin(radians(mag_theta)) * sin(radians(mag_phi))

    return rho_NSF,rho_SF



def normalize(raw_intensity, cell, Q, lattice):
    '''
    Overview:
        This normalizes the scattering based on the surface area that has been
    probed in the calculation.


    Parameters:

    raw_intensity:(cplx,3D array|angstroms^2) = The unnormalized scattering
    intensity.

    cell(Unit_Cell) = a Unit_Cell Object that describes the unit cell that
    created Intensity

    Q:(q_space) = A Q_space object that holds all of the information about the
    desired q space output.

    lattice:(Lattice) = A lattice object that holds all of the information
    needed to solve the structure factor of the scattering.
    '''

    qz_array = Q.q_list[2].reshape(1,1,Q.points[2])

    raw_intensity *= (4.0 * pi / (qz_array * lattice.repeat[0] *
                                  cell.Dxyz[0] * lattice.repeat[1]*
                                  cell.Dxyz[1]))**2

    raw_intensity = sum(raw_intensity,axis=1)


    return raw_intensity



def BA_FT(total_rho,step,Q):
    '''
    Overview:
        Calculates the BA scattering off of a Unit_Cell object. This can be used
    for both magnetic and non-magnetic BA calculations

    cell(float:3D array|angstroms^2) = The structural scattering potential
    of the feature being scattered off of.

    Q:(q_space) = A Q_space object that holds all of the information about the
    desired q space output.

    lattice:(Lattice) = A lattice object that holds all of the information
    needed to solve the structure factor of the scattering.
    '''
    form_factor = trip_czt(total_rho, step, Q.points, Q.minimums, Q.maximums)
    return form_factor


def QxQyQz_to_k(qx,qy,qz,wavelength):
    '''
    Overview:
        A Python implementation of the Q space to k_in k_out conversion.
    It includes the qy component in the magnitude


    Parameters:


    '''

    qx = asarray(qx)
    qy = asarray(qy)
    qz = asarray(qz)
    wavelength = asarray(wavelength)

    if qy == None: qy = asarray([0.0])

    Qmag = sqrt((qx**2) + (qy**2) + (qz**2))

    if Qmag.size == 1:
        if Qmag==0.0: Qmag = 1.0
    else:
        Qmag[Qmag == 0.0] = 1.0

    k0 = 2.0*pi/wavelength
    twoth = 2.0 * arcsin(Qmag/(2.0*k0))
    tilt = arctan2(qx,qz)

    th_in = (twoth/2.0) + tilt
    th_out = (twoth/2.0) - tilt

    kz_in = k0 * sin(th_in)
    kz_out = -k0 * sin(th_out)
    '''
    from numpy import *
    from pylab import *
    tilttest = ((abs(twoth/2.0-tilt)/(abs(twoth/2.0)+abs(tilt))))
    #print tilttest[qx < 0.0]
    imshow(flipud(transpose(log10(sum(tilttest, axis = 1)/tilttest.shape[1]))))
    title('Precision Loss for -')
    colorbar()
    show()
    '''
    return kz_in, kz_out



def complete_formula(czt_result, step, q_list):
    '''
    the scattering from the unit cell is not just the fft of the unit cell.
    There is another mathematical piece that must be solved for. This
    calculation solves for this additional math and multiples the given fft
    by it.

    czt_results = the resulting Fourier transform of the unit cell for the Q
    values requested by the user
    step = the real space step size represented by the discretized units of the
    unit cell
    '''

    qx = q_list[0].reshape(size(q_list[0]),1,1)
    qy = q_list[1].reshape(1,size(q_list[1]),1)
    qz = q_list[2].reshape(1,1,size(q_list[2]))




    I = complex64(1j)
    if (qx[0].dtype == float64):
        delta_x = step[0]
        delta_y = step[1]
        delta_z = step[2]
        x_comp = (( 1. - exp(1j*qx*delta_x))*(-1j / qx))
        y_comp = (( 1. - exp(1j*qy*delta_y))*(-1j / qy))
        z_comp = (( 1. - exp(1j*qz*delta_z))*(-1j / qz))

    elif (qx[0].dtype == float32):
        print 'float32 calculation'
        delta_x = float32(step[0])
        delta_y = float32(step[1])
        delta_z = float32(step[2])
        x_comp = (( 1. - exp(I*qx*delta_x))*(-I / qx))
        y_comp = (( 1. - exp(I*qy*delta_y))*(-I / qy))
        z_comp = (( 1. - exp(I*qz*delta_z))*(-I / qz))


    x_comp[qx == 0.0] = delta_x
    y_comp[qy == 0.0] = delta_y
    z_comp[qz == 0.0] = delta_z

    #Creates individual values of the delta_realspace so that the formalism is
    #easier to read.

    czt_result *= x_comp * y_comp * z_comp

    return czt_result


def trip_czt(unit,step,q_points,q_mins,q_maxs):
    '''
    Given a 3d unit with a set of three limits, this solves the 3d fft using the
    Chirp-z transform.

    unit = The matrix representation of the unit cell
    q_points = the number of points that q will be solved for in the x,y and z
    direction
    q_mins = the minimum q values for the x, y, and z direction
    q_max = the maximum q values for the x, y, and z direction
    '''

    frequancy = (2*pi)/step

    intermediate_matrix_one = czt.zoomfft(unit,q_mins[0],q_maxs[0],q_points[0],
                                          frequancy[0],axis = 0)

    intermediate_matrix_two = czt.zoomfft(intermediate_matrix_one,q_mins[1],
                                          q_maxs[1],q_points[1],frequancy[1],
                                          axis=1)

    czt_result = czt.zoomfft(intermediate_matrix_two,q_mins[2],
                             q_maxs[2],q_points[2],frequancy[2],axis=2)

    return czt_result


def xy_czt(unit,step,q_points,q_mins,q_max):
    '''
    Given a 3d unit with a set of three limits, this solves the 3d fft using the
    Chirp-z transform.

    unit = The matrix representation of the unit cell
    q_points = the number of points that q will be solved for in the x,y and z
    direction
    q_mins = the minimum q values for the x, y, and z direction
    q_max = the maximum q values for the x, y, and z direction
    '''
    max_qxyz = (2*pi)/step
    intermediate = czt.zoomfft(unit,q_mins[0],q_max[0],q_points[0],
                               max_qxyz[0],axis = 0)
    xy_czt_result = czt.zoomfft(intermediate,q_mins[1],q_max[1],
                                q_points[1],max_qxyz[1],axis=1)

    return xy_czt_result


def calc_struc(cell,Q,lattice):
    '''
    This module calculates the structure factor of a repeating structure
    in  the case of a parallelapiped lattice

    cell = the Unit_Cell object that holds all of the information about the
    unit cell
    Q = the Q space requested to be solved for by the user
    lattice = the lattice parameters of the repeat structure

    '''
    if size(Q)== 3:
        qx = asarray([Q[0]])
        qy = asarray([Q[1]])
        qz = asarray([Q[2]])

    else:
        qx,qy,qz = Q.vectorize()

    struc_x =((sin((qx*cell.Dxyz[0]*lattice.repeat[0])/2) /
               sin((qx*cell.Dxyz[0])/2)))
    struc_y =((sin((qy*cell.Dxyz[1]*lattice.repeat[1])/2) /
               sin((qy*cell.Dxyz[1])/2)))
    struc_z =((sin((qz*cell.Dxyz[2]*lattice.repeat[2])/2) /
               sin((qz*cell.Dxyz[2])/2)))

    struc_x[qx == 0.0] = lattice.repeat[0]
    struc_y[qy == 0.0] = lattice.repeat[1]
    struc_z[qz == 0.0] = lattice.repeat[2]

    structure_factor = struc_x * struc_y * struc_z
    return structure_factor




def wavefunction_format(unit, z_step, absorbtion = None):
    '''
    from a given 3 dimensional matrix, this takes the average of the xy plane
    and creates a variable to be returned in the form:
    SLD_WF_format [:,0] = the average SLD from the top of the sample to the
    bottom
    SLD_WF_format [:,1] = The depth of the layer from  the top of the sample
    to the bottom
    SLD_WF_format [:,2] = the average absorption from the top of the sample
    to the bottom

    unit = a three dimensional array that represents the unit cell
    z_values = an array of values that represent the real space locations of
    each layer
    '''

    SLD_list = (apply_over_axes(average, unit, [0,1])).reshape(size(
                                                                unit[1,1,:]))

    layer_count = size(SLD_list)
    SLD_WF_format = empty([layer_count,3])

    SLD_WF_format [:,0] = SLD_list[::-1]
    SLD_WF_format [:,1] = z_step

    if absorbtion == None:
        SLD_WF_format [:,2] = 0.0
    else:
        SLD_WF_format [:,2] = absorbtion[::-1]

    return SLD_WF_format


def _test():
    '''
    import sample_prep
    import scatter
    import pylab


    #thetaBA test
    Au = sample_prep.Ellipse(SLD = 4.506842e-6,dim=[3.75e4,3.75e4,600.0])
    Cr = sample_prep.Layer(SLD = 3.01e-6,thickness_value = 20.0)

    #Au.on_top_of(Cr)
    scene = sample_prep.Scene([Au])

    GeoUnit = sample_prep.GeomUnit(Dxyz = [10.0e4,10.0e4, 2000.0], n = [50,51,52],scene = scene)
    unit = GeoUnit.buildUnit()


    q_space = sample_prep.Q_space([-.0006,-0.001,0.00002],[.0002,0.001,0.025],[150,10,150])

    theta = sample_prep.Theta_space([-1.0,0.0],[1.0,1.2],[300,300])
    lattice = sample_prep.Rectilinear([20,20,1],unit)

    beam = sample_prep.Beam(5.0,None,None,0.05,None)
    sample = scatter.Calculator(lattice,beam,theta,unit)
    sample.BA()


    pylab.figure()


    extent = [theta.minimums[1],theta.maximums[1],theta.minimums[0],theta.maximums[0]]
    pylab.imshow((flipud(log(sample.results))),extent=extent,aspect = 'auto')
    colorbar()
    q = theta.q_calc(5.0)
    pylab.figure()
    pylab.pcolormesh(q[0],q[1],(log(sample.results)))
    colorbar()
    show()
    '''
    from numpy import shape
    from sample_prep import *
    from scatter import *
    from pylab import *
    from omfLoader import *
    #magneticBA test
    mag = Omf('/home/mettingc/Documents/test.omf')
    #mag.view()



    Au = (sample_prep.Parallelapiped(SLD = 4.506842e-6,
                                     Ms = 8.6e5,dim=[5.0e4,5.0e4,2.0e4]))

    Cr = (sample_prep.Layer(SLD = 3.01e-6,Ms = 8.6e5,
                            thickness_value = 1000.0))

    Au.on_top_of(Cr)
    scene = sample_prep.Scene([Au,Cr])

    GeoUnit = (sample_prep.GeomUnit(Dxyz = [10.0e4,10.0e4,2.5e4],
                                    n = [70,70,30],scene = scene))

    unit = GeoUnit.buildUnit()
    unit.generateMIF()
    #unit.viewSlice()
    space = Q_space([-.0001,-0.001,0.00002],[.0001,0.001,0.04],[50,5,50])
    lattice = Rectilinear([20,20,1],unit)

    beam = Beam(5.0,None,None,0.05,None)


    magtest = Calculator(lattice, beam, space, unit, mag)

    magtest.magneticBA()
    magtest.view_uncorrected()
    show()
if __name__=="__main__":_test()
