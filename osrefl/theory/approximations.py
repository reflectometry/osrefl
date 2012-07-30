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
from osrefl.theory import wavefunction_kernel
import osrefl.viewers.view
from osrefl.theory import czt
import numpy as np

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

        form_factor[i] = complete_formula(form_factor[i],struc_cell.step,Q)

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
        form_factor[i] = complete_formula(form_factor[i],struc_cell.step,Q)


        intensity[i] = abs(structure_factor)**2 * abs(form_factor[i])**2
        intensity[i] = (normalize(intensity[i], struc_cell, Q, lattice)).real

    return intensity

def cudaMagBA(struc_cell, Q, lattice, beam, magVec, rho_m, precision = 'float32',
              refract = False):

    from smba_driver import magCudaBA_form
    intensity = [None]*4

    form_factor = magCudaBA_form(struc_cell,Q,lattice,beam,magVec,rho_m,
                              precision = precision,refract = refract)

    #structure_factor = lattice.struc_calc(Q)
    structure_factor = lattice.gauss_struc_calc(Q)

    for i in range(4):
        #intensity[i] = sum(abs(form_factor[i])**2,axis = 1)
        intensity[i] = (structure_factor)**2 * abs(form_factor[i])**2
        intensity[i] = normalize(intensity[i],struc_cell, Q, lattice).real

    return intensity

def magneticBA(struc_cell,Q,lattice,beam,magVec,rho_m):
    '''
    **Overview:**

        This calculation solves the Born Approximation for a magnetic sample
        using a(n) :class:`~omfLoader.Omf` object which hold the information about
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
    #rho_m = omf.rhoM
    #magVec = [omf.mx, omf.my, omf.mz]

    for i in range (size(Q.q_list[0])):
        print i
        for ii in range (size(Q.q_list[1])):
            for iii in range (size(Q.q_list[2])):

                q = [Q.q_list[0][i],Q.q_list[1][ii],Q.q_list[2][iii]]
                qn = [qxn[i,0,0],qyn[0,ii,0],qzn[0,0,iii]]

                (form_factor[0][i,ii,iii],form_factor[1][i,ii,iii],
                 form_factor[2][i,ii,iii], form_factor[3][i,ii,iii]) = mag_func(
                                        struc_cell, magVec, q, qn, rho_m)

    structure_factor = lattice.gauss_struc_calc(Q)

    for i in range(4):
        form_factor[i] = complete_formula(form_factor[i],struc_cell.step,Q)
        #intensity[i] = sum(abs(form_factor[i])**2,axis = 1)
        intensity[i] = (structure_factor)**2 * abs(form_factor[i])**2
        intensity[i] = normalize(intensity[i],struc_cell, Q, lattice).real

    return intensity, form_factor, structure_factor


def mag_func(strucUnit, magVec, q, qn, rho_m):
    mx, my, mz = magVec
    recip = [None,None,None,None]
    #form_factor = [None,None,None,None]

    QdotM = (qn[0] * mx) + (qn[1] * my) + (qn[2] * mz)

    qcx = mx - asarray(qn[0]*QdotM)
    qcy = my - asarray(qn[1]*QdotM)
    qcz = mz - asarray(qn[2]*QdotM)
    
    qca = qcx * 0.0 + qcy * 1.0 + qcz * 0.0;
    qcb = qcx * 0.0 + qcy * 0.0 + qcz * 1.0;
    qcc = qcx * 1.0 + qcy * 0.0 + qcz * 0.0;
    
    x,y,z = strucUnit.value_list
    
    e = exp(1j * (q[0] *  x[:, newaxis, newaxis] + q[1] * y[newaxis, :, newaxis] + q[2] * z[newaxis, newaxis, :]))
    #e = exp_part(q[0], q[1], q[2], strucUnit.value_list)
    
    ruu = (strucUnit.unit) + (qca*rho_m)

    rdd = (strucUnit.unit) - (qca*rho_m)


    rud = (qcb + (1j*qcc)) * rho_m

    rdu = (qcb - (1j*qcc)) * rho_m

    recip[0] = ruu * e
    recip[1] = rdd * e
    recip[2] = rud * e
    recip[3] = rdu * e


    for i in range(4):
        form_factor[i] = sum(recip[i])

    return form_factor[0],form_factor[1],form_factor[2],form_factor[3]


def cudaBA(cell,Q,lattice,beam,precision = 'float32',refract = True):
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

    from smba_driver import cudaBA_form
    form_factor = cudaBA_form(cell,Q,lattice,beam,
                              precision = precision,refract = refract)
    if refract == True:
        strucRefract = True
    else:
        strucRefract = False

    structure_factor = lattice.gauss_struc_calc(Q,strucRefract)
    #structure_factor = lattice.struc_calc(Q)

    intensity = (structure_factor)**2 * abs(form_factor)**2

    return normalize(intensity, cell, Q, lattice).real

def SMBAfft(cell,Q,lattice,beam,precision = 'float32',refract = True):
    '''
    Overview:
        It may be true that the SMBA can be solved by using a czt transform and
    perturbing this by the wavefunctions. This calculation will solve the czt
    and then perturb it by the wave function.

    
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

    from scipy.interpolate import RectBivariateSpline, interp1d
    from ..model.sample_prep import Q_space
    from smba_wave_driver import wave
    #from pylab import *

    stack = cell.inc_sub
    #stack = wavefunction_format(cell.unit, cell.step[2], absorbtion = None)

    #This is stupid. I don't have a function that produces a refraction
    #so I reuse the psi calculation twice. This needs to be fixed
    psi_in_one,psi_in_two,psi_out_one,psi_out_two,qx_refract = (
                wave(stack, Q.q_list[0], Q.q_list[1],
                     Q.q_list[2],beam.wavelength,cell.step[2],
                     precision=precision))
    Q.qx_refract = qx_refract

    if refract == True:
        scat = zeros(Q.points, dtype = 'complex')

        qxMinTemp = Q.minimums[0]-(stack[0,-1]*beam.wavelength)-2*Q.q_step[0]
        qxMaxTemp = Q.maximums[0]+(stack[0,-1]*beam.wavelength)+2*Q.q_step[0]

        newX = arange(qxMinTemp,qxMaxTemp,Q.q_step[0])

        newQ = Q_space([qxMinTemp,Q.minimums[1],Q.minimums[2]],
                       [qxMaxTemp,Q.maximums[1],Q.maximums[2]],
                       [size(newX),Q.points[1],Q.points[2]])

        form_factor = zeros(newQ.points, dtype = 'complex')

        psi_in_one,psi_in_two,psi_out_one,psi_out_two,qx_refract = (
                        wave(stack, newQ.q_list[0], newQ.q_list[1],
                             newQ.q_list[2],beam.wavelength,cell.step[2],
                             precision=precision))

        temp_form_factor = BA_FT(cell.unit,cell.step,newQ)

        scatProc = [None]*4
        scatProc[0] = psi_in_one*temp_form_factor*psi_out_one
        scatProc[1] = psi_in_one*temp_form_factor*psi_out_two
        scatProc[2] = psi_in_two*temp_form_factor*psi_out_one
        scatProc[3] = psi_in_two*temp_form_factor*psi_out_two

        temp_form_factor = scatProc[0]+scatProc[1]+scatProc[2]+scatProc[3]

        temp_form_factor = complete_formula(temp_form_factor,cell.step,newQ)

        structure_factor = lattice.gauss_struc_calc(newQ)

        intensity = abs(structure_factor)**2 * abs(temp_form_factor)**2

        for ii in range (size(Q.q_list[1])):
            for iii in range(size(Q.q_list[2])):
                realSplineFunc = interp1d(newQ.q_list[0],intensity.real[:,ii,iii])
                imagSplineFunc = interp1d(newQ.q_list[0],intensity.imag[:,ii,iii])

                interpReal = realSplineFunc(Q.qx_refract[:,ii,iii])
                interpImag = imagSplineFunc(Q.qx_refract[:,ii,iii])

                scat[:,ii,iii].real = interpReal
                scat[:,ii,iii].imag = interpImag
    else:

        form_factor = BA_FT(cell.unit,cell.step,Q)

        scatProc = [None]*4
        scatProc[0] = psi_in_one*form_factor*psi_out_one
        scatProc[1] = psi_in_one*form_factor*psi_out_two
        scatProc[2] = psi_in_two*form_factor*psi_out_one
        scatProc[3] = psi_in_two*form_factor*psi_out_two

        form_factor = scatProc[0]+scatProc[1]+scatProc[2]+scatProc[3]

        form_factor = complete_formula(form_factor,cell.step,Q)

        structure_factor = lattice.gauss_struc_calc(Q)

        scat = abs(structure_factor)**2 * abs(form_factor)**2

    return normalize(scat, cell, Q, lattice).real



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

    form_factor = cudaSMBA_form(cell,Q,lattice,beam, precision = precision,
                                refract = refract)

    #structure_factor = lattice.struc_calc(Q)
    if refract == True:
        strucRefract = True
    else:
        strucRefract = False

    structure_factor = lattice.gauss_struc_calc(Q,strucRefract)

    intensity = abs(structure_factor)**2 * abs(form_factor)**2

    return normalize(intensity, cell, Q, lattice)


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
    qx_refract = zeros(Q.points,dtype='complex')
    #structure_factor = calc_struc(cell,Q,lattice)

    for i in range (size(Q.q_list[2])):
        print i
        for ii in range (size(Q.q_list[1])):
            for iii in range (size(Q.q_list[0])):
                arg = [Q.q_list[2][i], Q, cell,lattice,beam]
                form_factor[iii,ii,i],qx_refract[iii,ii,i] = SMBA_form(
                                                           Q.q_list[0][iii],
                                                           Q.q_list[1][ii],arg)
    Q.qx_refract = qx_refract
    complete_formula(form_factor,cell.step,Q)
    #structure_factor = lattice.struc_calc(Q)
    structure_factor = lattice.gauss_struc_calc(Q)
    intensity = (structure_factor)**2 * abs(form_factor)**2

    intensity = normalize(intensity,cell, Q, lattice)

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
            misc = complete_formula(1,cell.step,Q)
            form_factor += sum(cell.unit * e*misc*t_i *t_o )

    #intensity = (complete_formula(form_factor,cell.step,[qx,qy,qz])[0,0,0])

    return form_factor,qx

def interpWaveCalc(cell,Q,beam,proc='gpu'):

    '''
    UNDER CONSTRUCTION
    '''
    from scipy.interpolate import RectBivariateSpline
    import wavefunction_kernel
    from smba_wave_driver import wave
    import time

    Q.getKSpace(beam.wavelength)

    klim = array([amin(Q.kin),amin(Q.kout),amax(Q.kin),amax(Q.kout)])
    kn = array([5*max(Q.points),5*max(Q.points)])
    #kn = array([2000,2000])

    kinVec = linspace(klim[0],klim[2],kn[0])
    koutVec = linspace(klim[1],klim[3],kn[1])

    kstep = array([kinVec[1]-kinVec[0],koutVec[1]-koutVec[0]])

    pioTab=ones(kn,dtype = 'complex')
    pitTab=ones(kn,dtype = 'complex')
    pooTab=ones(kn,dtype = 'complex')
    potTab=ones(kn,dtype = 'complex')
    qxRefractTrax = ones(kn)

    pio = zeros(Q.points, dtype = 'complex')
    pit = zeros(Q.points, dtype = 'complex')
    poo = zeros(Q.points, dtype = 'complex')
    pot = zeros(Q.points, dtype = 'complex')

    from pylab import imshow,show,colorbar,figure,subplot
    '''
    #Create lookup table
    for i in range(size(kinVec)):
        print i
        for ii in range(size(koutVec)):

            pioTab[i,ii],pitTab[i,ii],pooTab[i,ii],potTab[i,ii],qxRefractTrax[i,ii] = (
                                   SMBA_wavecalc(qxRefractTrax[i,ii],
                                                 cell.step[2],
                                                 cell.inc_sub, beam.wavelength,
                                                 kinVec[i], koutVec[ii]))

    realPsiInterp = [None]*4
    imagPsiInterp = [None]*4
    psiTabGroup = [pioTab,pitTab,pooTab,potTab]
    psiGroup = [pio,pit,poo,pot]

    print id(pioTab)
    print id(psiGroup[i])
    show()

    for i in range(4):

        realPsiInterp[i] = RectBivariateSpline(kinVec,koutVec,
                                               psiTabGroup[i].real)

        imagPsiInterp[i] = RectBivariateSpline(kinVec,koutVec,
                                               psiTabGroup[i].imag)


        for v in range(Q.points[1]):
            for vv in range (Q.points[2]):
                psiGroup[i][:,v,vv].real =realPsiInterp[i].ev(Q.kin[:,v,vv],Q.kout[:,v,vv])
                psiGroup[i][:,v,vv].imag =imagPsiInterp[i].ev(Q.kin[:,v,vv],Q.kout[:,v,vv])
    '''


    psi_in_one,psi_in_two,psi_out_one,psi_out_two,qx_refract = (
                    wave(cell.inc_sub, Q.q_list[0], Q.q_list[1],
                         Q.q_list[2],beam.wavelength,cell.step[2], proc=proc))

    pioTab=ones(Q.points,dtype = 'complex')
    pitTab=ones(Q.points,dtype = 'complex')
    pooTab=ones(Q.points,dtype = 'complex')
    potTab=ones(Q.points,dtype = 'complex')

    '''
    for i in range(Q.points[0]):
        print i
        for ii in range(Q.points[1]):
            for iii in range(Q.points[2]):
                pioTab[i,ii,iii],pitTab[i,ii,iii],pooTab[i,ii,iii],potTab[i,ii,iii],temp = (
                                       SMBA_wavecalc(qxRefractTrax[i,ii],
                                                     cell.step[2],
                                                     cell.inc_sub, beam.wavelength,
                                                     Q.kin[i,ii,iii], Q.kout[i,ii,iii]))
    '''
    qxRefractTrax = zeros(Q.points)
    t = time.time()
    pioTab,pitTab,pooTab,potTab,qxRefractTrax = (SMBA_wavecalcMultiK(qxRefractTrax,
                                         cell.step[2],
                                         cell.inc_sub, beam.wavelength,
                                         Q.kin, Q.kout))
    print time.time()-t
    '''
    pio = psiGroup[0]
    pit = psiGroup[1]
    poo = psiGroup[2]
    pot = psiGroup[3]
    '''
    '''
    imshow(sum(psi_in_one,axis = 1).real)
    colorbar()
    figure()
    imshow(sum(pioTab,axis=1).real)
    colorbar()
    figure()
    imshow(sum(abs(psi_in_one-pioTab)/sqrt(abs(psi_in_one)+abs(pioTab)),axis = 1).real)
    colorbar()
    show()

    subplot(2,2,1)
    imshow(rot90(log10(sum(psi_in_one,axis=1).real)))
    subplot(2,2,2)
    imshow(rot90(log10(sum(psi_in_two,axis=1).real)))
    subplot(2,2,3)
    imshow(rot90(log10(sum(psi_out_one,axis=1).real)))
    subplot(2,2,4)
    imshow(rot90(log10(sum(psi_out_two,axis=1).real)))
    colorbar()

    figure()
    subplot(2,2,1)
    imshow(rot90(log10(sum(pioTab,axis=1).real)))
    subplot(2,2,2)
    imshow(rot90(log10(sum(pitTab,axis=1).real)))
    subplot(2,2,3)
    imshow(rot90(log10(sum(pooTab,axis=1).real)))
    subplot(2,2,4)
    imshow(rot90(log10(sum(potTab,axis=1).real)))
    colorbar()

    figure()
    subplot(2,2,1)
    imshow(rot90(sum(abs(psi_in_one-pioTab)/sqrt(abs(psi_in_one)+abs(pioTab)),axis = 1).real))
    subplot(2,2,2)
    imshow(rot90(sum(abs(psi_in_two-pitTab)/sqrt(abs(psi_in_one)+abs(pitTab)),axis = 1).real))
    subplot(2,2,3)
    imshow(rot90(sum(abs(psi_out_one-pooTab)/sqrt(abs(psi_in_one)+abs(pooTab)),axis = 1).real))
    subplot(2,2,4)
    imshow(rot90(sum(abs(psi_out_two-potTab)/sqrt(abs(psi_in_one)+abs(potTab)),axis = 1).real))
    colorbar()

    show()
    '''
    qxrefract = ((Q.vectorize()[0]).repeat(Q.points[1],axis=1)).repeat(
                                                           Q.points[2],axis=2)

    qxrefract[Q.kin < 0.0] += beam.wavelength * cell.inc_sub[-1,0]
    qxrefract[Q.kout > 0.0] -= beam.wavelength * cell.inc_sub[-1,0]

    #return psi_in_one,psi_in_two,psi_out_one,psi_out_two,qxrefract
    return pioTab,pitTab,pooTab,potTab,qxrefract


def SMBA_wavecalcMultiK(qx,deltaz,stack,wavelength, ki_z,kf_z):

    pio = zeros(shape(ki_z),dtype = 'complex')
    pit = zeros(shape(ki_z),dtype = 'complex')
    poo = zeros(shape(ki_z),dtype = 'complex')
    pot = zeros(shape(ki_z),dtype = 'complex')

    wfinn = smbaWavefunction(-ki_z, flipud(stack))
    wfinp = smbaWavefunction(ki_z, stack)
    wfoutn = smbaWavefunction(kf_z , stack)
    wfoutp = smbaWavefunction(kf_z , stack)

    pio[ki_z<0.0] =  wfinn.t[ki_z<0.0]*exp(1j*(
                                   -wfinn.kz_transmitted[ki_z<0.0])*deltaz)
    pit[ki_z<0.0] = 0.0

    pio[ki_z>=0.0] = 1.0*exp(1j*ki_z[ki_z>0.0]*deltaz)
    pit[ki_z>=0.0] = wfinp.r[ki_z>=0.0]*exp(-1j*ki_z[ki_z>0.0]*deltaz)

    poo[kf_z<0.0] = 1.0*exp(-1j*kf_z[kf_z<0.0]*deltaz)
    pot[kf_z<0.0] =  wfoutn.r[kf_z<0.0]*exp(1j*kf_z[kf_z<0.0]*deltaz)

    poo[kf_z>=0.0] = (wfoutp.t[kf_z>=0.0]*exp(
                                  -1j*wfoutp.kz_transmitted[kf_z>=0.0]*deltaz))
    pot[kf_z>=0.0] = 0.0

    qx[ki_z < 0.0] += wavelength * stack[-1,0]
    qx[kf_z > 0.0] -= wavelength * stack[-1,0]


    return pio, pit, poo, pot,qx


def SMBA_wavecalc(qx,deltaz,stack,wavelength, ki_z,kf_z):


    if ki_z < 0.0:

        wf1 = wavefunction_BBM.neutron_wavefunction(-ki_z,flipud(stack))
        r_in = wf1.r
        t_in = wf1.t
        kz_in_sample = -wf1.kz_transmitted
        qx -= wavelength * stack[-1,0]

        pio = (t_in*exp(1j*kz_in_sample*deltaz))
        pit = 0.0

    else:

        wf1 = wavefunction_BBM.neutron_wavefunction(ki_z,stack)
        r_in = wf1.r
        t_in = wf1.t
        kz_in_sample = complex(ki_z)

        pio = (1.0*exp(1j*kz_in_sample*deltaz))
        pit = (r_in*exp(-1j*kz_in_sample*deltaz))


    if kf_z < 0.0: #leaving
        wf2 = wavefunction_BBM.neutron_wavefunction(kf_z , stack)
        r_out = wf2.r

        kz_out_sample = complex(kf_z)

        poo = (1.0*exp(-1j*kz_out_sample*deltaz))
        pot = (r_out*exp(1j*kz_out_sample*deltaz))

    else:
        wf2 = wavefunction_BBM.neutron_wavefunction(kf_z ,stack)
        r_out = wf2.r
        t_out = wf2.t

        kz_out_sample = wf2.kz_transmitted
        qx += wavelength * stack[-1,0]

        poo = (t_out*exp(-1j*kz_out_sample*deltaz))
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

    flo64 = complete_formula(flo64,cell.step,Q)
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

def BAres(cell,q,lattice,beam):
    from DWBA import dwbaWavefunction
    from scipy.integrate import quad
    m = 1.674e-27
    h_bar = 6.62607e-14

    Vfac = -m/(2*pi*h_bar**2)
    ftwRef = zeros([q.points[0],q.points[0],cell.n[2]],dtype = 'complex')
    scat = zeros(q.points,dtype = 'complex')
    x = cell.value_list[0].reshape((cell.n[0],1,1))
    y = cell.value_list[1].reshape((1,cell.n[1],1))
    z = cell.value_list[2].reshape((1,1,cell.n[2]))
    
    flipCell = zeros(shape(cell.unit))
    SLDArray = wavefunction_format(cell.unit, cell.step[2], absorbtion = None)
    
    
    
    for i in range(cell.n[2]):
        flipCell[:,:,i] = cell.unit[:,:,shape(cell.unit)[2]-i-1]
    
    flipCell2 = flipud(cell.unit)

    Vres = flipCell - (SLDArray[:,0]).reshape((1,1,cell.n[2]))

    rhoTilOverRho = Vres/(SLDArray[:,0]).reshape((1,1,cell.n[2]))
    rhoTilOverRho[isnan(rhoTilOverRho)] = 0.0

    SF = lattice.gauss_struc_calc(q)

    for i in range(size(q.q_list[0])):
        print 'qx number: ', i, ' calculating'

        for ii in range(size(q.q_list[1])):
            
            #Equation 20-
            laux = ((-1j / q.q_list[0][i]) * 
                    (exp(1j *q.q_list[0][i] * cell.step[0]) - 1.0))
            lauy = ((-1j / q.q_list[1][ii]) * 
                    (exp(1j *q.q_list[1][ii] * cell.step[1]) - 1.0))

            if isnan(laux):
                laux = cell.step[0]
            if isnan(lauy):
                lauy = cell.step[1]
                
            #ftwRef = (Vfac)*sum(sum(Vres * exp(1j*q.q_list[0][i]*x)*exp(1j*q.q_list[1][ii]*y),axis = 0),axis=0)
            ftHoldx = zeros([cell.n[0]], dtype = 'complex')
            ftHoldy = zeros([cell.n[1]], dtype = 'complex')
            ftwRef = zeros(cell.n,dtype = 'complex')
            ftwCollector = zeros([size(q.q_list[0]),size(q.q_list[1]),cell.n[2]],dtype = 'complex')
            for nz in range(cell.n[2]):
                for nx in range(cell.n[0]):
                    for ny in range(cell.n[1]):
                        args = [x[nx,0,0],0,rhoTilOverRho[nx,ny,nz]]
                        a = q.q_list[0][i] - (q.q_step[0]/2.)
                        b = q.q_list[0][i] + (q.q_step[0]/2.)
                        
                        ftHoldx = quad(ftIntegral,a,b,args)[0]/q.q_step[0]
                        
                        args = [y[0,ny,0],0,rhoTilOverRho[nx,ny,nz]]
                        a = q.q_list[1][ii] - (q.q_step[1]/2.)
                        b = q.q_list[1][ii] + (q.q_step[1]/2.)
                        ftHoldy = quad(ftIntegral,a,b,args)[0]/q.q_step[1]

                        ftwRef[nx,ny,nz] = ftHoldx * ftHoldy
            
            ftwCollector[i,ii,:] = (sum(sum(ftwRef,axis = 0),axis = 0)).reshape((1,1,cell.n[2]))

            '''   
            ftwRef[i,ii,:] = (Vfac*sum(sum(rhoTilOverRho * exp(1j*q.q_list[0][i]*x)*
                                   exp(1j*q.q_list[1][ii]*y),axis = 0),axis=0))
            '''
            ftwCollector[i,ii,:] *= laux
            ftwCollector[i,ii,:] *= lauy

            
            ftwCollector[i,ii,:] *=SF[i,ii,0]

            ftwCollector[i,ii,:] = ftwCollector[i,ii,:]/(lattice.repeat[0]*lattice.repeat[1])
            
            ftwCollector[i,ii,:] = ((SLDArray[:,0]).reshape((1,1,cell.n[2]))*
                      ftwCollector[i,ii,:].reshape((1,1,cell.n[2])))

    max_qxyz = (2*pi)/cell.step
    scat = czt.zoomfft(ftwRef,q.minimums[2],q.maximums[2],q.points[2],
                               max_qxyz[2],axis = 2)
    
    vecq = q.vectorize()
    lauz = ((-1j / vecq[2]) * (exp(1j *vecq[2] * cell.step[2]) - 1.0))
    scat *= lauz
    intensity =sum(abs(scat)**2,axis=1)
    
    k_spec = q.q_list[2]/2.0
    dwba_spec = dwbaWavefunction(k_spec,SLDArray)

    locx = q.q_list[0].searchsorted(0.0)
    locy = q.q_list[1].searchsorted(0.0)

    scat[locx,locy,:] = dwba_spec.r
    return intensity

def ftIntegral(q,args):
    d = args[0]
    axis = args[1]
    SLD = args[2]
    return sum(exp(1j*q*d).real,axis = axis)


def BA(cell,Q,lattice,beam):
    '''
    Overview:
        Solves the Born Approximation for the off-specular scattering using the
    chirp-z transform which allows for the direct selection of areas and
    spacings in reciprocal space.+-

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
    #form_factor = complete_formula(form_factor,cell.step,Q)

    #structure_factor = describes the scattering off of the lattice
    structure_factor = lattice.gauss_struc_calc(Q)

    intensity = abs(form_factor)**2 * abs(structure_factor)**2
    return normalize(intensity, cell, Q, lattice)


def BA_FF(cell,Q):
    '''
    Overview:
        Solves the Born Approximation for the off-specular scattering using the
    chirp-z transform. Calculates only the form factor.

    Parameters:

    cell(float:3D array|angstroms^2) = The structural scattering potential
    of the feature being scattered off of.

    Q:(q_space) = A Q_space object that holds all of the information about the
    desired q space output.

    beam:(Beam) = Holds all of the information about the experimental beam
    needed to apply beam dependent corrections to the data.
    '''
    
    raw_intensity = approximations.BA_FT(cell.unit, cell.step, Q)
    raw_intensity = abs(approximations.complete_formula(raw_intensity, cell.step, Q))**2

    qx_array = Q.q_list[0].reshape(Q.points[0],1,1)
    qy_array = Q.q_list[1].reshape(1,Q.points[1], 1)
    qz_array = Q.q_list[2].reshape(1,1,Q.points[2])

    raw_intensity *= (4.0 * pi / (qz_array * qx_array * qy_array))**2

    #raw_intensity = sum(raw_intensity,axis=1).astype('float64')  

    return raw_intensity


def thetaBA(cell,theta,lattice,beam):
    formfactor = zeros(theta.points)
    q = theta.q_calc(beam.wavelength)

    for i in range(size(theta.theta_list[0])):
        print 'at theta_i number ',i
        for ii in range(size(theta.theta_list[1])):

            formfactor[i,ii] = longBAForm(q[0][i,ii],asarray(0.0),q[1][i,ii],
                                          cell)

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
    return abs(formfactor)**2 * abs(structure_factor)**2

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

    '''

    form_factor = trip_czt(total_rho, step, Q.points, Q.minimums, Q.maximums)

    return form_factor


def old_QxQyQz_to_k(qx,qy,qz,wavelength):
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

    #Qmag = sqrt((qx**2) + (qy**2) + (qz**2))
    Qmag = sqrt((qx**2) + zeros_like(qy)**2 +  (qz**2))

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

    return kz_in, kz_out

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
    #Qmag = sqrt((qx**2) + zeros_like(qy)**2 +  (qz**2))

    if Qmag.size == 1:
        if Qmag==0.0: Qmag = 1.0
    else:
        Qmag[Qmag == 0.0] = 1.0

    k0 = 2.0*pi/wavelength
    kfy = abs(qy) # just looking for length here
    kfxz = sqrt(k0**2 - kfy**2)
    Qxzmag = sqrt(qx**2 + qz**2)
    tilt = arctan2(qx,qz)
    
    cos_gamma = Qmag**2 / ( 2 * Qxzmag * k0 )
    sin_gamma = sqrt(1 - cos_gamma**2)
    
    #th_in = pi/2.0 - tilt - arccos(Qmag**2 / (2*Qxzmag*k0))
    #th_out = -th_in + arccos((2*k0**2 - Qmag**2)/(2*kfxz*k0))

    #kz_in = k0 * sin(th_in)
    #kz_out = -k0 * sin(th_out)
    
    kz_in = k0 * (cos(tilt) * cos_gamma - sin(tilt) * sin_gamma)
    kz_out = kz_in - qz

    return kz_in, kz_out
    
def QxQyQz_to_gisans_k(qx,qy,qz,wavelength,incident_angle):
    '''
    Overview:
        A Python implementation of the Q space to k_in k_out conversion.
    using the assumptions common to GISANS:
    the incoming beam is fixed in angle at itheta_inplane = 0, 
    itheta_outofplane = constant
    and the outgoing beam is just referenced to that!

    Parameters:

    '''
    
    qx = asarray(qx)
    qy = asarray(qy)
    qz = asarray(qz)
    outshape = (size(qx), size(qy), size(qz))
    wavelength = asarray(wavelength)

    if qy == None: qy = asarray([0.0])

    Qmag = sqrt((qx**2) + (qy**2) + (qz**2))

    if Qmag.size == 1:
        if Qmag==0.0: Qmag = 1.0
    else:
        Qmag[Qmag == 0.0] = 1.0

    k0 = 2.0*pi/wavelength
    kz_in = ones(outshape) * k0 * sin(incident_angle)
    
    # still have to fix this...
    twoth = 2.0 * arcsin(Qmag/(2.0*k0))
    tilt = arctan2(qx,qz)

    th_in = (twoth/2.0) + tilt
    th_out = (twoth/2.0) - tilt

    kz_in = k0 * sin(th_in)
    kz_out = -k0 * sin(th_out)

    return kz_in, kz_out

def QxQyQz_to_angle(space, alphai, intensity, wavelength):  
    
    '''
    Overview:
        A Python implementation of the Q space to real (angular) space conversion.
        
    '''
    
    # convert angle to radians
    alphai = alphai * (pi / 180)
    
    # determine wave vector (k)
    kvec = 2.0*pi/wavelength
    
    x = space.q_list[0]
    y = space.q_list[1]
    z = space.q_list[2]
    
    # upper and lowerbounds for reflected angle
    alphaf_min = alphai
    alphaf_max = 25 * alphai

    # upper and lowerbounds for in-plane angle 
    iptheta_max = arcsin((y[y.argmax()] / kvec))
    iptheta_min = -iptheta_max
    
    # grab equally spaced intervals between upper and lowerbound angles
    alphaf = linspace(alphaf_min, alphaf_max, size(x))
    iptheta = linspace(iptheta_min, iptheta_max, size(y))
    
    # calculate the vertical axis q values
    qx = kvec * (cos( alphaf ) - cos(  alphai) )
    qz = kvec * (sin( alphaf ) - sin(  alphai) )
    
    # calculate the horizontal axis q values
    qy = kvec * sin( iptheta )

    xvals = x.copy()
    yvals = y.copy()
    zvals = z.copy()
    
    # Finds closest x, y, and z values for each of the qx, qy, and qz values
    for i in range(size(qx)):
        closestx = x[0]
        closesty = y[0]
        closestz = z[0]

        for ii in range(size(x)):
            if abs(x[ii] - qx[i]) < abs(closestx - qx[i]): 
                closestx = x[ii]
                xvals[i] = ii
                
        for ii in range(size(y)):
            if abs(y[ii] - qy[i]) < abs(closesty - qy[i]): 
                closesty = y[ii]
                yvals[i] = ii
               
        for ii in range(size(z)):
            if abs(z[ii] - qz[i]) < abs(closestz - qz[i]): 
                closestz = z[ii]
                zvals[i] = ii
    
    xvals = xvals.astype(int)
    yvals = yvals.astype(int) 
    zvals = zvals.astype(int)
    
    # Converts axis values to degrees
    y_values = degrees(iptheta)
    z_values = degrees(alphaf)
    
    angular_intensity = np.empty((size(qx),size(qz)))
    
    # Pulls intensity from the born approximation array and stores it in the correct
    # angle space
    for i in range(size(iptheta)):
        for j in range(size(alphaf)): 
            angular_intensity[i][j] = intensity[xvals[5]][yvals[i]][zvals[j]]
         
    return angular_intensity, y_values, z_values

    
def complete_formula(czt_result, step, Q):
    '''
    the scattering from the unit cell is not just the fft of the unit cell:
    the scattering is the integral of e^{i q r} * rho(r), which can be
    approximated by a Fourier sum if rho(r) is taken to be constant over the 
    step size in r, then the integral over a single block dr can be calculated
    exactly.  The blocks are then all summed to give the complete integral - 
    that last summation is equivalent to the FFT mathematically.  The constant
    integration factor for a block is pulled out of the sum and is calculated 
    below...  and then multiplied by the FFT.

    czt_results = the resulting Fourier transform of the unit cell for the Q
    values requested by the user
    step = the real space step size represented by the discretized units of the
    unit cell
    '''
    if Q.qx_refract == None:
        qx = Q.q_list[0].reshape(size(Q.q_list[0]),1,1)
    else:
        qx = Q.qx_refract

    qy = Q.q_list[1].reshape(1,size(Q.q_list[1]),1)
    qz = Q.q_list[2].reshape(1,1,size(Q.q_list[2]))

    x_comp = zeros_like(qx)
    y_comp = zeros_like(qy)
    z_comp = zeros_like(qz)
    
    delta_x = step[0]
    delta_y = step[1]
    delta_z = step[2]
        
    x_comp[qx != 0.0] = (( 1. - exp(1j*qx[qx != 0]*delta_x))*(-1j / qx[qx != 0]))
    y_comp[qy != 0.0] = (( 1. - exp(1j*qy[qy != 0]*delta_y))*(-1j / qy[qy != 0]))
    z_comp[qz != 0.0] = (( 1. - exp(1j*qz[qz != 0]*delta_z))*(-1j / qz[qz != 0]))
        
    x_comp[qx == 0.0] = delta_x
    y_comp[qy == 0.0] = delta_y
    z_comp[qz == 0.0] = delta_z

    #Creates individual values of the delta_realspace so that the formalism is
    #easier to read.

    czt_result *= x_comp * y_comp * z_comp

    return czt_result


def trip_czt(unit,stepSpace,q_points,q_mins,q_maxs):
    '''
    Given a 3d unit with a set of three limits, this solves the 3d fft using the
    Chirp-z transform.

    unit = The matrix representation of the unit cell
    
    q_points = the number of points that q will be solved for in the x,y and z
    direction
    
    q_mins = the minimum q values for the x, y, and z direction
    
    q_max = the maximum q values for the x, y, and z direction
    
    '''

    frequency = (2*pi)/stepSpace

    intermediate_matrix_one = czt.zoomfft(unit,q_mins[0],q_maxs[0],q_points[0],
                                          frequency[0],axis = 0)

    intermediate_matrix_two = czt.zoomfft(intermediate_matrix_one,q_mins[1],
                                          q_maxs[1],q_points[1],frequency[1],
                                          axis=1)

    czt_result = czt.zoomfft(intermediate_matrix_two,q_mins[2],
                             q_maxs[2],q_points[2],frequency[2],axis=2)

    return czt_result



def xy_czt(unit,step,q_points,q_mins,q_max):
    '''
    Given a 3d unit with a set of three limits, this solves the 2d fft using the
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

class smbaWavefunction:

    def __init__(self, kz_in, array_of_sld):
        from numpy import exp,shape

        self.kz_in = kz_in
        self.array_of_sld = array_of_sld

        layer_num_total = array_of_sld.shape[0]
        self.layer_num_total = layer_num_total
        self.total_thickness = sum(array_of_sld[1:-1,1])

        SLD_incident = array_of_sld[0,0]
        SLD_substrate = array_of_sld[-1,0]

        nz =[None]*layer_num_total

        k0z = sqrt(array(kz_in**2 + 4 * pi * SLD_incident,dtype = 'complex'))

        nz[0] = sqrt( complex(1) - 4 * pi * SLD_incident / k0z**2 )

        nz[-1] = sqrt( complex(1) - 4 * pi * SLD_substrate / k0z**2 )

        B11 = ones(shape(kz_in),dtype='complex')
        B22 = ones(shape(kz_in),dtype='complex')
        B21 = zeros(shape(kz_in),dtype='complex')
        B12 = zeros(shape(kz_in),dtype='complex')

        for layer_num in range(1, layer_num_total-1):

            #leaving off the incident medium and substrate from sum
            SLD,thickness,mu = array_of_sld[layer_num]

            nz[:][layer_num] = sqrt(complex( 1 - 4 * pi * SLD/ k0z**2 ))
            # THIS WAS INCORRECT:
            # I believe that this should be kz = nz * k0z,
            # not kz = nz * kz_in.
            # This was a mistake carried over from my wavefunction.py library, 
            # which was also incorrect. - BBM 2/10/2012
            #kz = nz[:][layer_num] * kz_in
            kz = nz[:][layer_num] * k0z
            n = nz[:][layer_num]

            M11 = cos(kz * thickness)
            M12 = 1/n * sin(kz * thickness)
            M21 = -n * sin(kz * thickness)
            M22 = cos(kz * thickness)

            C1 = B11*M11 + B21*M12
            C2 = B11*M21 + B21*M22
            B11 = C1
            B22 = C2

        self.nz = nz

        r = (B11 + (1j * nz[:][0] * B12) + (1/(1j * nz[:][-1])*(
            -B21 - 1j * nz[:][0] * B22))) / ((-B11 + 1j * nz[:][0] * B12) + (
                             1/(1j * nz[:][-1])*( B21 - 1j * nz[:][0] * B22)))

        self.r = r

        self.t = zeros(shape(r),dtype = 'complex')
        self.t[nz[:][-1].real != 0.0] = 1.0 + self.r[nz[:][-1].real != 0.0]

        self.kz_transmitted = nz[:][-1] * k0z

        return None

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
    import osrefl.model.sample_prep
    from osrefl.model.sample_prep import Parallelapiped, Layer, Scene, GeomUnit, Rectilinear, Beam, Q_space
    #from scatter import *
    #from pylab import *
    #from osrefl.loaders.omfLoader import Omf
    #magneticBA test
    import os
    DATA_PATH = os.path.join('..', os.path.dirname(osrefl.__file__), 'examples', 'data')
    #print open(DATA_PATH).read()
    import wx
    filepath = wx.FileSelector('test.omf')
    mag = Omf(filepath)
    #mag.view()



    Au = (Parallelapiped(SLD = 4.506842e-6,
                                     Ms = 8.6e5,dim=[5.0e4,5.0e4,2.0e4]))

    Cr = (Layer(SLD = 3.01e-6,Ms = 8.6e5,
                            thickness_value = 1000.0))

    Au.on_top_of(Cr)
    scene = Scene([Au,Cr])

    GeoUnit = (GeomUnit(Dxyz = [10.0e4,10.0e4,2.5e4],
                                    n = [70,70,30],scene = scene))

    unit = GeoUnit.buildUnit()
    unit.generateMIF(filename='test.mif')
    #unit.viewSlice()
    space = Q_space([-.0001,-0.001,0.00002],[.0001,0.001,0.04],[50,5,50])
    lattice = Rectilinear([20,20,1],unit)

    beam = Beam(5.0,None,None,0.05,None)

    from osrefl.theory.scatter import Calculator
    magtest = Calculator(lattice, beam, space, unit, mag)

    magtest.magneticBA()
    magtest.view_uncorrected()
    show()
    
    print 'gpu result:'
    interpWaveCalc(unit,space,beam, proc='gpu')
    print 'cpu result:'
    return space, interpWaveCalc(unit,space,beam, proc='cpu')
    

if __name__=="__main__":_test()
