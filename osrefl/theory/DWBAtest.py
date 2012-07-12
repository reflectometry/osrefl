# Copyright (C) 2008 University of Maryland
# All rights reserved.
# See LICENSE.txt for details.
# Author: Christopher Metting

#Starting Date:6/12/2009

from pylab import *
from numpy import *
from time import time
from  ..model.sample_prep import Q_space
from .approximations import wavefunction_format

def DWBA_form(cell,lattice,beam,q,refract = True):
    '''
    The scattering is calculated in scatCalc because we need to open up the
    possibility for qx refraction on the interpolation.
    '''

    if refract == True:

        from scipy.interpolate import interp1d

        scat = zeros(q.points, dtype = 'complex')
        qvec = q.vectorize()

        q.getKSpace(beam.wavelength)
        qx_refract = qvec[0].repeat(q.points[1],axis=1)
        qx_refract = qx_refract.repeat(q.points[2],axis=2)

        qx_refract[q.kin <= 0.0] += beam.wavelength*cell.inc_sub[1,0]
        qx_refract[q.kout >= 0.0] -= beam.wavelength*cell.inc_sub[1,0]

        q.qx_refract = qx_refract

        qxMinTemp = qx_refract.min()-3*q.q_step[0]
        qxMaxTemp = qx_refract.max()+3*q.q_step[0]

        #doubles  the interpolation q for a more accurate interpolation
        newX = arange(qxMinTemp,qxMaxTemp,q.q_step[0]/2.0)

        newQ = Q_space([qxMinTemp,q.minimums[1],q.minimums[2]],
                       [qxMaxTemp,q.maximums[1],q.maximums[2]],
                       [size(newX),q.points[1],q.points[2]])

        largScat = scatCalc(cell,lattice,beam,newQ)

        for ii in range (size(q.q_list[1])):
            for iii in range(size(q.q_list[2])):
                realSplineFunc = interp1d(newQ.q_list[0],largScat.real[:,ii,iii])
                imagSplineFunc = interp1d(newQ.q_list[0],largScat.imag[:,ii,iii])

                interpReal = realSplineFunc(qx_refract[:,ii,iii])
                interpImag = imagSplineFunc(qx_refract[:,ii,iii])

                scat[:,ii,iii].real = interpReal
                scat[:,ii,iii].imag = interpImag

    else:
        scat = scatCalc(cell,lattice,beam,q)
    '''
    imshow(log10(rot90(sum(((abs(scat)**2)).real,axis=1))), extent = q.getExtent(), aspect = 'auto')
    show()
    '''
    return(scat)

def scatCalc(cell,lattice,beam,q):
    '''
    Math from Kentzinger et al. in Physical Review B, 77, 1044335(2008)
    '''
    
    #Front of Eq (20)
    m = 1.674e-27
    h_bar = 6.62607e-14

    Vfac = -m/(2*pi*h_bar**2)

    q.getKSpace(beam.wavelength)
    
    wfc = zeros(q.points,dtype = 'complex')
    wfd = zeros(q.points,dtype = 'complex')
    
    pio = [None]*cell.n[2]
    pit = [None]*cell.n[2]
    poo = [None]*cell.n[2]
    pot = [None]*cell.n[2]
    
    SLDArray = wavefunction_format(cell.unit, cell.step[2], absorbtion = None)
    
    for i in range(size(q.q_list[0])):
        print 'qx number: ', i, ' calculating'

        for ii in range(size(q.q_list[1])):
            
            poskiWavePar = dwbaWavefunction(q.kin[i,ii,:],SLDArray)
            negkfWavePar = dwbaWavefunction(-q.kout[i,ii,:],(SLDArray))
            
            pio = poskiWavePar.c
            pit = poskiWavePar.d
            poo = negkfWavePar.c
            pot = negkfWavePar.d

            wfc[i,ii,:] = sum(pio + poo)
            wfd[i,ii,:] = sum(pit + pot)
    
    return (wfc * wfd)

class dwbaWavefunction:

    def __init__(self, kz, SLDArray):

        self.kz = kz
        self.SLDArray = SLDArray

        self.layerCount = SLDArray.shape[0]
        self.thickness = sum(SLDArray[1:-1,1])

        SLD_inc = SLDArray[0,0]
        SLD_sub = SLDArray[-1,0]

        B11 = ones(shape(kz),dtype='complex')
        B22 = ones(shape(kz),dtype='complex')
        B21 = zeros(shape(kz),dtype='complex')
        B12 = zeros(shape(kz),dtype='complex')

        M11 = [None]*self.layerCount
        M12 = [None]*self.layerCount
        M21 = [None]*self.layerCount
        M22 = [None]*self.layerCount

        Bl11 = [None]*self.layerCount
        Bl12 = [None]*self.layerCount
        Bl21 = [None]*self.layerCount
        Bl22 = [None]*self.layerCount

        Bl11[0] = B11
        Bl12[0] = B22
        Bl21[0] = B21
        Bl22[0] = B12

        self.c = [None]*self.layerCount
        self.d = [None]*self.layerCount

        nz =[None]*self.layerCount

        k0z = sqrt(asarray(kz**2 + 4 * pi * SLD_inc,dtype = 'complex'))

        nz[0] = sqrt( complex(1) - 4 * pi * SLD_inc / k0z**2 )
        nz[-1] = sqrt( complex(1) - 4 * pi * SLD_sub / k0z**2 )

        for l in range(1, self.layerCount-1):

            #leaving off the incident medium and substrate from sum
            SLD,thickness,mu = self.SLDArray[l]

            nz[l] = sqrt(complex(1) - 4 * pi * SLD/ k0z**2 )
            kzl =( nz[l] * k0z ) # edit: BBM 02/10/2012
            n = nz[l]

            M11[l] = asarray(cos(kzl * thickness),dtype = 'complex')
            M12[l] = asarray(1/n * sin(kzl * thickness),dtype = 'complex')
            M21[l] = asarray((-n) * sin(kzl * thickness),dtype = 'complex')
            M22[l] = asarray(cos(kzl * thickness),dtype = 'complex')

            C1 = B11*M11[l] + B21*M12[l]
            C2 = B11*M21[l] + B21*M22[l]
            B11 = C1
            B21 = C2


            C1 = B12*M11[l] + B22*M12[l]
            C2 = B12*M21[l] + B22*M22[l]
            B12 = C1
            B22 = C2

            Bl11[l] = B11
            Bl21[l] = B21
            Bl12[l] = B12
            Bl22[l] = B22

        self.kz_l = nz * k0z

        r = (B11 + (1j * nz[0] * B12) + (1/(1j * nz[-1])*(
            -B21 - 1j * nz[0] * B22))) / ((-B11 + 1j * nz[0] * B12) + (
                             1/(1j * nz[-1])*( B21 - 1j * nz[0] * B22)))

        Bl11[-1] = ones(shape(kz))
        Bl12[-1] = zeros(shape(kz))
        Bl21[-1] = ones(shape(kz))
        Bl22[-1] = zeros(shape(kz))

        self.r = r

        self.t = zeros(shape(r),dtype = 'complex')
        self.t[nz[-1].real != 0.0] = 1.0 + self.r[nz[-1].real != 0.0]

        self.c[0] = ones(shape(kz),dtype='complex') # incident beam has intensity 1
        self.d[0] = r # reflected beam has intensity |r|**2

        p = asarray(1.0 + r,dtype ='complex') #psi
        pp = asarray(1j * nz[0] * (1 - r),dtype='complex') #psi prime / k0z

        M11[0] = ones(shape(kz),dtype='complex')
        M12[0] = ones(shape(kz),dtype='complex')
        M21[0] = ones(shape(kz),dtype='complex')
        M22[0] = ones(shape(kz),dtype='complex')

        M11[-1] = zeros(shape(kz),dtype='complex')
        M12[-1] = ones(shape(kz),dtype='complex')
        M21[-1] = ones(shape(kz),dtype='complex')
        M22[-1] = zeros(shape(kz),dtype='complex')

        z_interface = 0.0

        for l in range(1,self.layerCount-1):
            ## this algorithm works all the way into the substrate

            pForDot = copy(p)
            ppForDot = copy(pp)

            p = (M11[l]*pForDot) + (M12[l]*ppForDot)
            pp = (M21[l]*pForDot) + (M22[l]*ppForDot)


            #Fine, This is c and d
            self.c[l] = (.5* exp(-1j*k0z*nz[l]*(z_interface))*
                         (p + (pp/(1j*nz[l]))))
            self.d[l] = (.5* exp(1j*k0z*nz[l]*(z_interface))*
                         (p - (pp/(1j*nz[l]))))

            z_interface += thickness

        # fill final c,d

        self.c[-1] = self.t
        self.d[-1] = zeros(shape(kz),dtype='complex')
        return

