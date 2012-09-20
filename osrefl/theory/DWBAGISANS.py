# Copyright (C) 2008 University of Maryland
# All rights reserved.
# See LICENSE.txt for details.
# Author: Christopher Metting

#Starting Date:6/12/2009

from pylab import *
from numpy import *
import numpy as np
from time import time
from  ..model.sample_prep import Q_space
from .approximations import wavefunction_format

SIMULATE_BA = False

def DWBA_form(cell,lattice,beam,q, angle_in):

    scat = scatCalc(cell,lattice,beam,q, angle_in)

    return scat[0], scat[1], scat[2]

def scatCalc(cell,lattice,beam,q, angle_in):
    
    '''
    Math from Kentzinger et al. in Physical Review B, 77, 1044335(2008)
    
    '''
    
    #Front of Eq (20)
    m = 1.674e-27
    h_bar = 6.62607e-14

    Vfac = -m/(2*pi*h_bar**2)

    wavelength = beam.wavelength

    # convert angle to radians
    angle_in = angle_in * (pi / 180)
    
    # determine wave vector (k)
    kvec = 2.0*pi/wavelength

    # upper and lowerbounds for reflected angle
    alphaf_min = angle_in
    alphaf_max = 100 * angle_in

    # upper and lowerbounds for in-plane angle 
    iptheta_max = arcsin((q.q_list[1][q.q_list[1].argmax()] / kvec))
    iptheta_min = -iptheta_max
    
    # grab equally spaced intervals between upper and lowerbound angles
    angle_out = linspace(alphaf_min, alphaf_max, size(q.q_list[2]))
    iptheta = linspace(iptheta_min, iptheta_max, size(q.q_list[1]))
    #alphai = angle_in
    #angle_in = np.zeros_like(angle_out)
    #angle_in.fill(alphai)

    kz_out = kvec * sin( angle_out )
    kx_out = kvec * cos( angle_out )
    ky_out = -kvec * sin( iptheta )
    
    kz_in = kvec * sin( angle_in )
    kx_in = kvec * cos( angle_in )
    ky_in = zeros_like( ky_out )

    scat = zeros((q.points[1], q.points[2]), dtype = 'complex')
        
    # PSI in one
    # PSI in two
    # PSI out one
    # PSI out two
    
    pio = [None]*cell.n[2]
    pit = [None]*cell.n[2]
    poo = [None]*cell.n[2]
    pot = [None]*cell.n[2]

    pil = [None]*cell.n[2]
    pfl = [None]*cell.n[2]

    q_piopoo = [None]*cell.n[2]
    q_piopot = [None]*cell.n[2]
    q_pitpoo = [None]*cell.n[2]
    q_pitpot = [None]*cell.n[2]

    x = cell.value_list[0].reshape((cell.n[0],1,1))
    y = cell.value_list[1].reshape((1,cell.n[1],1))
    z = cell.value_list[2].reshape((1,1,cell.n[2]))
    
    #Averages the in-plane scattering length density and formats the new
    #object as [SLD,Thickeness,Absorbtion] for each z layer
    SLDArray = wavefunction_format(cell.unit, cell.step[2], absorbtion = None)
    
    #This is the calculation of the critical edge. It is needed for the
    #calculation of p.
    pcl = sqrt(4*pi*SLDArray[:,0])
    
    #The cell is originally oriented so that the the bottom of the unit cell
    #is located at the origin. This flips the cell so that the stack is ordered
    #in the opposite direction.
    flipCell = zeros(shape(cell.unit))

    for i in range(cell.n[2]):
        flipCell[:,:,i] = cell.unit[:,:,shape(cell.unit)[2]-i-1]

    #This calculates the residual potential by taking the difference between
    #the reference potential and the actual potential
    Vres = flipCell - (SLDArray[:,0]).reshape((1,1,cell.n[2]))
    
    #This is the rho used in eq. 20. The integration is the residual potential
    #relative to the reference potential.
    rhoTilOverRho = Vres/(SLDArray[:,0]).reshape((1,1,cell.n[2]))
    rhoTilOverRho[isnan(rhoTilOverRho)] = 0.0

    #The next few lines calculate the c and d values for each layer.
    #This is done by calculating the specular reflectivity and then
    #tracing the final reflected intensity back into the sample.
       
    for i in range(size(iptheta)):
        print 'iptheta:', degrees(iptheta[i]), 'calculating (', i+1, 'of', size(iptheta), ')' 
            
        for ii in range(size(angle_out)):
            if SIMULATE_BA:
                pio = ones((cell.n[2],), dtype='complex') 
                pit = zeros((cell.n[2],), dtype='complex')
                poo = zeros((cell.n[2],), dtype='complex')
                pot = ones((cell.n[2],), dtype='complex')
            else:  
                poskiWavePar = dwbaWavefunction(kz_in,SLDArray)
                negkfWavePar = dwbaWavefunction(kz_out[ii],SLDArray)
                 
                pio = poskiWavePar.c
                pit = poskiWavePar.d
                poo = negkfWavePar.c
                pot = negkfWavePar.d
        
            for l in range(cell.n[2]):
                    
                #Solves the equation shown after eq. 11 on page 5.
                pil[l]=sqrt(asarray((kz_in**2)-(pcl[l]**2),
                                        dtype = 'complex'))
                pfl[l]=sqrt(asarray((kz_out[ii]**2)-(pcl[l]**2),
                                        dtype = 'complex'))
                
                #Equations directly after eq (18).
                q_piopoo[l] = -pfl[l] - pil[l]
                q_piopot[l] = -pfl[l] + pil[l]
                q_pitpoo[l] = pfl[l] - pil[l]
                q_pitpot[l] = pfl[l] + pil[l]

            pil = asarray(pil)
            pfl = asarray(pfl)
        
            q_piopoo = asarray(q_piopoo)
            q_piopot = asarray(q_piopot)
            q_pitpoo = asarray(q_pitpoo)
            q_pitpot = asarray(q_pitpot)
        
            pio = asarray(pio)
            pit = asarray(pit)
            poo = asarray(poo)
            pot = asarray(pot)
    
            ######## 
            # EDIT: bbm 07/20/2012
            # this is not Eq. 18, which refers only to the out-of-plane (z) Laue factor
            # this is the necessary Laue factor to do the integral in eq. 20
            # as a finite sum over blocks of constant rho in the x-y plane
            ########f (mask.all() != False): 
            qx = kx_in - kx_out[ii]
            if qx != 0:
                laux = ((-1j / qx) * (exp(1j * qx * cell.step[0]) - 1.0))
            else:
                laux = complex(cell.step[0])
                
            qy = -ky_out[i]
            if qy != 0:
                lauy = ((-1j / qy) * (exp(1j * qy * cell.step[1]) - 1.0))
            else:
                lauy = complex(cell.step[1])                 
            
            #Eq. 20 (including only rhoN - rhoM is assumed to be zero)
            ftwRef = (Vfac*sum(sum(rhoTilOverRho * exp(1j*qx*x)*
                       exp(1j*qy*y),axis=0),axis=0))
            
            # finite-sum corrections for the x and y directions
            ftwRef *= laux
            ftwRef *= lauy
    
            #Eq. 19
            ftwRef = ((SLDArray[:,0]).reshape((1,1,cell.n[2]))*
                      ftwRef.reshape((1,1,cell.n[2])))
              
            ft = ftwRef.copy()

            pioSel = pio.reshape((1,1,cell.n[2]))
            pitSel = pit.reshape((1,1,cell.n[2]))
            pooSel = poo.reshape((1,1,cell.n[2]))
            potSel = pot.reshape((1,1,cell.n[2]))

            q_piopoo_sel = q_piopoo.reshape((1,1,cell.n[2]))
            q_piopot_sel = q_piopot.reshape((1,1,cell.n[2]))
            q_pitpoo_sel = q_pitpoo.reshape((1,1,cell.n[2]))
            q_pitpot_sel = q_pitpot.reshape((1,1,cell.n[2]))

            pil_sel = pil.reshape((1,1,cell.n[2]))
            pfl_sel = pfl.reshape((1,1,cell.n[2]))

            #equation 15
            scat_PioPoo = (pioSel * exp(1j*pil_sel*z)*ft*
                           exp(1j*pfl_sel*z) * pooSel)
            scat_PioPot = (pioSel * exp(1j*pil_sel*z)*ft*
                           exp(-1j*pfl_sel*z)*potSel)
            scat_PitPoo = (pitSel * exp(-1j*pil_sel*z)*ft*
                           exp(1j*pfl_sel*z) *pooSel)
            scat_PitPot = (pitSel * exp(-1j*pil_sel*z)*ft*
                           exp(-1j*pfl_sel*z)* potSel)
            
            mask = (q_piopoo_sel != 0)
            scat_PioPoo[mask] *= ((-1j / q_piopoo_sel[mask]) * 
                            (exp(1j *q_piopoo_sel[mask] * cell.step[2]) - 1.0))
            scat_PioPoo[q_piopoo_sel == 0] *= cell.step[2]
            
            mask = (q_piopot_sel != 0)
            scat_PioPot[mask] *= ((-1j / q_piopot_sel[mask]) * 
                            (exp(1j *q_piopot_sel[mask] * cell.step[2]) - 1.0))
            scat_PioPot[q_piopot_sel == 0] *= cell.step[2]
            
            mask = (q_pitpoo_sel != 0)
            scat_PitPoo[mask] *= ((-1j / q_pitpoo_sel[mask]) * 
                            (exp(1j *q_pitpoo_sel[mask] * cell.step[2]) - 1.0))
            scat_PitPoo[q_pitpoo_sel == 0] *= cell.step[2]
            
            mask = (q_pitpot_sel != 0)
            scat_PitPot[mask] *= ((-1j / q_pitpot_sel[mask]) * 
                            (exp(1j *q_pitpot_sel[mask] * cell.step[2]) - 1.0))
            scat_PitPot[q_pitpot_sel == 0] *= cell.step[2]
            
            #Exactly equation15
            scat[i, ii]= sum(scat_PioPoo + scat_PioPot + scat_PitPoo + scat_PitPot)
            
    xvals = degrees(iptheta)
    yvals = degrees(angle_out)

    return scat, xvals, yvals

class dwbaWavefunction:

    def __init__(self, kz, SLDArray):
        
        if not isinstance(kz, ndarray):
            kz = array([kz], dtype=complex)
        #kz = array([kz]).flatten().astype(complex)
        self.kz = kz
        kzlen = kz.shape
        sldlen = len(SLDArray)
        self.SLDArray = SLDArray
        self.r = zeros(kzlen, dtype=complex)
        self.kz_l = zeros((sldlen,) + kzlen, dtype=complex)
        self.c = zeros((sldlen,) + kzlen, dtype=complex)
        self.d = zeros((sldlen,) + kzlen, dtype=complex)
        
        neg_k_mask = (self.kz < 0)
        pos_k_mask = logical_not(neg_k_mask)
        kz_ln, cn, dn, rn = self.calc_r_cd(self.kz[neg_k_mask], kz_neg=True)
        self.kz_l[:,neg_k_mask] = kz_ln
        self.c[:,neg_k_mask] = cn
        self.d[:,neg_k_mask] = dn
        self.r[neg_k_mask] = rn
        
        kz_l, c, d, r = self.calc_r_cd(self.kz[pos_k_mask], kz_neg=False)
        self.kz_l[:,pos_k_mask] = kz_l
        self.c[:,pos_k_mask] = c
        self.d[:,pos_k_mask] = d
        self.r[pos_k_mask] = r
        
    
    def calc_r_cd(self, kz, kz_neg=False):
        if kz_neg==True:
            workingSLD = self.SLDArray[::-1]
        else:
            workingSLD = self.SLDArray
            
        layerCount = workingSLD.shape[0]
        thickness = sum(workingSLD[1:-1,1])

        SLD_inc = workingSLD[0,0] + 1j * workingSLD[0,2]
        SLD_sub = workingSLD[-1,0] + 1j * workingSLD[-1,2]

        B11 = ones(shape(kz),dtype='complex')
        B22 = ones(shape(kz),dtype='complex')
        B21 = zeros(shape(kz),dtype='complex')
        B12 = zeros(shape(kz),dtype='complex')

        M11 = [None]*layerCount
        M12 = [None]*layerCount
        M21 = [None]*layerCount
        M22 = [None]*layerCount

        Bl11 = [None]*layerCount
        Bl12 = [None]*layerCount
        Bl21 = [None]*layerCount
        Bl22 = [None]*layerCount

        Bl11[0] = B11
        Bl12[0] = B22
        Bl21[0] = B21
        Bl22[0] = B12

        c = [None]*layerCount
        d = [None]*layerCount
        nz =[None]*layerCount

        k0z = sqrt(asarray(kz**2 + 4 * pi * SLD_inc,dtype = 'complex'))# always positive!

        nz[0] = sqrt( complex(1) - 4 * pi * SLD_inc / k0z**2 )
        nz[-1] = sqrt( complex(1) - 4 * pi * SLD_sub / k0z**2 )

        for l in range(1, layerCount-1):

            #leaving off the incident medium and substrate from sum
            SLD,thickness,mu = workingSLD[l]

            nz[l] = sqrt(complex(1) - 4 * pi * (SLD+1j*mu)/ k0z**2 )
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

        kz_l = nz * k0z

        r = (B11 + (1j * nz[0] * B12) + (1/(1j * nz[-1])*(
            -B21 - 1j * nz[0] * B22))) / (-B11 + (1j * nz[0] * B12) + (
                             1/(1j * nz[-1])*( B21 - 1j * nz[0] * B22)))

        Bl11[-1] = ones(shape(kz))
        Bl12[-1] = zeros(shape(kz))
        Bl21[-1] = ones(shape(kz))
        Bl22[-1] = zeros(shape(kz))

        #self.r = r

        self.t = zeros(shape(r),dtype = 'complex')
        #self.t[nz[-1].real != 0.0] = 1.0 + self.r[nz[-1].real != 0.0]

        c[0] = ones(shape(kz),dtype='complex') # incident beam has intensity 1
        d[0] = r # reflected beam has intensity |r|**2

        p = asarray(1.0 + r,dtype ='complex') #psi
        pp = asarray(1j * kz_l[0] * (1 - r),dtype='complex') #psi prime

        #M11[0] = ones(shape(kz),dtype='complex')
        #M12[0] = ones(shape(kz),dtype='complex')
        #M21[0] = ones(shape(kz),dtype='complex')
        #M22[0] = ones(shape(kz),dtype='complex')

        z_interface = 0.0

        for l in range(1,layerCount-1):
            ## this algorithm works all the way into the substrate
            SLD,thickness,mu = workingSLD[l]
            pForDot = copy(p)
            ppForDot = copy(pp)
            
            #Fine, This is c and d
            kzl =( nz[l] * k0z ) 

            c[l] = (.5* exp(-1j*kzl*(z_interface))*(p + (pp/(1j*kzl))))
            d[l] = (.5* exp( 1j*kzl*(z_interface))*(p - (pp/(1j*kzl))))
            ## Moved ^ above v to model wavefunction.js WRT 7/16/12
            
            p = (M11[l]*pForDot) + (M12[l]*ppForDot/k0z)
            pp = (k0z*M21[l]*pForDot) + (M22[l]*ppForDot) 

            z_interface += thickness

        # fill final c,d
        l=-1
        kzl =( nz[l] * k0z ) 
        c[l] = (.5* exp(-1j*kzl*(z_interface))*(p + (pp/(1j*kzl))))
        #self.c[-1] =  (.5* exp(-1j*kzl*(z_interface))*(p + (pp/(1j*kzl))))
        d[-1] = zeros(shape(kz),dtype='complex')
                
        if kz_neg==True:
            print "neg_kz!"
            return [-kz_l[::-1], d[::-1], c[::-1], r[::-1]]
        else:
            return [kz_l, c, d, r]

def _test():
    # run from ipython by starting in root osrefl directory, 
    # from osrefl.theory.DWBA import _test
    # test()
    # ...
    from osrefl.model.sample_prep import Parallelapiped, Layer, Scene, GeomUnit, Rectilinear, Beam
    
    Au = Parallelapiped(SLD = 4.506842e-6,dim=[3.7e4,3.7e4,630.0])#, curve = .56)
    Cr = Layer(SLD = 3.01e-6,thickness_value = 48.0)

    #Au.on_top_of(Cr)

    #scene = Scene([Au,Cr])
    scene = Scene([Au])

    GeoUnit = GeomUnit(Dxyz = [10.0e4,10.0e4,700.0], n = [20,21,40],
                       #scene = scene, inc_sub = [0.0,0.0])
                       scene = scene, inc_sub = [0.0,2.07e-6])
    unit = GeoUnit.buildUnit()
    unit.add_media()

    lattice = Rectilinear([20.0,20.0,1.0],unit)
    beam = Beam(5.0,.02,None,0.02,None)
    q = Q_space([-.0002,-0.002,0.00002],[.0002,.002,0.1],[100,5,150])


    SLDArray = wavefunction_format(unit.unit, unit.step[2], absorbtion = None)

    '''
    kiWavePar = dwbaWavefunction(q.kin,SLDArray)

    test = 2

    bbmTest = neutron_wavefunction(q.kin[test,2,50],SLDArray)

    cCollect = zeros(shape(kiWavePar.c)[0])
    dCollect = zeros(shape(kiWavePar.d)[0])

    c = asarray(kiWavePar.c)
    d = asarray(kiWavePar.d)
    for i in range(shape(kiWavePar.c)[0]):

        temp = kiWavePar.c[i]
        cCollect[i] = temp[test,2,50]
        temp = kiWavePar.d[i]
        dCollect[i] = temp[test,2,50]

        cCollect=c[:,test,2,50]
        dCollect=d[:,test,2,50]

    plot(bbmTest.c,label = 'bbm')
    plot(cCollect,label = 'py')
    legend()
    figure()
    plot(bbmTest.d,label = 'bbm')
    plot(dCollect,label = 'py')
    legend()
    figure()
    diff = abs(bbmTest.c.real-cCollect.real)/((abs(bbmTest.c.real)+abs(cCollect.real))/2.0)
    plot(diff,label = 'diff')
    show()
    '''
    DWBA_form(unit,lattice,beam,q)

if __name__=="__main__": _test()

