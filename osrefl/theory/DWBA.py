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

def print_timing(func):
    def wrapper(*arg):
        t1 = time()
        res = func(*arg)
        t2 = time()
        print '%s took %0.3f ms' % (func.func_name, (t2-t1)*1000.0)
        return res
    return wrapper

@print_timing
def scatCalc(cell,lattice,beam,q):
    '''
    Math from Kentzinger et al. in Physical Review B, 77, 1044335(2008)
    '''
    
    #Front of Eq (20)
    m = 1.674e-27
    h_bar = 6.62607e-14

    Vfac = -m/(2*pi*h_bar**2)

    
    q.getKSpace(beam.wavelength)
    
    scat = zeros(q.points,dtype = 'complex')

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
    
    #calculates the structure factor using the gaussian convolution.
    if lattice != None:
        SF = lattice.gauss_struc_calc(q)

    for i in range(size(q.q_list[0])):
        print 'qx number: ', i, ' calculating (', i+1, 'of', size(q.q_list[0]), ')'  

        for ii in range(size(q.q_list[1])):
            
            #The next few lines calculate the c and d values for each layer.
            #This is done by calculating the specular reflectivity and then
            #tracing the final reflected intensity back into the sample.
           
            poskiWavePar = dwbaWavefunction(q.kin[i,ii,:],SLDArray)
            negkfWavePar = dwbaWavefunction(-q.kout[i,ii,:],SLDArray)
            
            pio = poskiWavePar.c
            pit = poskiWavePar.d
            k_inl =poskiWavePar.kz_l
            poo = negkfWavePar.c
            pot = negkfWavePar.d
            k_outl =negkfWavePar.kz_l

            for l in range(cell.n[2]):
                
                #Solves the equation shown after eq. 11 on page 5.
                pil[l]=sqrt(asarray((q.kin[i,ii,:]**2)-(pcl[l]**2),
                                    dtype = 'complex'))
                pfl[l]=sqrt(asarray((q.kout[i,ii,:]**2)-(pcl[l]**2),
                                    dtype = 'complex'))
                #print pil[l], poskiWavePar.nz[l] * pio[0]
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

            k_inl = asarray(k_inl)
            k_outl = asarray(k_outl)

            ######## 
            # EDIT: bbm 07/20/2012
            # this is not Eq. 18, which refers only to the out-of-plane (z) Laue factor
            # this is the necessary Laue factor to do the integral in eq. 20
            # as a finite sum over blocks of constant rho in the x-y plane
            ########
            qx = q.q_list[0][i]
            if qx != 0:
                laux = ((-1j / qx) * (exp(1j * qx * cell.step[0]) - 1.0))
            else:
                laux = complex(cell.step[0])
                
            qy = q.q_list[1][ii]
            if qy != 0:
                lauy = ((-1j / qy) * (exp(1j * qy * cell.step[1]) - 1.0))
            else:
                lauy = complex(cell.step[1])       
            
            
            #Eq. 20 (including only rhoN - rhoM is assumed to be zero)
            ftwRef = (Vfac*sum(sum(rhoTilOverRho * exp(1j*q.q_list[0][i]*x)*
                       exp(1j*q.q_list[1][ii]*y),axis = 0),axis=0))
            
            # finite-sum corrections for the x and y directions
            ftwRef *= laux
            ftwRef *= lauy
            

            #Eq.18 with the added structure factor.
            if lattice != None:
                ftwRef *=SF[i,ii,0]

            #ftwRef = ftwRef*((lattice.repeat[0]*cell.Dxyz[0]*lattice.repeat[1]*cell.Dxyz[1]))
            #ftwRef = ftwRef*(lattice.repeat[0]*cell.Dxyz[0])\
            
            #Eq. 19
            ftwRef = ((SLDArray[:,0]).reshape((1,1,cell.n[2]))*
                      ftwRef.reshape((1,1,cell.n[2])))
            
            for iii in range(size(q.q_list[2])):

                ft = ftwRef.copy()

                pioSel = pio[:,iii].reshape((1,1,cell.n[2]))
                pitSel = pit[:,iii].reshape((1,1,cell.n[2]))
                pooSel = poo[:,iii].reshape((1,1,cell.n[2]))
                potSel = pot[:,iii].reshape((1,1,cell.n[2]))

                q_piopoo_sel = q_piopoo[:,iii].reshape((1,1,cell.n[2]))
                q_piopot_sel = q_piopot[:,iii].reshape((1,1,cell.n[2]))
                q_pitpoo_sel = q_pitpoo[:,iii].reshape((1,1,cell.n[2]))
                q_pitpot_sel = q_pitpot[:,iii].reshape((1,1,cell.n[2]))

                pil_sel = pil[:,iii].reshape((1,1,cell.n[2]))
                pfl_sel = pfl[:,iii].reshape((1,1,cell.n[2]))
                
                #equation 15
                scat_PioPoo = (pioSel * exp(1j*pil_sel*z)*ft*
                               exp(1j*pfl_sel*z) * pooSel)
                scat_PioPot = (pioSel * exp(1j*pil_sel*z)*ft*
                               exp(-1j*pfl_sel*z)*potSel)
                scat_PitPoo = (pitSel * exp(-1j*pil_sel*z)*ft*
                               exp(1j*pfl_sel*z) *pooSel)
                scat_PitPot = (pitSel * exp(-1j*pil_sel*z)*ft*
                               exp(-1j*pfl_sel*z)* potSel)

                #equation 15
                # edit 7/23/2012, bbm: 
                # the integration over z is taken care of by eq. 17 and 18, 
                # giving the Laue factor - 
                # the mu and nu sum comes out to 1/4 * 4 * g for unpolarized
                # NO - Wait - changing my mind.
                # 
                # looks like Chris was right - the S propagator in eq. 11
                # is for a wavefunction referenced to the boundary of the 
                # current layer, while our c and d are calculated with respect
                # to z = 0 (not z=z_l), so the extra factor of e^{ikz_l} might
                # be necessary.
#                scat_PioPoo = (pioSel * ft * pooSel)
#                scat_PioPot = (pioSel * ft * potSel)
#                scat_PitPoo = (pitSel * ft * pooSel)
#                scat_PitPot = (pitSel * ft * potSel)

                #equation 18
                # edit 7/23/12, bbm:
                # scat_ was incorrectly set to = cell.step[2] for q==0 case,
                # instead of multiplying (should be *= )
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
                scat[i,ii,iii]= sum(scat_PioPoo + scat_PioPot + 
                                    scat_PitPoo + scat_PitPot)
            
                    
    k_spec = q.q_list[2]/2.0
    dwba_spec = dwbaWavefunction(k_spec,SLDArray)

    locx = q.q_list[0].searchsorted(0.0)
    locy = q.q_list[1].searchsorted(0.0)

    #scat[locx,locy,:] = dwba_spec.r
    
    #semilogy(q.q_list[2],(abs(dwba_spec.r)**2))
    #semilogy(q.q_list[2],sum((abs(scat)**2).real,axis=1)[locx+5,:])
    #figure()
    
    return(scat)


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
            -B21 - 1j * nz[0] * B22))) / (-B11 + (1j * nz[0] * B12) + (
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
        pp = asarray(1j * kz[0] * (1 - r),dtype='complex') #psi prime

        M11[0] = ones(shape(kz),dtype='complex')
        M12[0] = ones(shape(kz),dtype='complex')
        M21[0] = ones(shape(kz),dtype='complex')
        M22[0] = ones(shape(kz),dtype='complex')

        #M11[-1] = zeros(shape(kz),dtype='complex')
        #M12[-1] = ones(shape(kz),dtype='complex')
        #M21[-1] = ones(shape(kz),dtype='complex')
        #M22[-1] = zeros(shape(kz),dtype='complex')

        z_interface = 0.0

        for l in range(1,self.layerCount-1):
            ## this algorithm works all the way into the substrate

            pForDot = copy(p)
            ppForDot = copy(pp)
            
            #Fine, This is c and d
            kzl =( nz[l] * k0z ) 
            self.c[l] = (.5* exp(-1j*kzl*(z_interface))*
                         (p + (pp/(1j*kzl))))
            self.d[l] = (.5* exp(1j*kzl*(z_interface))*
                         (p - (pp/(1j*kzl))))
            
            ## Moved ^ above v to model wavefunction.js WRT 7/16/12
            
            p = (M11[l]*pForDot) + (M12[l]*ppForDot/k0z)
            pp = (k0z*M21[l]*pForDot) + (M22[l]*ppForDot) 

            z_interface += thickness

        # fill final c,d

        self.c[-1] = self.t
        self.d[-1] = zeros(shape(kz),dtype='complex')
        self.nz = nz
        return

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

