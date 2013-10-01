from numpy import *
import numpy as np

class bornWavefunction:

    def __init__(self, kz, SLDArray, sigmax = 1e6, sigmaz=1e3):
        
        if not isinstance(kz, ndarray):
            kz = array([kz], dtype=complex)
        #kz = array([kz]).flatten().astype(complex)
        self.kz = kz
        self.sigmax = sigmax
        self.sigmaz = sigmaz
        kzlen = kz.shape
        sldlen = len(SLDArray)
        self.SLDArray = SLDArray
        self.r = zeros(kzlen, dtype=complex)
        self.r = self.calc_r(self.kz)
           
    def calc_r(self, kz):
        qz = kz*2.0
        sld = self.SLDArray
        r = zeros_like(kz)
        zs = cumsum(sld, axis=0)
        for i in range(len(sld)-1):
            for j in range(len(sld)-1):
                zi = zs[i][1] # interface i location
                zj = zs[j][1]
                dsldi = sld[i+1][0] - sld[i][0]
                dsldj = sld[j+1][0] - sld[j][0]
                r += 16*pi**2/qz**4*dsldi*dsldj*cos(qz*(zj-zi))*exp(-(zj-zi)**2/(self.sigmaz**2))
        return sqrt(r)
        
    def calc_r_refract(self, kz):
        qz = kz*2.0
        sld = self.SLDArray
        rho = sld[:,0]
        #qz_l_sq = 4.0*((kz**2) - 4*pi*rho[:,None])
        qz_l = sqrt(4.0*((kz**2) - 4*pi*rho[:,None]))
        #print "qz_l_sq shape:", qz_l_sq.shape
        r = zeros_like(kz)
        zs = cumsum(sld, axis=0)
        #qzs = sqrt(qz_l_sq) * (sld[:,1,None])
        qzs = qz_l * (sld[:,1,None])
        qzs = cumsum(qzs, axis=0)
        for i in range(len(sld)-1):
            for j in range(len(sld)-1):
                qzi = qzs[i] 
                qzj = qzs[j]
                zi = zs[i][1] # interface i location
                zj = zs[j][1]
                #dsldi = sld[i+1][0]/qz_l_sq[i+1] - sld[i][0]/qz_l_sq[i]
                #dsldj = sld[j+1][0]/qz_l_sq[j+1] - sld[j][0]/qz_l_sq[j]
                dsldi = sld[i+1][0]/qz_l[i+1] - sld[i][0]/qz_l[i]
                dsldj = sld[j+1][0]/qz_l[j+1] - sld[j][0]/qz_l[j]
                #print dsldi.shape, qzj.shape, zj.shape
                #print i,j,qzi-qzj, qz*(zi-zj)
                r += 16*pi**2/qz**2*dsldi*dsldj*cos(qzj-qzi)*exp(-(zj-zi)**2/(self.sigmaz**2))
        return sqrt(r)
                
