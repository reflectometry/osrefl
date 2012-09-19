#from greens_thm_form import greens_form_line, greens_form_shape
from greens_thm_form import div_form_shape as greens_form_shape
from numpy import arange, linspace, float64, indices, zeros_like, ones_like, pi, sin, complex128, array, exp, newaxis, cumsum, sum, cos, sin, log, log10, zeros
from osrefl.theory.DWBAGISANS import dwbaWavefunction

class Shape:
    def __init__(self, name):
        self.name = name
        self.points = []
        self.sld = 0.0
        self.sldi = 0.0

class GISANS_problem(object):
    def __init__(self, 
                 sublayers,
                 matrix,
                 front_sld, front_sldi, 
                 back_sld, back_sldi, 
                 wavelength, 
                 qx, qy, qz, 
                 autoFT=True):
        self.sublayers = sublayers
        self.matrix = matrix
        self.front_sld = front_sld
        self.front_sldi = front_sldi
        self.back_sld = back_sld
        self.back_sldi = back_sldi
        self.wavelength = wavelength
        self._qx = qx
        self._qy = qy
        self._qz = qz
        self.update_SLDArray()
        self.alpha_in = None
        self.FTs = []
        if autoFT == True: self.update_FTs()
    
    def get_qx(self):
        return self._qx
    def set_qx(self, value):
        self._qx = value
        self.update_FTs()
    def del_qx(self):
        del self._qx
    qx = property(get_qx, set_qx, del_qx, "I'm the qx property.")
    
    def get_qy(self):
        return self._qy
    def set_qy(self, value):
        self._qy = value
        self.update_FTs()
    def del_qy(self):
        del self._qy
    qy = property(get_qy, set_qy, del_qy, "I'm the qy property.")
    
    def get_qz(self):
        return self._qz
    def set_qz(self, value):
        self._qz = value
        self.update_FTs()
    def del_qz(self):
        del self._qz
    qz = property(get_qz, set_qz, del_qz, "I'm the qz property.")
        
    def update_SLDArray(self):
        SLDArray = [ [self.front_sld, 0, self.front_sldi] ] # [sld.real, thickness, sld.imag]
        for sl in self.sublayers:
            SLDArray.append([sl[1], sl[3], sl[2]])
        SLDArray.append([self.back_sld, 0, self.back_sldi])
        self.SLDArray = array(SLDArray) 
    
    def update_sublayers(self, sublayers):
        self.sublayers = sublayers
        self.update_SLDArray()
        
    def update_FTs(self):
        dFTs = [] # differential = SLD - (avg. SLD)
        FTs = []
        for sl in self.sublayers:
            dFT = zeros((self.qx.shape[0], self.qy.shape[0]), dtype=complex128)
            #FT = zeros((self.qx.shape[0], self.qy.shape[0]), dtype=complex128)
            qx = self.qx
            qy = self.qy
            shapes = sl[0]
            for shape in shapes:
                dFT += greens_form_shape(shape.points, qx, qy) * (shape.sld)
            dFT += greens_form_shape(self.matrix.points, qx, qy) * (self.matrix.sld)
            FT = dFT.copy()
            FTs.append(FT) # do this before subtracting avg. SLD
            dFT += greens_form_shape(self.matrix.points, qx, qy) * (-sl[1]) # subtract FT of average SLD
            dFTs.append(dFT)
        self.FTs = FTs
        self.dFTs = dFTs
        
    def calc_gisans(self, alpha_in, show_plot=True):
        kz_in_0 = 2*pi/self.wavelength * sin(alpha_in * pi/180.0)
        kz_out_0 = kz_in_0 - self.qz

        wf_in = dwbaWavefunction(kz_in_0, self.SLDArray)
        wf_out = dwbaWavefunction(-kz_out_0, self.SLDArray) # solve 1d equation for time-reversed state
        
        kz_in_l = wf_in.kz_l # inside the layers
        kz_in_p_l = -kz_in_l # prime
        kz_out_l = -wf_out.kz_l # inside the layers
        kz_out_p_l = -kz_out_l    # kz_f_prime in the Sinha paper notation

        dz = self.SLDArray[1:-1,1][:,newaxis]
        zs = cumsum(self.SLDArray[1:-1,1]) - self.SLDArray[1,1] # start at zero with first layer
        z_array = array(zs)[:,newaxis]

        qrt_inside = -kz_in_l[1:-1] - kz_out_l[1:-1]
        qtt_inside = -kz_in_l[1:-1] + kz_out_l[1:-1]
        qtr_inside = +kz_in_l[1:-1] + kz_out_l[1:-1]
        qrr_inside = +kz_in_l[1:-1] - kz_out_l[1:-1]
               
        # the overlap is the forward-moving amplitude c in psi_in multiplied by 
        # the forward-moving amplitude in the time-reversed psi_out, which
        # ends up being the backward-moving amplitude d in the non-time-reversed psi_out
        # (which is calculated by the wavefunction calculator)
        # ... and vice-verso for d and c in psi_in and psi_out
        overlap  = wf_out.c[1:-1] * wf_in.c[1:-1] / (1j * qtt_inside) * (exp(1j * qtt_inside * dz) - 1.0)*exp(1j*qtt_inside*z_array)
        overlap += wf_out.d[1:-1] * wf_in.d[1:-1] / (1j * qrr_inside) * (exp(1j * qrr_inside * dz) - 1.0)*exp(1j*qrr_inside*z_array)
        overlap += wf_out.c[1:-1] * wf_in.d[1:-1] / (1j * qtr_inside) * (exp(1j * qtr_inside * dz) - 1.0)*exp(1j*qtr_inside*z_array)
        overlap += wf_out.d[1:-1] * wf_in.c[1:-1] / (1j * qrt_inside) * (exp(1j * qrt_inside * dz) - 1.0)*exp(1j*qrt_inside*z_array)
        self.overlap = overlap
        overlap_BA  = 1.0 / (1j * self.qz) * (exp(1j * self.qz * dz) - 1.0) * exp(1j*self.qz*z_array)
        self.overlap_BA = overlap_BA
        gisans = sum(sum(overlap[:,newaxis,newaxis,:] * array(self.dFTs)[:,:,:,newaxis], axis=0), axis=0) # first over layers, then Qx
        gisans_BA = sum(sum(overlap_BA[:,newaxis,newaxis,:] * array(self.FTs)[:,:,:,newaxis], axis=0), axis=0) 
        extent = [self.qy.min(), self.qy.max(), self.qz.min(), self.qz.max()]
        
        if show_plot == True:
            from pylab import imshow, figure, colorbar
            zmax = max(log10(abs(gisans)**2).max(), log10(abs(gisans_BA)**2).max())
            zmin = min(log10(abs(gisans)**2).min(), log10(abs(gisans_BA)**2).min())
            figure()
            imshow(log10(abs(gisans)**2).T, origin='lower', extent=extent, aspect='auto', vmax=zmax, vmin=zmin)
            colorbar()
            
            figure()
            imshow(log10(abs(gisans_BA)**2).T, origin='lower', extent=extent, aspect='auto', vmax=zmax, vmin=zmin)
            colorbar()
        
        self.alpha_in = alpha_in
        self.gisans = gisans
        self.gisans_BA = gisans_BA


