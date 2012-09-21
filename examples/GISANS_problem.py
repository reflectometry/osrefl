#from greens_thm_form import greens_form_line, greens_form_shape
from greens_thm_form import div_form_shape as greens_form_shape
from numpy import arange, linspace, float64, indices, zeros_like, ones_like, pi, sin, complex128, array, exp, newaxis, cumsum, sum, cos, sin, log, log10, zeros, sqrt, ones
from osrefl.theory.DWBAGISANS import dwbaWavefunction
from gaussian_envelope import FWHM_to_sigma, normgauss

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
                 Lx,Ly,
                 autoFT=True,
                 name='grazing_incidence'):
        self.name = name
        self.sublayers = sublayers
        self.matrix = matrix
        self.Lx = Lx
        self.Ly = Ly
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
    
    def update_Qs(self, alpha_in=None):
        if alpha_in is not None:
            self.alpha_in = alpha_in
        k0 = 2*pi/self.wavelength
        kz_in = array([[[k0 * sin(self.alpha_in * pi/180.0)]]], dtype=complex128)
        kx_in = array([k0 * cos(self.alpha_in * pi/180.0)], dtype=complex128)
        
        kz_out = kz_in - self.qz
        ky_out = -self.qy
        kx_out = sqrt(k0**2 - kz_out**2 - ky_out**2)
        
        self.kz_in = kz_in
        self.kz_out = kz_out
        self.qx = kx_in - kx_out
        
    def update_FTs(self):
        dFTs = [] # differential = SLD - (avg. SLD)
        FTs = []
        for sl in self.sublayers:
            dFT = zeros((self.qx.shape[0],self.qy.shape[1]), dtype=complex128)
            #FT = zeros((self.qx.shape[0], self.qy.shape[0]), dtype=complex128)
            qx = self.qx[:,:,0]
            qy = self.qy[:,:,0]
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
    
    def calc_overlap(self):
        wf_in = dwbaWavefunction(self.kz_in, self.SLDArray)
        wf_out = dwbaWavefunction(-self.kz_out, self.SLDArray) # solve 1d equation for time-reversed state
        self.wf_in = wf_in
        self.wf_out = wf_out
        
        kz_in_l = wf_in.kz_l # inside the layers
        kz_out_l = -wf_out.kz_l # inside the layers
        
        dz = self.SLDArray[1:-1,1][:,newaxis,newaxis,newaxis]
        zs = cumsum(self.SLDArray[1:-1,1]) - self.SLDArray[1,1] # start at zero with first layer
        z_array = array(zs)[:,newaxis,newaxis,newaxis]
        thickness = sum(self.SLDArray[1:-1,1])
        
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
        return overlap
    
    def calc_overlap_BA(self):
        dz = self.SLDArray[1:-1,1][:,newaxis,newaxis,newaxis]
        zs = cumsum(self.SLDArray[1:-1,1]) - self.SLDArray[1,1] # start at zero with first layer
        z_array = array(zs)[:,newaxis,newaxis,newaxis]
        
        overlap_BA  = 1.0 / (1j * self.qz) * (exp(1j * self.qz * dz) - 1.0) * exp(1j*self.qz*z_array)
        self.overlap_BA = overlap_BA
        return overlap_BA
    
    def calc_gisans(self, alpha_in=None, show_plot=True, add_specular=False):
        if alpha_in is not None:
            self.update_Qs(alpha_in)
        overlap = self.calc_overlap()
        gisans = sum(sum(overlap * array(self.dFTs)[:,:,:,newaxis], axis=0), axis=0) # first over layers, then Qx
        # now if you want to add specular back in...
        if add_specular == True:
            specular = ones((self.qx.shape[0], self.qy.shape[1], self.qz.shape[2]), dtype=complex128)
            specular *= complex128(2)*pi/self.Lx * normgauss(self.qx, FWHM_to_sigma(2.0*pi/self.Lx), x0=0.0)
            specular *= complex128(2)*pi/self.Ly * normgauss(self.qy, FWHM_to_sigma(2.0*pi/self.Ly), x0=0.0)
            specular *= 2.0*1j*self.kz_in*self.wf_in.r*self.Lx*self.Ly        
            specular = sum(specular, axis=0)/self.qx.shape[0] # sum over Qx, taking average
            self.specular = specular
            gisans += specular
        self.gisans = gisans
        if show_plot == True:
            self.plot_gisans()            
    
    def calc_gisans_BA(self, show_plot=True):
        overlap_BA = self.calc_overlap_BA()   
        gisans_BA = sum(sum(overlap_BA * array(self.FTs)[:,:,:,newaxis], axis=0), axis=0)
        self.gisans_BA = gisans_BA 
        if show_plot == True: 
            self.plot_gisans_BA()
        
    def calc_both(self, show_plot=True, add_specular=False):
        self.calc_gisans(show_plot=False, add_specular=add_specular)
        self.calc_gisans_BA(show_plot=False)
        if show_plot == True: self.plot_both()
        
    def plot_gisans(self, vmax=None, vmin=None):    
        from pylab import imshow, figure, colorbar
        extent = [self.qy.min(), self.qy.max(), self.qz.min(), self.qz.max()]
        figure()
        imshow(log10(abs(self.gisans)**2).T, origin='lower', extent=extent, aspect='auto', vmax=vmax, vmin=vmin)
        colorbar()
        
    def plot_gisans_BA(self, vmax=None, vmin=None):
        from pylab import imshow, figure, colorbar
        extent = [self.qy.min(), self.qy.max(), self.qz.min(), self.qz.max()]
        figure()
        imshow(log10(abs(self.gisans_BA)**2).T, origin='lower', extent=extent, aspect='auto', vmax=vmax, vmin=vmin)
        colorbar()
        
    def plot_both(self):
        vmax = max(log10(abs(self.gisans)**2).max(), log10(abs(self.gisans_BA)**2).max())
        vmin = min(log10(abs(self.gisans)**2).min(), log10(abs(self.gisans_BA)**2).min())
        self.plot_gisans(vmax=vmax, vmin=vmin)
        self.plot_gisans_BA(vmax=vmax, vmin=vmin)

class GISANS_angle_problem(GISANS_problem):
    def __init__(self, 
                 sublayers,
                 matrix,
                 front_sld, front_sldi, 
                 back_sld, back_sldi, 
                 wavelength, 
                 angle_in, angle_out, inplane_angle,
                 Lx,Ly,
                 autoFT=True,
                 name='grazing_incidence'):
        self.name = name
        self._qx = self._qy = self._qz = None
        self.sublayers = sublayers
        self.matrix = matrix
        self.Lx = Lx
        self.Ly = Ly
        self.front_sld = front_sld
        self.front_sldi = front_sldi
        self.back_sld = back_sld
        self.back_sldi = back_sldi
        self.wavelength = wavelength
        self._angle_in = angle_in
        self._angle_out = angle_out
        self._inplane_angle = inplane_angle
        self.update_SLDArray()
        self.alpha_in = None
        self.FTs = []
        if autoFT == True: 
            self.update_Qs()    
            self.update_FTs()
    
    def get_angle_in(self):
        return self._angle_in
    def set_angle_in(self, value):
        self._angle_in = value
        self.update_Qs()
        self.update_FTs()
    def del_angle_in(self):
        del self._angle_in
    angle_in = property(get_angle_in, set_angle_in, del_angle_in, "I'm the angle_in property.")
    
    def get_angle_out(self):
        return self._angle_out
    def set_angle_out(self, value):
        self._angle_out = value
        self.update_Qs()
        self.update_FTs()
    def del_angle_out(self):
        del self._angle_out
    angle_out = property(get_angle_out, set_angle_out, del_angle_out, "I'm the angle_out property.")
    
    def get_inplane_angle(self):
        return self._inplane_angle
    def set_inplane_angle(self, value):
        self._inplane_angle = value
        self.update_Qs()
        self.update_FTs()
    def del_inplane_angle(self):
        del self._inplane_angle
    inplane_angle = property(get_inplane_angle, set_inplane_angle, del_inplane_angle, "I'm the inplane_angle property.")
        
    def update_Qs(self):
        wavelength = self.wavelength

        # convert angle to radians
        angle_in = self.angle_in * pi / 180.
        angle_out = self.angle_out * pi/180.
        iptheta = self.inplane_angle * pi/180.
        
        # determine wave vector (k)
        kvec = 2.0*pi/wavelength

        kz_out = kvec * sin( -angle_out )[newaxis,newaxis,:]
        kx_out = kvec * cos( -angle_out )[newaxis,newaxis,:]
        ky_out = -kvec * cos( -angle_out ) * sin( iptheta )[newaxis,:,newaxis]
        
        kz_in = kvec * sin( angle_in )
        kx_in = kvec * cos( angle_in )
        ky_in = zeros_like( ky_out )
        
        self._qx = kx_in - kx_out
        self._qy = ky_in - ky_out
        self._qz = kz_in - kz_out
        self.kz_in = kz_in
        self.kz_out = kz_out
        self.update_FTs()
        
    def plot_gisans(self, vmax=None, vmin=None):    
        from pylab import imshow, figure, colorbar, xlabel, ylabel, title
        extent = [self.inplane_angle.min(), self.inplane_angle.max(), self.angle_out.min(), self.angle_out.max()]
        figure()
        imshow(log10(abs(self.gisans)**2).T, origin='lower', extent=extent, aspect='auto', vmax=vmax, vmin=vmin)
        title('%s GISANS, angle_in = %g degrees' % (self.name, self.angle_in))
        colorbar()
        
    def plot_gisans_BA(self, vmax=None, vmin=None):
        from pylab import imshow, figure, colorbar, xlabel, ylabel, title
        extent = [self.inplane_angle.min(), self.inplane_angle.max(), self.angle_out.min(), self.angle_out.max()]
        figure()
        imshow(log10(abs(self.gisans_BA)**2).T, origin='lower', extent=extent, aspect='auto', vmax=vmax, vmin=vmin)
        title('%s GISANS (Born Approximation), angle_in = %g degrees' % (self.name, self.angle_in))
        colorbar()
        
class GISANS_problem_old(object):
    def __init__(self, 
                 sublayers,
                 matrix,
                 front_sld, front_sldi, 
                 back_sld, back_sldi, 
                 wavelength, 
                 qx, qy, qz,
                 Lx,Ly,
                 autoFT=True):
        self.sublayers = sublayers
        self.matrix = matrix
        self.Lx = Lx
        self.Ly = Ly
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
    
    def update_Qs(self, alpha_in=None):
        if alpha_in is not None:
            self.alpha_in = alpha_in
        k0 = 2*pi/self.wavelength
        kz_in = array([k0 * sin(self.alpha_in * pi/180.0)], dtype=complex128)
        kx_in = array([k0 * cos(self.alpha_in * pi/180.0)], dtype=complex128)
        
        kz_out = kz_in - self.qz
        ky_out = -self.qy
        kx_out = sqrt(k0**2 - kz_out[newaxis,newaxis,:]**2 - ky_out[newaxis,:,newaxis]**2)
        
        self.kz_in = kz_in
        self.kz_out = kz_out
        self.qx = kx_in - kx_out
        
    def update_FTs(self):
        dFTs = [] # differential = SLD - (avg. SLD)
        FTs = []
        for sl in self.sublayers:
            dFT = zeros((self.qx.shape[0], self.qy.shape[1]), dtype=complex128)
            #FT = zeros((self.qx.shape[0], self.qy.shape[0]), dtype=complex128)
            qx = self.qx[:,:,0]
            qy = self.qy[:,:,0]
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
        
    def calc_gisans(self, alpha_in, show_plot=True, add_specular=False):
        k0 = 2*pi/self.wavelength
        kz_in_0 = array([k0 * sin(alpha_in * pi/180.0)], dtype=complex128)
        kx_in_0 = array([k0 * cos(alpha_in * pi/180.0)], dtype=complex128)
        kz_out_0 = kz_in_0 - self.qz
        self.kz_out_0 = kz_out_0
        
        ky_out_0 = -self.qy
        kx_out_0 = sqrt(k0**2 - kz_out_0[newaxis,newaxis,:]**2 - ky_out_0[newaxis,:,newaxis]**2)
        qx = kx_in_0 - kx_out_0
        self.qx_derived = qx

        kz_out_neg = kz_out_0 < 0
        kz_in_neg = kz_in_0 < 0
        
        wf_in = dwbaWavefunction((kz_in_0), self.SLDArray)
        wf_out = dwbaWavefunction((-kz_out_0), self.SLDArray) # solve 1d equation for time-reversed state
        self.wf_in = wf_in
        self.wf_out = wf_out
        
        kz_in_l = wf_in.kz_l # inside the layers
        #kz_in_l[:, kz_in_neg] *= -1.0     
        kz_in_p_l = -kz_in_l # prime
        kz_out_l = -wf_out.kz_l # inside the layers
        #kz_out_l[:, kz_out_neg] *= -1.0
        kz_out_p_l = -kz_out_l    # kz_f_prime in the Sinha paper notation
        
        dz = self.SLDArray[1:-1,1][:,newaxis]
        zs = cumsum(self.SLDArray[1:-1,1]) - self.SLDArray[1,1] # start at zero with first layer
        z_array = array(zs)[:,newaxis]
        thickness = sum(self.SLDArray[1:-1,1])

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
        gisans = sum(sum(overlap * array(self.dFTs)[:,:,:,newaxis], axis=0), axis=0) # first over layers, then Qx
        # now if you want to add specular back in...
        if add_specular == True:
            specular = complex128(2)*pi/self.Lx * normgauss(qx, FWHM_to_sigma(2.0*pi/self.Lx), x0=0.0)
            specular *= complex128(2)*pi/self.Ly * normgauss(self.qy[newaxis,:,newaxis], FWHM_to_sigma(2.0*pi/self.Ly), x0=0.0)
            specular *= 2.0*1j*kz_in_0*wf_in.r[newaxis,newaxis,:]*self.Lx*self.Ly        
            specular = sum(specular, axis=0)/self.qx.shape[0] # sum over Qx, taking average
            self.specular = specular
            gisans += specular
        
        gisans_BA = sum(sum(overlap_BA * array(self.FTs)[:,:,:,newaxis], axis=0), axis=0) 
        extent = [self.qy.min(), self.qy.max(), self.qz.min(), self.qz.max()]
        
        self.alpha_in = alpha_in
        self.gisans = gisans
        self.gisans_BA = gisans_BA
        
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
