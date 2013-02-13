#from greens_thm_form import greens_form_line, greens_form_shape
from greens_thm_form import div_form_shape as greens_form_shape
from numpy import arange, linspace, float64, indices, zeros_like, ones_like, pi, sin, complex128, array, exp, newaxis, cumsum, sum, cos, sin, log, log10, zeros, sqrt, ones
from osrefl.theory.DWBAGISANS import dwbaWavefunction
from osrefl.theory.approximations import QxQyQz_to_k
from GISANS_problem import Shape, GISANS_problem
from gaussian_envelope import FWHM_to_sigma, normgauss
from pylab import find

class OFFSPEC_problem(GISANS_problem):
    
    def calc_gisans(self, alpha_in, show_plot=True, **kwargs):
        print "this class is for offspec, not gisans"
        return
        
    def calc_both(self, show_plot=True, add_specular=False):
        self.calc_offspec(show_plot=False, add_specular=add_specular)
        self.calc_offspec_BA(show_plot=False, add_specular=add_specular)
        if show_plot == True: self.plot_both()
                
    def plot_both(self):
        vmax = max(log10(abs(self.offspec)**2).max(), log10(abs(self.offspec_BA)**2).max())
        vmin = min(log10(abs(self.offspec)**2).min(), log10(abs(self.offspec_BA)**2).min())
        self.plot_offspec(vmax=vmax, vmin=vmin)
        self.plot_offspec_BA(vmax=vmax, vmin=vmin)
        
    def update_Qs(self):
        kz_in_0, kz_out_0 = QxQyQz_to_k(self.qx,self.qy,self.qz,self.wavelength)
        self.kz_in = kz_in_0
        self.kz_out = kz_out_0
        
    def calc_offspec(self, show_plot=True, add_specular=True):
        overlap = self.calc_overlap()
        offspec = sum(sum(overlap * array(self.dFTs)[:,:,:,newaxis], axis=0), axis=1) # first over layers, then Qy
        
        if add_specular == True:
            specular = ones((self.qx.shape[0], self.qy.shape[1], self.qz.shape[2]), dtype=complex128)
            specular *= complex128(2)*pi/self.Lx * normgauss(self.qx, FWHM_to_sigma(2.0*pi/self.Lx), x0=0.0)
            specular *= complex128(2)*pi/self.Ly * normgauss(self.qy, FWHM_to_sigma(2.0*pi/self.Ly), x0=0.0)
            specular *= 2.0*1j*self.kz_in*self.wf_in.r*self.Lx*self.Ly        
            specular = sum(specular, axis=1)/self.qy.shape[1] # sum over Qy, taking average
            self.specular = specular
            offspec += specular
        self.offspec = offspec
        
        if show_plot == True:
            self.plot_offspec()
    
    def calc_offspec_BA(self, show_plot=True, add_specular=True):
        overlap = self.calc_overlap_BA()
        offspec = sum(sum(overlap * array(self.dFTs)[:,:,:,newaxis], axis=0), axis=1) # first over layers, then Qy
        
        if add_specular == True:
            specular = ones((self.qx.shape[0], self.qy.shape[1], self.qz.shape[2]), dtype=complex128)
            specular *= complex128(2)*pi/self.Lx * normgauss(self.qx, FWHM_to_sigma(2.0*pi/self.Lx), x0=0.0)
            specular *= complex128(2)*pi/self.Ly * normgauss(self.qy, FWHM_to_sigma(2.0*pi/self.Ly), x0=0.0)
            specular *= 2.0*1j*self.kz_in*self.wf_in.r*self.Lx*self.Ly        
            specular = sum(specular, axis=1)/self.qy.shape[1] # sum over Qy, taking average
            self.specular = specular
            offspec += specular
        self.offspec_BA = offspec
        
        if show_plot == True:
            self.plot_offspec_BA()
    
    def plot_offspec(self, vmax=None, vmin=None):    
        from pylab import imshow, figure, colorbar, xlabel, ylabel, title
        extent = [self.qx.min(), self.qx.max(), self.qz.min(), self.qz.max()]
        figure()
        imshow(log10(abs(self.offspec)**2).T, origin='lower', extent=extent, aspect='auto', vmax=vmax, vmin=vmin)
        title('%s offspecular scattering' % (self.name,))
        colorbar()
        
    def plot_offspec_BA(self, vmax=None, vmin=None):    
        from pylab import imshow, figure, colorbar, xlabel, ylabel, title
        extent = [self.qx.min(), self.qx.max(), self.qz.min(), self.qz.max()]
        figure()
        imshow(log10(abs(self.offspec_BA)**2).T, origin='lower', extent=extent, aspect='auto', vmax=vmax, vmin=vmin)
        title('%s offspecular scattering (Born)' % (self.name,))
        colorbar()
        
    def calc_offspec_old(self, show_plot=True, add_specular=True):
        kz_in_0, kz_out_0 = QxQyQz_to_k(self.qx[:,newaxis,newaxis],self.qy[newaxis,:,newaxis],self.qz[newaxis,newaxis,:],self.wavelength)
        kzi_shape = kz_in_0.shape
        kzf_shape = kz_out_0.shape
        
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

        dz = self.SLDArray[1:-1,1][:,newaxis,newaxis,newaxis]
        zs = cumsum(self.SLDArray[1:-1,1]) - self.SLDArray[1,1] # start at zero with first layer
        z_array = array(zs)[:,newaxis,newaxis,newaxis]
        thickness = sum(self.SLDArray[1:-1,1])

        qrt_inside = -kz_in_l[1:-1] - kz_out_l[1:-1]
        qtt_inside = -kz_in_l[1:-1] + kz_out_l[1:-1]
        qtr_inside = +kz_in_l[1:-1] + kz_out_l[1:-1]
        qrr_inside = +kz_in_l[1:-1] - kz_out_l[1:-1]
        
        self.qrt = qrt_inside
        self.qtt = qtt_inside
        self.qtr = qtr_inside
        self.qrr = qrr_inside
               
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
        offspec = sum(sum(overlap * array(self.dFTs)[:,:,:,newaxis], axis=0), axis=1) # first over layers, then Qy
        # now need to add specular back in...
        if add_specular == True:
            specular = 2.0*1j*kz_in_0*wf_in.r*self.Lx*self.Ly
        
            #specular *= 2*pi/self.Lx * normgauss(self.qx[:,newaxis,newaxis], FWHM_to_sigma(2.0*pi/self.Lx), x0=0.0)
            #specular *= 2*pi/self.Ly * normgauss(self.qy[newaxis,:,newaxis], FWHM_to_sigma(2.0*pi/self.Ly), x0=0.0)
            specular = sum(specular, axis=1)/kz_in_0.shape[1] # sum over Qy, taking average
            self.specular = specular
            qx_min = find(abs(self.qx) == abs(self.qx).min())[0]
            offspec[qx_min] += specular[qx_min]
        
            self.R_tilde_o_sq = abs(offspec)**2 / (4.0 * self.Lx**2 * self.Ly**2 * kz_in_0[:,0,:]**2)
        
        offspec *= 1.0/(self.Lx * self.Ly)

        offspec_BA = sum(sum(overlap_BA * array(self.FTs)[:,:,:,newaxis], axis=0), axis=1)
        offspec_BA *= 1.0/(self.Lx * self.Ly)
        extent = [self.qx.min(), self.qx.max(), self.qz.min(), self.qz.max()]
        self.offspec = offspec
        self.offspec_BA = offspec_BA
        
        if show_plot == True:
            from pylab import imshow, figure, colorbar
            zmax = max(log10(abs(offspec)**2).max(), log10(abs(offspec_BA)**2).max())
            zmin = min(log10(abs(offspec)**2).min(), log10(abs(offspec_BA)**2).min())
            figure()
            imshow(log10(self.R_tilde_o_sq.real).T, origin='lower', extent=extent, aspect='auto')
            colorbar()
            
            figure()
            imshow(log10(abs(offspec)**2).T, origin='lower', extent=extent, aspect='auto', vmax=zmax, vmin=zmin)
            colorbar()
            
            figure()
            imshow(log10(abs(offspec_BA)**2).T, origin='lower', extent=extent, aspect='auto', vmax=zmax, vmin=zmin)
            colorbar()
        
        


