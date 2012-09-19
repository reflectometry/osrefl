#from greens_thm_form import greens_form_line, greens_form_shape
from greens_thm_form import div_form_shape as greens_form_shape
from numpy import arange, linspace, float64, indices, zeros_like, ones_like, pi, sin, complex128, array, exp, newaxis, cumsum, sum, cos, sin, log, log10, zeros
from osrefl.theory.DWBAGISANS import dwbaWavefunction
from osrefl.theory.approximations import QxQyQz_to_k
from GISANS_problem import Shape, GISANS_problem

class OFFSPEC_problem(GISANS_problem):
    
    def calc_gisans(self, alpha_in, show_plot=True):
        print "this class is for offspec, not gisans"
        return
        
    def calc_offspec(self, show_plot=True):
        kz_in_0, kz_out_0 = QxQyQz_to_k(self.qx[:,newaxis,newaxis],self.qy[newaxis,:,newaxis],self.qz[newaxis,newaxis,:],self.wavelength)
        kzi_shape = kz_in_0.shape
        kzf_shape = kz_out_0.shape
        
        kz_out_neg = kz_out_0 < 0
        kz_in_neg = kz_in_0 < 0
        
        wf_in = dwbaWavefunction(abs(kz_in_0), self.SLDArray)
        wf_out = dwbaWavefunction(abs(kz_out_0), self.SLDArray) # solve 1d equation for time-reversed state
        
        kz_in_l = wf_in.kz_l # inside the layers
        kz_in_l[:, kz_in_neg] *= -1.0     
        kz_in_p_l = -kz_in_l # prime
        kz_out_l = wf_out.kz_l # inside the layers
        kz_out_l[:, kz_out_neg] *= -1.0
        kz_out_p_l = -kz_out_l    # kz_f_prime in the Sinha paper notation

        dz = self.SLDArray[1:-1,1][:,newaxis,newaxis,newaxis]
        zs = cumsum(self.SLDArray[1:-1,1]) - self.SLDArray[1,1] # start at zero with first layer
        z_array = array(zs)[:,newaxis,newaxis,newaxis]

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
        offspec = sum(sum(overlap * array(self.FTs)[:,:,:,newaxis], axis=0), axis=1) # first over layers, then Qy
        offspec_BA = sum(sum(overlap_BA * array(self.FTs)[:,:,:,newaxis], axis=0), axis=1) 
        extent = [self.qx.min(), self.qx.max(), self.qz.min(), self.qz.max()]
        self.offspec = offspec
        self.offspec_BA = offspec_BA
        
        if show_plot == True:
            from pylab import imshow, figure, colorbar
            zmax = max(log10(abs(offspec)**2).max(), log10(abs(offspec_BA)**2).max())
            zmin = min(log10(abs(offspec)**2).min(), log10(abs(offspec_BA)**2).min())
            figure()
            imshow(log10(abs(offspec)**2).T, origin='lower', extent=extent, aspect='auto', vmax=zmax, vmin=zmin)
            colorbar()
            
            figure()
            imshow(log10(abs(offspec_BA)**2).T, origin='lower', extent=extent, aspect='auto', vmax=zmax, vmin=zmin)
            colorbar()
        
        


