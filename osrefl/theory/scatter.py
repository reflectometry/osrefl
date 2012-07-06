# Copyright (C) 2008 University of Maryland
# All rights reserved.
# See LICENSE.txt for details.
# Author: Christopher Metting

#Starting Date:6/19/2009
from numpy import *

from pylab import figure, show, subplot, imshow

#import approximations, view, scale, sample_prep, magPlotSlicer
from . import approximations
import osrefl.viewers.view as view
from osrefl.viewers.plot_slicer import MultiView
import osrefl.loaders.scale
import osrefl.model.sample_prep
from . import resolution
import numpy as np


class Calculator(object):
    '''
    **Overview:**

        This holds all of the information for calculation of scattering for
        reflectometry. This allows a user to build a sample, request an output
        and based on an approximation choice, produce scattering.

    **Parameters:**

        *lattice:* (Lattice)
            see :class:`~sample_prep.Lattice` for more information.

        *probe:* (Beam)
            see :class:`~sample_prep.Beam` for more information.

        *space* (space)
            see :class:`~sample_prep.Q_space` or `~sample_prep.theta_space` for
            more information.

        *feature:* (Unit_Cell)
            see :class:`~sample_prep.Unit_Cell` for more information.

        *omf:* (Omf)
            This is an object which holds the magnetic moment information
            about the sample. It contains three arrays of the the same size
            as the unit cell which hold each of the x, y, and z components
            of the magnetic moment.

    '''
    def __init__(self, lattice, probe, space, feature, omf = None):

        self.feature = feature
        self.omf = omf
        self.lattice = lattice
        self.probe = probe
        self.space = space
        self.results = None
        self.corrected_results = None
        self.anglexvals = None
        self.anglezvals = None
        
        return
    def BAres(self):

        self.results = approximations.BAres(self.feature,self.space,
                                             self.lattice, self.probe)
        print self.results
        return

    def BA(self):
        r'''

        **Overview:**

            This Born Approximation calculation is written entirely in Python
            and assumes that the scattered beam is so small that the transmitted
            beam is essentially t=1. This makes for a simple calculation,
            however, it does not allows for the capturing of the dynamic
            effects seen in real scattering.

            Because of the simplistic nature of this calculation. Some tricks
            can be used to speed up the calculation. This version of the BA
            calculation utilizes the chirp-z transform (CZT) to solve the Form Factor.
            The chirp-z is essentially an FFT which allows for solving the
            transform anywhere on the sphere. With this, we can solve for any Q
            range without wasting any resources calculating for areas we don't
            need.

            The Form Factor calculation is:

            .. math::

                FF = abs \biggl(\frac{-i}{q_{x,y,z}}*(1.0 - exp^{i*q_{x,y,z}*
                \Delta d_{x,y,z}})*CZT\lbrace \rho_{unit}\rbrace\biggr)^{2}

            .. _normalization:

            It is also normalized by the surface area:

            .. math::

               Norm factor = \biggl(\frac{4*\pi}{q_z*M_{x,y}*D_{x,y}}
               \biggr)^{2}


            For the formalism to the structure factor see
            :meth:`~sample_prep.Rectilinear` or :meth:`~sample_prep.Hexagonal`

        '''
        if self.space.type == 'Q_space':

            self.results = approximations.BA(self.feature,self.space,
                                             self.lattice, self.probe)
        elif self.space.type == 'Theta_space':
            self.results = approximations.thetaBA(self.feature,
                                          self.space,self.lattice, self.probe)
        else:
            raise Exception, ("This calculation does not ",
                              self.space.type,"space type")
        return
    
    def BA_FormFactor(self):
        
        raw_intensity = approximations.BA_FT(self.feature.unit, self.feature.step, self.space)
        raw_intensity = abs(approximations.complete_formula(raw_intensity, self.feature.step, self.space))**2
    
        qx_array = self.space.q_list[0].reshape(self.space.points[0],1,1)
        qy_array = self.space.q_list[1].reshape(1,self.space.points[1], 1)
        qz_array = self.space.q_list[2].reshape(1,1,self.space.points[2])
        # * qx_array * qy_array
        raw_intensity *= (4.0 * pi / (qz_array))**2
    
        #raw_intensity = sum(raw_intensity,axis=1).astype('float64')  
    
        return raw_intensity

    def longBA(self):
        r'''
        **Overview:**

            For testing and validation, it can be handy to have a long-hand
            version of the Born Approximation. This BA is written entirely in
            Python and solves the form factor using an explicit sum instead of
            the FFT or CZT modules. It is slow! The form factor is:

            .. math::

                FF = {abs \biggl( \frac{-i}{q_{x,y,z}}*(1.0 -
                exp^{i*q_{x,y,z}* \Delta d_{x,y,z}})*\sum^{D_{x,y,z}}_{n=0}
                \lbrace rho_{unit}*exp^{i*Q_{x,y,z}*D_{x,y,z}}\rbrace
                \biggr)}^{2}

            This form factor is also :ref:`normalized <normalization>`

            For the formalism to the structure factor see
            :meth:`~sample_prep.Rectilinear` or :meth:`~sample_prep.Hexagonal`

        '''
        self.results = approximations.longBA(self.feature,self.space,
                         self.lattice, self.probe)
        return

    def cudaBA(self,precision = 'float32', refract = True):
        r'''
        **Overview:**

            This version of the Born Approximation (BA) uses a scattering kernel
            that is written in C++ and was developed for solving the Substrate
            Modified Born Approximation(SMBA) (see :meth:`cudaSMBA`). This
            kernel normally takes in a set of incoming and outgoing wave
            functions to perturb the probing wave with. Because the BA assumes
            that the wavefunction does not change as a function of sample
            penetration, The incoming and outgoing wavefunctions are set so that
            t = 1, effectively solving the long hand version of the BA (see
            :meth:`longBA`).

            The advantage of this method is that it can distribute the
            calculation across multiple GPU devices solving the problem
            significantly faster.

            This form factor is also :ref:`normalized <normalization>`

            For the formalism to the structure factor see
            :meth:`~sample_prep.Rectilinear` or :meth:`~sample_prep.Hexagonal`


        **Parameters**
            *precision:* (str|precision)

                This parameters allows the user to toggle between float32 and
                float64 precision. For most nvidia graphics cards, the float32
                is handled better and makes for significantly faster
                calculations.

            *refract:* (bool)

                This parameters toggles the refractive shift calculation.
                Generally, the refractive index of the substrate of a sample
                cause a shift in effective Q below the horizons. setting refract
                to TRUE will cause a shift of:

                .. math::

                    q_x+\lambda * \rho_{substrate}

                at :math:`-q_x` values and:

                .. math::

                    q_x-\lambda * \rho_{substrate}

        '''
        self.results = approximations.cudaBA(self.feature,self.space,
                                     self.lattice, self.probe,precision,refract)
        return

    def cudaSMBA(self,precision = 'float32', refract = True):
        r'''
        **Overview**

            The Substrate Modified Born Approximation (SMBA) is a variation of
            the Born Approximation (BA) where by the scattering is perturbed by
            the wavefunction of the income and outgoing wavefunction as it
            interacts with the incident media/substrate interface. This
            perturbation gives rise to the horizons of the sample where the beam
            enters directly from the side face of the substrate.

            This calculation uses C++ calculation kernels on the nvidia GPUs.
            It uses binders from pyCuda to simplify the kernel parallelization.
            There are two C++ calculation kernels used for this calculation.

            The first kernel can be found in :mod:`wavefunction_kernel.cc`. it
            solves for the wavefunction for a wave interacting with the incident
            media/substrate interface. it first calculates the scattering
            vector:

            .. math::

                    k_j = n_jk_0 = \sqrt{1 - \frac{4\pi\rho_j}{k_0}}

            Once this is solved for, the reflection and transmission is
            calculated for the substrate/incident media stack. This is a matrix
            equation:

            .. math::
               \overline{M_l(\Delta z)}=
                \begin{pmatrix}
                cos(k_l\Delta z) & \frac{1}{k_l}sin(k_l\Delta z) \\
                -k_l sin(k_l\Delta z) & cos(k_l\Delta z)
                \end{pmatrix}

            and the reflection is:

            .. math::

                r = \frac{M_{1,1}+(i*n_0*M_{0,1}) +
                \frac{-i}{n_f}*(-M_{1,0}-i*n_0*M_{0,0})}
                {-M_{1,1}+i*n_0*M_{0,1}\frac{-i}{n_f}(M_{1,0}-i*n_0*M_{0,0})}

            and the transmitted beam is:

            .. math::
                t = 1.0 + r

            Now the wavefunction for the perturbation is solved for. The wave
            function used for the perturbation is dependent on the direction of
            the incoming and outgoing beam:

                For :math:`k_{in} < 0.0`:

                    .. math::

                        \Psi_{in_1} = t * \exp^{-ik_{\parallel} \Delta q_z}

                    .. math::

                        \Psi_{in_2} = 0.0

                For :math:`k_{in} > 0.0`:

                    .. math::
                        \Psi_{in_1} = 1 * \exp^{-ik_{\parallel} \Delta q_z}

                    .. math::

                        \Psi_{in_2} = r * \exp^{-ik_{\parallel} \Delta q_z}

                For :math:`k_{out} < 0.0`:

                    .. math::
                        \Psi_{out_1} = 1 * \exp^{-ik_{\parallel} \Delta q_z}

                    .. math::

                        \Psi_{out_2} = r * \exp^{-ik_{\parallel} \Delta q_z}

                For :math:`k_{out} > 0.0`:

                    .. math::

                        \Psi_{out_1} = t * \exp^{-ik_{\parallel} \Delta q_z}

                    .. math::

                        \Psi_{out_2} = 0.0

        With these pieces of information, the SMBA can be solved for. The final
        form factor is:


        .. math::

            FF = {\Biggl( \frac{-i}{q_{x,y,z}}*(1.0 - exp^{i*q_{x,y,z}* \Delta
            d_{x,y,z}})* \biggl[\Psi_{in}* \sum^{D_{x,y,z}}_{n=0} \lbrace
            rho_{unit}*exp^{i*Q_{x,y,z}*D_{x,y,z}}\rbrace *\Psi_{out}\biggr]
            \Biggr)}^{2}

        This form factor is also :ref:`normalized <normalization>`

        For the formalism to the structure factor see
        :meth:`~sample_prep.Rectilinear` or :meth:`~sample_prep.Hexagonal`


        **Parameters**
            *precision:* (str|precision)

                This parameters allows the user to toggle between float32 and
                float64 precision. For most nvidia graphics cards, the float32
                is handled better and makes for significantly faster
                calculations.

            *refract:* (bool)

                This parameters toggles the refractive shift calculation.
                Generally, the refractive index of the substrate of a sample
                cause a shift in effective Q below the horizons. setting refract
                to TRUE will cause a shift of:

                .. math::

                    q_x+\lambda * \rho_{substrate}

                at :math:`-q_x` values and:

                .. math::

                    q_x-\lambda * \rho_{substrate}

        '''
        self.results = approximations.cudaSMBA(self.feature,self.space,
                                           self.lattice, self.probe,precision,
                                           refract)
        return

    def SMBA(self):
        '''
        **Overview:**

            This is a Python implementation of the :meth:`cudaSMBA`. It is
            significantly slower and was only really used for testing and
            validation purposes. Still, it may be useful in the future and is
            available in this package.
        '''
        self.results = approximations.SMBA(self.feature,
                                        self.space,self.lattice, self.probe)
        return

    def SMBAfft(self,precision = 'float32',refract = True):
        '''
        **Overview:**


        '''

        self.results = approximations.SMBAfft(self.feature,
                                        self.space,self.lattice,
                                        self.probe,precision,refract)
        return

    def cudaMagBA(self,precision = 'float32', refract = True):
        self.results = approximations.cudaMagBA(self.feature,self.space,
                                     self.lattice, self.probe,self.omf,
                                     precision,refract)
        return

    def magneticBA(self):
        '''
        **Overview:**

            This calculation solves the Born Approximation for a magnetic sample
            using an :class:`~omfLoader.Omf`.

            For any magnetic scattering, four cross-sections must be solved for.
            First, the magnetic scattering length density must be obtained. This
            can be found


        **Parameters:**

            *struc_cell:* (float:3D array|angstroms^2)
                The structural scattering potential of the feature being
                scattered off of.

            *Q:* (q_space)
                A Q_space object that holds all of the information about the
                desired q space output.

            *lattice:* (Lattice)
                A lattice object that holds all of the information needed to
                solve the structure factor of the scattering.

            *space:* (Space)
                Holds all of the information about the experimental beam needed
                to apply beam dependent corrections to the data.

            *omf:* (Omf)
                This is an object which holds the magnetic moment information
                about the sample. It contains three arrays of the the same size
                as the unit cell which hold each of the x, y, and z components
                of the magnetic moment.

        '''


        self.results = approximations.magneticBA(self.feature,self.space,
                                          self.lattice, self.probe, self.omf)
        return

    def partial_magnetic_BA(self):

        '''
        **Overview:**

            This calculation does the magnetic born approximation but assumes
            that the contribution to the magnetic SLD from the qx and qy
            components of the magnetic moment are negligible and the whole
            system can be estimated as only containing magnetic contribution in
            the qz direction.

            .. warning::
                This method is not accurate for magnetic moments aligned in
                the q directions and should not be used!



        **Parameters:**

            *struc_cell:* (float:3D array|angstroms^2)
                The structural scattering potential of the feature being
                scattered off of.

            *mag_cell:* (float:3D array|angstroms^2)
                The magnetic scattering potential of the feature being scattered
                off of.

            *Q:* (q_space)
                A Q_space object that holds all of the information about the
                desired q space output.

            *lattice:* (Lattice)
                A lattice object that holds all of the information needed to
                solve the structure factor of the scattering.

            *beam:* (Beam)
                Holds all of the information about the experimental beam needed
                to apply beam dependent corrections to the data.

            '''

        self.results = approximations.partial_magnetic_BA(self.feature,
                                                  self.mag_feature, self.space,
                                                  self.lattice, self.probe)

        return

    def partial_magnetic_BA_long(self):
        '''
        '''
        self.results = approximations.partial_magnetic_BA_long(self.feature,
                                                   self.mag_feature, self.space,
                                                   self.lattice, self.probe)
        return

    def DWBA(self,refract = True):
        '''
        **Overview:**


        '''
        from DWBA import DWBA_form
        from numpy import sum
        results = asarray(DWBA_form(self.feature,self.lattice,
                                 self.probe,self.space,refract = refract))
        #print shape(results)
        self.results = sum((abs(results)**2),axis=1)
        return
    
    def DWBA_FormFactor(self,refract = True):
        '''
        **Overview:**
            Returns the formfactor from the DWBA calculation.

        '''
        from DWBA import DWBA_form
        from numpy import sum
        results = asarray(DWBA_form(self.feature,None,
                                 self.probe,self.space,refract = refract))
        
        return (abs(results)**2)
    
    def toAngular(self, incident_angle, intensity):
        
        data = approximations.QxQyQz_to_angle(self.space, incident_angle, intensity, self.probe.wavelength)
        self.results = data[0]
        self.anglexvals = data[1]
        self.anglezvals = data[2]

    def resolution_correction(self):
        '''
        **Overview:**

            Applies a resolution correction to the data using the beam
            information from a :class:`sample_prep.Beam` included in
            :class:`scatter.Calculator` class. It applies a gaussian correction
            for the beam's angular divergence and the divergence in beam energy.

        **Returns**

            *self.corrected_data* (float,array|angstroms^-2)

                Fills in the the values for this attribute of the
                :class:`scatter.Calculator`. If this method is not run, the
                value of this attribute is None.

        '''
        if self.space.minimums[2] < 0.0:

            min = (array([self.space.minimums[0],
                          self.space.minimums[1],self.space.maximums[2] *-1]))

            max = (array([self.space.maximums[0],
                          self.space.maximums[1],self.space.minimums[2] *-1]))

            temp_Q = sample_prep.Q_space(min, max, self.space.points)


            temp_data = fliplr(self.results)

            corrected_results = resolution.Resolution_correction(
                            temp_data, [temp_Q.minimums[0],temp_Q.maximums[0]],
                            [temp_Q.minimums[2],temp_Q.maximums[2]],
                            self.probe.angular_div, self.probe.wavelength_div,
                            self.probe.wavelength)

            self.corrected_results = fliplr(corrected_results.map_out)

        else:
            if len(self.results)==4:
                corrected_results = [None]*4
                self.corrected_results = [None]*4
                for i in range(4):
                    corrected_results[i] = resolution.Resolution_correction(
                        self.results[i], [self.space.minimums[0],
                                       self.space.maximums[0]],
                        [self.space.minimums[2],self.space.maximums[2]],
                        self.probe.angular_div, self.probe.wavelength_div,
                        self.probe.wavelength)
                    self.corrected_results[i] = corrected_results[i].map_out
            else:     
                corrected_results = resolution.Resolution_correction(
                            self.results, [self.space.minimums[0],
                                           self.space.maximums[0]],
                            [self.space.minimums[2],self.space.maximums[2]],
                            self.probe.angular_div, self.probe.wavelength_div,
                            self.probe.wavelength)

                self.corrected_results = corrected_results.map_out

        return


    def view_uncorrected(self, lbl = None, vmin = None, vmax = None):
        '''
        **Overview:**

            This plots the resulting scattering without any resolution
            correction applied to the scattering calculation. This can be useful
            for studying the effects that the resolution has on the data
            measured from the instrument.

        **Parameters:**

            *lbl:* (str)
                This parameter can be used to change the title of the plot.

            *vmin:* (float|Angstroms^-2)
                This is the minimum intensity value plotted on the 2D plot. Any
                intensity value below this value is plotted as the minimum.

            *vmin:* (float|Angstroms^-2)
                This is the maximum intensity value plotted on the 2D plot. Any
                intensity value below this value is plotted as the maximum.

        '''
        if lbl !=None:
            lbl = 'Uncorrected - ', + lbl
        else:
            lbl = 'Uncorrected'
        if size(shape(self.results)) == 2:
            (view.intensity_plot(self.results,self.space.minimums,
                     self.space.maximums,header = lbl,vmin = vmin, vmax = vmax))

        elif size(shape(self.results)) > 2:

            mins = zeros(4)
            maxs = zeros(4)

            for i in range(4):
                mins[i] = (min(log((abs(self.results[i]
                                        [isfinite(self.results[i])]).real))))

                maxs[i] = (max(log((abs(self.results[i]
                                        [isfinite(self.results[i])]).real))))


            for i in range(4):

                subplot(2,2,i+1)
                (view.intensity_plot(self.results[i],
                                     self.space.minimums,self.space.maximums,
                                     'Uncorrected', vmin = vmin, vmax = vmax))

        return

    def viewCorUncor(self):
        '''
        **Overview:**

            Uses the magPlotSlicer.py module to view both the resolution
            corrected and uncorrected models. This module includes tools for:

            * Slice averaging for the data vertically and horizontally
            * Viewing linear and log plots of both 2D slices and 3D image plots
            * Altering of the color axis scale

        '''
        if size(self.results) == 4:
            titles = ['++,uncor','--,uncor','+-,uncor','-+,uncor',
                      '++,cor','--,cor','+-,cor','-+,cor']
            data = [[self.results[0],'Theory'],[self.results[1],'Theory'],
                    [self.results[2],'Theory'],[self.results[3],'Theory'],
                    [self.corrected_results[0],'Theory'],
                    [self.corrected_results[1],'Theory'],
                    [self.corrected_results[2],'Theory'],
                    [self.corrected_results[3],'Theory']]
        else:
            titles = ['uncorrected','corrected']

            data = [[self.results,'Theory1'],[self.corrected_results,'Theory2']]


        extent = self.space.getExtent()
        MultiView(data,[self.space.q_step[0],self.space.q_step[2]],
                  [self.space.points[0],self.space.points[2]],
                  titles=titles,extent= extent,
                  axisLabel = ['qx(A^-1)','qz(A^-1)'])
        return

    def viewUncor(self):
        '''
        **Overview:**

            Uses the magPlotSlicer.py module to view the uncorrected models.
            This module includes tools for:

            * Slice averaging for the data vertically and horizontally
            * Viewing linear and log plots of both 2D slices and 3D image plots
            * Altering of the color axis scale

        '''
        if (shape(self.results)[0]) == 4:
            titles = ['++ (uncor)','-- (uncor)','+- (uncor)','-+ (uncor)']
            data = [[self.results[0],'Theory'],
                    [self.results[1],'Theory'],
                    [self.results[2],'Theory'],
                    [self.results[3],'Theory']]
        else:
            titles = ['uncorrected']
            data = [[self.results,'Theory']]

        extent = self.space.getExtent()
        MultiView(data,[self.space.q_step[0],self.space.q_step[2]],
          [self.space.points[0],self.space.points[2]],
          titles=titles,extent= extent,
          axisLabel = ['qx(A^-1)','qz(A^-1)'])

    def viewCor(self):
        '''
        **Overview:**

            Uses the magPlotSlicer.py module to view the resolution
            corrected models. This module includes tools for:

            * Slice averaging for the data vertically and horizontally
            * Viewing linear and log plots of both 2D slices and 3D image plots
            * Altering of the color axis scale

        '''
        print (len(self.results))
        if (shape(self.results)[0]) == 4:
            print 'there are 4 crosssections'
            titles = ['++ (cor)','-- (cor)','+- (cor)','-+ (cor)']
            data = [[self.corrected_results[0],'Theory'],
                    [self.corrected_results[1],'Theory'],
                    [self.corrected_results[2],'Theory'],
                    [self.corrected_results[3],'Theory']]
        else:
            titles = ['Corrected']
            data = [[self.corrected_results,'Theory']]

        extent = self.space.getExtent()
        MultiView(data,[self.space.q_step[0],self.space.q_step[2]],
          [self.space.points[0],self.space.points[2]],
          titles=titles,extent= extent,
          axisLabel = ['qx(A^-1)','qz(A^-1)'])
        return
    
    def viewAngular(self):
        '''
        **Overview:**

            Uses the magPlotSlicer.py module to view the uncorrected models in real space.
            This module includes tools for:

            * Slice averaging for the data vertically and horizontally
            * Viewing linear and log plots of both 2D slices and 3D image plots
            * Altering of the color axis scale

        '''
        
        x_values = self.anglexvals
        z_values = self.anglezvals
        
        titles = ['uncorrected']
        data = [[self.results,'Theory']]

        xstep = (x_values[size(x_values)-1] - x_values[0]) / size(x_values)
        zstep = (z_values[size(z_values)-1] - z_values[0]) / size(z_values)
        
        xmin = x_values[0]
        xmax = x_values[size(x_values)-1]              
        zmin = z_values[0]
        zmax = z_values[size(z_values)-1]
        
        xmin = xmin.tolist()
        xmax = xmax.tolist()      
        zmin = zmin.tolist()
        zmax = zmax.tolist()
    
        extent = asarray([xmin, xmax, zmin, zmax])
        
        MultiView(data,[xstep,zstep],
          [x_values,z_values],
          titles=titles,extent= extent,
          axisLabel = ['in-plane angle (degrees)','angle of reflection (degrees)'])

    def generalCompare(self,otherData,titles):
        extent = self.space.getExtent()
        data = self.corrected_data + otherData
        MultiView([self.corrected_results,],
          [self.space.q_step[0],self.space.q_step[2]],
          [self.space.points[0],self.space.points[2]],
          titles=titles,extent= extent,
          axisLabel = ['qx(A^-1)','qz(A^-1)'])
        return


    def fitCompare(self,other,extraCompare = None,titles = ['other','self'],cross_section = 'All'):
        '''
        **Overview:**

            This method plots two different data sets on the sample window for
            an easy visual comparison of the data.


        '''
        
        if extraCompare != None:
            data = [[other.data,'Measured']]
            for extraData in [extraCompare]:
                data.append([extraData.data,'Measured'])
            
        else:
            data = [[other.data,'Measured']]

        if len(self.corrected_results) == 4:

            for crossection in self.corrected_results:
                data.append([crossection,'Theory'])

        else:
            data.append([self.corrected_results,'Theory'])
        
        extent = self.space.getExtent()
        print 'data',data
        MultiView(data,[self.space.q_step[0],self.space.q_step[2]],
          [self.space.points[0],self.space.points[2]],
          titles=titles,extent= extent,
          axisLabel = ['qx(A^-1)','qz(A^-1)'])
        return

    def view_linear(self):
        '''
        **Overview:**

            Generally used for testing purposes, this view plots the intensity
            on a linear scale rather then a log scale. This can be useful when
            troubleshooting a calculation.

        '''
        view.linear_plot(self.corrected_results,self.space.minimums,
                         self.space.maximums,'Uncorrected')
        return

    def view_corrected(self,lbl = None, vmin = None, vmax = None):
        '''
        **Overview:**

            This plots the resulting scattering with the resolution correction.
            The user should make sure they have run the
            :meth:`scatter.resolution_correction` method before using this
            method.

        **Parameters:**

            *lbl:* (str)
                This parameter can be used to change the title of the plot.

            *vmin:* (float|Angstroms^-2)
                This is the minimum intensity value plotted on the 2D plot. Any
                intensity value below this value is plotted as the minimum.

            *vmin:* (float|Angstroms^-2)
                This is the maximum intensity value plotted on the 2D plot. Any
                intensity value below this value is plotted as the maximum.

        '''

        view.intensity_plot(self.corrected_results,self.space.minimums,
                            self.space.maximums,'Corrected',
                            vmin = vmin, vmax = vmax)
        return

    def qz_slice(self,qz = 0.0):
        '''
        **Overview:**

            This plots a slice in the y direction designated by the qz. in the
            case where there is corrected resolution data, it also will give the
            qz slice from the resolution corrected data.

        **Parameter:**

            *qz:* (float,angstroms^-1)
                The qz value that is being sliced over. This does not average
                over a range. It will choose the closest qz value that was
                solved for by the calculation and produce a log plot of the
                results.
        '''
        print self.space.minimums
        print self.space.maximums
        view.qz_slice(self.results,self.space.minimums[2],self.space.maximums[2],
                      qz,self.corrected_results)
        return
