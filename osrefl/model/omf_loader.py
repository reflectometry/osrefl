# Copyright (C)
# All rights reserved.
# See LICENSE.txt for details.
# Author: nwright, Brian Maranville
# edited: Christopher Metting
#Starting Date:7/14/2009


from numpy import array, arange, indices, zeros, reshape, float, sqrt, int
from numpy import float64, fromfile, product, ma, isinf, isnan,inf,nan

import wx
DEBUG = True

class Omf(object):
    """
    **Overview:**
    
        When the Object Oriented Micro Magnetic Framework `(OOMMF)
        <http://math.nist.gov/~MDonahue/misc/oommf12a4pre-20091216.tar.gz>`_
        solves the magnetic minimization, it saves the results in a .omf file.
        This class allows the user to load the information from a .omf file
        about the magnetic moments in the sample and save them as a python
        array.
        
        It also works for the oommf12a4pre-20080627 version of OOMMF
    
    **Parameters**
        
        *M* (array,float|angstrom):
            The total magnetic moment vector of the scattering.
        
        *mx* (array,float|angstrom):
            The x component of the magnetic moment vector.
            
        *my* (array,float|angstrom):
            The y component of the magnetic moment vector.
            
        *mz* (array,float|angstrom):
            The z component of the magnetic moment vector.
        
        *parameters* (dictionary,str)
            Holds a dictionary which is generated from the header of the .omf
            file. This is useful for obtaining other information about the model
            run.
            
    .. Note::
        * This class contains other attributes which are not generally used for
          calculation purposes. The user should look in the code for information
          on these attributes.
          
    .. Warning::
        **The .omf file loaded by this module MUST be created from the mmDisplay
        screen.**
        

    **Creating a .omf file for loading**
    
        * Run oommf.tcl 
        
        * Select the appropriate server for processing the
          calculations(this is the local machine for non-distributed
          calculations.)

        * Select oxsii from the mmLaunch box.
        
        * In the Oxsii window select File>>load and select the .mif file created
          by the OsRefl software. (see
          :meth:`sample_prep.Unit_Cell.generateMIF`)
          
        * Run the magnetic minimization by pressing the "Run" button
        
        * Add a mmDisp from the mmLaunch menu
        
            **Selection Input:**
                
            +---------------------+----------------------------+---------------+
            |         Output      |      Destination           |Schedule       |
            +=====================+============================+===============+
            | Magnetization Output|mmDisp<*object for output*> |*Send* Button  |
            +---------------------+----------------------------+---------------+

        * In File >> Save As.. create a .omf file
        
        .. Note::
            * The omf loader supports *Text*, *Binary-4*, and *Binary-8* formats

    """
    def __init__(self, filename = None):
        if filename == None:
            filename = self._menuOpen()

        self.filename = filename
        self.file_object = open(self.filename, 'rb')
        
        self.parameters = {}
        self._parse_header() # this will set self.data_type and all parameters
        
        if self.parameters.has_key('valueunits'):
            self.parameters['valueunit'] = (
                                    self.parameters['valueunits'].split()[0])
            
        if not (self.parameters['meshtype'] == 'rectangular'):
            print 'this only works for rectangular meshes'
            return
        
        self.generate_coordinates() # this sets dimensions and nodes

        if self.data_type == 'Text': self._load_text_data()
        elif self.data_type == 'Binary 8': self._load_binary8_data()
        elif self.data_type == 'Binary 4': self._load_binary4_data()

        else:
            print 'unknown data type'
            return
        
        self.generate_normalized_m()
        self.file_object.close()


    def _menuOpen(self):
        """
        **Overview:**
            
            Used internally to load a .omf from the file system on the computer.
            
        """
        import wx, os
        
        App = wx.PySimpleApp()
        
        dirname = ''
        fileType = "OOMMF output files (.omf)|*.omf"
        dlg = wx.FileDialog(None, "Choose a file", dirname, "", 
                            fileType, wx.OPEN)
        
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
            full_path = os.path.join(dirname, filename)
            
        dlg.Destroy()
        
        return full_path


    def _parse_header_line(self, line):
        '''
        **Overview:**
            
            Used to parse each line of the header section of the .omf file. This
            holds the basic information about the run including the size and
            discretization count of the model.
            
        '''
        key_end = line.find(':')
        # comment_flag = '##'
        if (key_end > -1) and not (line[:2] == '##'):
            key = line[2:key_end]
            parameter = line[key_end + 2:-1]
            return key, parameter
        else:
            return None, None
          

    def _parse_header(self):
        '''
        **Overview:**
            
            Used to parse the header lines of the .omf file. This is used
            to search for extraction of important information from the header
            section.
            
        '''
        self.header = ''
        line = ''
        firstLine = self.file_object.readline()
        
        while not (line[:14] == '# Begin: Data '):
            #keys have a colon at the end, strip it off:
            line = self.file_object.readline()
            self.header += line
            key, param = self._parse_header_line(line)
            if key:
                self.parameters[key] = param
            
        self.data_type = line[14:-1]
        if DEBUG: print self.data_type
        self.parameters['omfType'] = ' '.join(firstLine.split()[2:])
        
    def _load_text_data(self):
        '''
        **Overview:**
            
            Loads the data if the data values are represented in text form.
            
        '''
        xlist = []
        ylist = []
        zlist = []
        #self.file_object.seek(self.data_begin_offset)
        while (1):
            line = self.file_object.readline()
            if line == '# End: Data Text\n':
                break

            xt, yt, zt = line.split()
            xlist.append(float(xt))
            ylist.append(float(yt))
            zlist.append(float(zt))
            
        self.x = array(xlist)
        self.y = array(ylist)
        self.z = array(zlist)

        dims = self.dims
        self.x = reshape(self.x, dims, order='F')
        self.y = reshape(self.y, dims, order='F')
        self.z = reshape(self.z, dims, order='F')


    def _load_binary8_data(self):
        '''
        **Overview:**
            
            These .omf files contain data in two bit formats. For the data that
            is in an 8-bit binary format, this loader can be used to extract
            the data.
            
        '''

        if self.parameters['omfType'] == 'OVF 2.0':
            bitDir = '<f8'
            
        elif self.parameters['omfType'] == 'rectangular mesh v1.0':
            bitDir = '>f8'
            
        else:
            print 'Error: This .omf format is not supported by this software.'
            
        test_entry = fromfile(self.file_object, bitDir, 1)
        
        if DEBUG: print test_entry
        data = fromfile(self.file_object, bitDir, product(self.dims) * 3)
        self.rawdata = data
        expanded_dims = (3, self.dims[0], self.dims[1], self.dims[2])
        self.data = reshape(data, expanded_dims, order='F')
        
        self.x = self.data[0]
        self.y = self.data[1]
        self.z = self.data[2]
        
        data_terminator = self.file_object.readline()
        if DEBUG: print 'terminator = ' + data_terminator


    def _load_binary4_data(self):
        '''
        **Overview:**
            
            These .omf files contain data in two bit formats. For the data that
            is in an 4-bit format, this loader can be used to extract the data.
            
        '''
        if self.parameters['omfType'] == 'OVF 2.0':
            bitDir = '<f4'
            
        elif self.parameters['omfType'] == 'rectangular mesh v1.0':
            bitDir = '>f4'
        else:
            print 'Error: This .omf format is not supported by this software.'
            
        test_entry = fromfile(self.file_object, bitDir, 1)
        if DEBUG: print test_entry
        data = fromfile(self.file_object, bitDir, product(self.dims) * 3)
        expanded_dims = (3, self.dims[0], self.dims[1], self.dims[2])
        self.data = reshape(data, expanded_dims, order='F')
        self.x = self.data[0]
        self.y = self.data[1]
        self.z = self.data[2]
        
        data_terminator = self.file_object.readline()
        if DEBUG: print 'terminator = ' + data_terminator
  
  
    def generate_coordinates(self):
        """
        **Overview:**
            
            Calculates  the x, y  and z values for each of the discretized units
            in the model from the information obtained from the header file.
            
        """
        dims = (int(self.parameters['xnodes']), int(self.parameters['ynodes']), 
                int(self.parameters['znodes']))
        
        self.dims = dims

        if self.parameters['meshtype'] == 'rectangular': 
            self.node_x, self.node_y, self.node_z = indices(
                                                    self.dims, dtype = float)
            
            self.node_x *= float(self.parameters['xstepsize'])
            self.node_x += float(self.parameters['xmin'])
            self.node_y *= float(self.parameters['ystepsize'])
            self.node_y += float(self.parameters['ymin'])
            self.node_z *= float(self.parameters['zstepsize'])
            self.node_z += float(self.parameters['zmin'])
        
        return

    def generate_normalized_m(self):
        '''
        **Overview:**
            
            The moments given in this file are the absolute magnitudes. This 
            method normalizes the data by the total moment.
            
        '''
        
        
        M = sqrt((self.x)**2 + (self.y)**2 + (self.z)**2)
        mask = M.nonzero()
        mx = self.x.copy()
        my = self.y.copy()
        mz = self.z.copy()

        mx[mask] = mx[mask]/M[mask]
        my[mask] = my[mask]/M[mask]
        mz[mask] = mz[mask]/M[mask]
  
        self.M = M
        self.mx = mx
        self.my = my
        self.mz = mz
        
    def downsample(self, down_factor=10):
        """
        **Overview:**
        
            This method resamples x,y data into bigger boxes. It does this by
            averaging the surrounding moments and assigning this weighted
            average to the rest of the data.
        
        .. Note::
            * Qx, Qy resolution is typically much worse than exchange length Qz
              resolution is pretty good, so this method does not resample in the
              z direction.
            
        .. Warning::
            **This module requires file *rebin_simple.py* which is not a common
            package.**
            
        """
        
        new_dims = (self.dims[0]/down_factor, self.dims[1]/down_factor, 
                                                                self.dims[2])
        
        from rebin_simple import congrid
        new_mx = congrid(self.mx, new_dims)
        new_my = congrid(self.my, new_dims)
        new_mz = congrid(self.mz, new_dims)
        new_node_x = congrid(self.node_x, new_dims)
        new_node_y = congrid(self.node_y, new_dims)
        new_node_z = congrid(self.node_z, new_dims)
        new_M = congrid(self.M, new_dims)
        
        self.mx = new_mx
        self.my = new_my
        self.mz = new_mz
        self.node_x = new_node_x
        self.node_y = new_node_y
        self.node_z = new_node_z
        self.M = new_M
        
    def ConvertRho(self):
        '''
        **Overview:**
            
            There is a factor :math:`C'` which is used to convert the magnetic
            moment to a magnetic scattering length density. Because the OOMMF
            software allows for different units, the :math:`C` must be chosen
            based on the OOMMF model.
            
        **Returns** 
            (array[3]|angstroms^-2)
            
        '''
        CPRIME_Am = 2.853e-12
        CPRIME_mT = 2.3164e-9
        
        if (self.parameters['valueunit'] == 'A/m'):
            C_prime = CPRIME_Am  
        elif (self.parameters['valueunit'] == 'mT'):
            C_prime = CPRIME_mT
        else:
            print "units not known \n"

        return self.M * C_prime 

    
    def viewFixedZ(self, plot_title = None, z_layer = 0):
        """
        **Overview:**
        
            This method shows a color plot of the angle between mx, my.
            
        """
        from numpy import arctan2, ma, pi
        from matplotlib.colors import LinearSegmentedColormap
        from pylab import figure,show
        
        wheel_cmap = LinearSegmentedColormap.from_list('wheel_rgby', 
                                   ['red', 'green', 'blue', 'yellow', 'red'])
        
        x,y = indices((100,100), dtype = float) - 50.
        wheel_angle = arctan2(x,y)
        
        fig1 = figure(frameon = False)
        ax1 = fig1.add_axes([0.1,0.1,0.7,0.8])
        
        angle = ma.array(arctan2(-self.my[:,:,z_layer], self.mx[:,:,z_layer]), 
                         mask=(self.M[:,:,z_layer] == 0))
        
        xmax = (float(self.parameters['xnodes']) * 
                float(self.parameters['xstepsize']) * 1.0e10)
        
        ymax = (float(self.parameters['ynodes']) * 
                float(self.parameters['ystepsize']) * 1.0e10)
        
        extent = [0., xmax, 0., ymax]
        
        im = ax1.imshow(angle.T, origin='lower', interpolation = 'nearest', 
                                            extent = extent, cmap = wheel_cmap, 
                                                        vmin = -pi, vmax = pi)
        
        ax1.set_title(plot_title)
        ax1.set_xlabel('$x (\AA)$', size = 'x-large')
        ax1.set_ylabel('$y (\AA)$', size = 'x-large')

        ax1w = fig1.add_axes([0.75,0.4,0.2,0.2], polar=True)
        ax1w.yaxis.set_visible(False)
        ax1w.imshow(wheel_angle, cmap=wheel_cmap, extent=[0,2*pi,0,pi])
        ax1w.set_title('M direction\n(in-plane)')
        
        show()

    def view(self):
        from wxzslice import momentView
        momentView(self)
            
def _test():
    from numpy import shape
    from sample_prep import *
    from scatter import *
    from pylab import *
    
    #mag = Omf('/home/mettingc/Documents/test.omf')
    #mag.view()
    mag = Omf('/home/mettingc/Documents/temp_mif-Oxs_MinDriver-Magnetization-01-0001400.omf')
    #mag.view()

    NiFe = (sample_prep.Parallelapiped(SLD = 9.02e-6,
                                     Ms = 8.6e5,dim=[5.0e4,5.0e4,2.0e4]))
    
    scene = sample_prep.Scene([NiFe])
    
    GeoUnit = (sample_prep.GeomUnit(Dxyz = [10.0e4,10.0e4,2.5e4],
                                    n = [70,70,30],scene = scene))
    
    unit = GeoUnit.buildUnit()
    #unit.generateMIF()
    #unit.viewSlice()
    space = Q_space([-.0001,-0.001,0.00002],[.0001,0.001,0.002],[500,20,500])
    lattice = Rectilinear([20,20,1],unit)
    
    beam = Beam(5.0,None,None,0.05,None)
    

    magtest = Calculator(lattice, beam, space, unit, mag)
    #magtest.magneticBA()
    magtest.cudaMagBA()
    print shape(magtest.results)
    #print magtest.results[isnan(magtest.results)]
    #magtest.resolution_correction()

    magtest.viewUncor()
    
if __name__=="__main__":_test()



    