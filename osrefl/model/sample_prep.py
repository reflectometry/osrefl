# Copyright (C) 2008 University of Maryland
# All rights reserved.
# See LICENSE.txt for details.
# Author: Christopher Metting

#Starting Date:3/12/2009

import copy, re, string

from numpy import array, asarray, arange, linspace, append, empty, zeros, ones
from numpy import rot90, flipud,fliplr, repeat, nonzero,radians
from numpy import float, int, round, float32
from numpy import shape, size, hstack, concatenate, reshape, abs,repeat
from numpy import abs, pi, min, max, exp, sum, log, cos,sin,arccos,arctan,arcsin
from numpy import round, average,sqrt
from scipy.integrate import quad
from . import calculations, omf_loader, image_util
from ..viewers import view


class Unit_Cell(object):
    '''
    **Overview:**

        Contains the information for the processed unit cell information. This
        class makes the implementation of scattering calculations easier by by
        creating one structure that only contains information necessary for the
        theory function calculation, regardless of how that information was
        obtained.


    **Producer Classes:**


    *GeomUnit*

    Uses geometric shape parameters to calculate the discretized unit cell.
    Magnetic support is included in the form of a unit cell which contains the
    magnetic SLD values.


    *OOMMFUnit*

    Uses the .omf output format produced by the Object Oriented MicroMagnetic
    Framework software. to created both
    the structure and magnetic scattering arrays. The Unit_Cell objected created
    by this class will contain the structural SLD as well as a list of arrays
    which contain the x,y, and z magnetic components. This can be used by the
    magnetic approximations to calculated the four magnetic cross- sections.


    *K3DUnit*

    Loads a raw data fill which is exported by the k3d modeling software `K-3D
    software <http://www.k-3d.org/>`_. This type of creation does not support
    magnetic representations.


    *GrayImgUnit*

    Loads a .png image and turns it into a 3D SLD array. It extends the image in
    the y direction where the sensitivity to scattering is low.

    '''
    def __init__(self, Dxyz, n, unit, value_list, step, mag_unit = None,
                 magVec = [None,None,None],inc_sub = [None,None],rawUnit = None,
                 ):

        self.Dxyz = Dxyz
        self.n = n
        self.unit = unit
        self.value_list = value_list
        self.step = step
        self.mag_unit = mag_unit
        self.magVec = magVec
        self.rawUnit = rawUnit
        self.inc_sub = inc_sub
        return

    def repeat(self,xy_repeat):
        '''
        **Overview:**

            Creates copies of the single unit cell array in the x-y direction as
            specified by the user.


        **Parameters:**

            *xy_repeat:* (int,[2]|count)
                The number of times the unit cell is repeated in the x and y
                direction (including the original unit).

        '''

        self.unit = concatenate([self.unit]*xy_repeat[0], axis=0)
        self.unit = concatenate([self.unit]*xy_repeat[1], axis=1)

        self.Dxyz[:2] = self.Dxyz[:2]*xy_repeat
        return


    def view(self):
        '''
        **Overview:**

            Outputs a 3D rendered viewing of the unit cell array using MayaVi.

        '''
        view.feature_plot(self.unit,self.n,self.step)
        return


    def mag_view(self):
        '''
        **Overview:**

            Outputs a 3D rendered viewing of the magnetic unit cell array
            using MayaVi.

        '''
        view.feature_plot(self.mag_unit,self.n,self.step)

        return

    def viewSlice(self):
        from .. viewers.wxzslice import unitView
        unitView(self.unit,
                 asarray([[0,self.Dxyz[0]],[0,self.Dxyz[1]],[0,self.Dxyz[2]]]),
                 self.step,self.n)


    def add_media(self):
        '''
        **Overview:**

            Adds a top and bottom layer to be the SLD of the incident medium and
            the substrate.


        .. Note::
            * This addition is important for the DWBA modeling where the
              calculation assumes that the top and bottom layer are
              semi-infinite.

        '''

        #initializes a sheet to be placed on the top and bottom of the matrix
        dummy_layer = ones([self.n[0],self.n[1],1])

        #creates the incident layer to be added
        dummy_layer *= self.inc_sub[0,0]

        #stacks the incident layer on top of the unit

        self.unit = concatenate([self.unit,dummy_layer],axis = 2)

        #creates the substrate
        dummy_layer = ones([self.n[0],self.n[1],1])
        dummy_layer *= self.inc_sub[1,0]

        #stacks the unit onto the substrate
        self.unit = concatenate([dummy_layer,self.unit],axis = 2)


        self.Dxyz[2] += self.step[2]*2

        self.n = [self.n[0],self.n[1],self.n[2]+2]

        self.value_list[2] = (arange(self.n[2]))*self.step[2]

        return

    def generateMIF(self,mifData = None,filename = None):

        import mif_creator
        mif_creator.create_mif(self,mifData,filename)
        return

class GrayImgUnit(object):
    '''
    **Overview:**

        This class creates a Unit_Cell object from a gray scale image by loading
        an image and the parameters that are correlated with that image.This
        class assumes that the direction into the picture is the same straight
        through. The user may choose this axis.


    **Parameters:**

        *newres* (float,[2]|count)
            Sometimes, the .png file is much higher in resolution than is needed
            for the scattering calculation. The user has the option of choosing
            a new resolution to down-scale the image file to.

    .. Note::
        * Currently, this only supports images that are colored to be on scale
          with the scattering length density values of the profile.

        * This file load system takes in .png files.

        * When loading the image, it assumes that the image is a three channel
          load but that each channel is equal. This can be much improved by
          counting channels and handling RGB files.


    '''

    def __init__(self,newres = None, filename = None):
        self.newres = newres
        if self.newres != None:
            image = imageUtil.imageload(filename = filename)[:,:,0]

            image = imageUtil.coursenRes(image,newres)

            self.image = flipud(image).T

        else:
            self.image = (flipud(imageUtil.imageload(filename = filename)
                                                                    [:,:,0]).T)
        return


    def unitBuild(self,Dxyz,scale,inc_sub = [None,None]):
        '''
        **Overview:**

            Producer method: This method produces a Unit_Cell object from the
            parameters that are carried over from the .omf loaded by the
            OOMMFUnit class. The object that this method produces is needed to
            work with the calculation API.


        **Parameters:**

            *Dxyz:* (float,[3]|angstroms)

                The real space length, width and height represented by the
                image. Because it is assumed that the image is the x z
                direction, the Dxyz[1] or Dy is meaningless here. It is left in
                only to allow for the consistency of Dxyz.


            *scale:* (float|factor)

                This parameter allows the user to uniformly scale the image
                values by a factor. This is to allow for directly scaling the
                image to SLD values.


        .. Note::
            * The user should be allowed to choose the axes on the image.

            * Some structure should be put together for manually assigning SLD
              values to colors on the map. This will only be needed if image
              loading becomes a highly used load system.

        '''
        YDIM = 8
        self.Dxyz = Dxyz
        self.scale = scale

        xdim,zdim = shape(self.image)
        dim = asarray([xdim,int(YDIM),zdim])
        d3image = zeros([xdim,1,zdim])
        d3image[:,0,:] = self.image

        self.d3image = asarray(repeat(d3image,dim[1],axis = 1))

        step = asarray([Dxyz[0]/dim[0], Dxyz[1]/dim[1], Dxyz[2]/dim[2]])

        self.inc_sub = asarray([[inc_sub[0],step[2],0.0],[inc_sub[1],
                                                        step[2],0.0]])

        valX = asarray(arange(dim[0])*step[0])
        valY = asarray(arange(dim[1])*step[1])
        valZ = asarray(arange(dim[2])*step[2])


        unitObject = Unit_Cell(Dxyz=Dxyz,n = dim,
                               unit = (self.d3image*self.scale),
                               value_list = [valX,valY,valZ],
                               step = step, inc_sub = self.inc_sub,
                               rawUnit = self)
        return unitObject

class OOMMFUnit(object):
    '''
    **Overview:**

        This class is used to create a Unit_Cell object from a .omf file output
        from the Object Oriented Micromagnetic Framework (OOMMF) software. It
        will produce a unit cell which has the structure unit array as well as a
        list of three arrays of the same size as the unit array which contains
        each of the three magnetic vector components. This information allows
        for the calculation of the four magnetic scattering cross-sections for
        the given system.


    .. Note::

        * Currently, this only supports systems that are constant in the
          z-direction.(eg. an SEM image of the feature is used as a mask to
          create a 2D image that OOMMF then calculates the minimized magnetic
          character for. The results are assumed to be consistent through the
          depth of the shape. This also only supports single feature unit cells.

    '''
    def __init__(self):
        omfImport = MIF.omf()
        omfImport.generate_coordinates()
        omfImport.generate_normalized_m()
        self.omfImport = omfImport
        return


    def unitBuild(self,SLD,inc_sub = [None,None]):
        '''
        **Overview:**

            Producer method: This method produces a Unit_Cell object from the
            parameters that are carried over from the .omf loaded by the
            OOMMFUnit class. The object that this method produces is needed to
            work with the calculation API.


        **Parameters:**

            *SLD:* (float|angstroms^2)
                When doing an OOMMF simulation, the software does not care about
                the structural SLD, however, to calculate the full off-specular
                scattering  this information is needed. At the time of Unit_Cell
                creation the user must specify the structural SLD for the
                magnetic feature being loaded by OOMMFUnit.


        .. Note::
            * Produces a Unit_Cell object.

            * This may not be the best place for the SLD input. This will have
              to be evaluated.

        '''
        if (self.omfImport.parameters['ValueRangeMinMag'] == 0.0):
            import warnings
            warnings.warn(
              'Non-magnetic material may be confused with no materials at all')

        unit = zeros(self.omfImport.dims)
        unit[self.omfImport.M.nonzero()] = SLD

        Dxyz = asarray([float(self.omfImport.parameters['xmax'])*1.0e10,
                      float(self.omfImport.parameters['ymax'])*1.0e10,
                      float(self.omfImport.parameters['zmax'])*1.0e10])


        value_list = (self.omfImport.node_x[:,0,0],
                      self.omfImport.node_y[0,:,0],
                      self.omfImport.node_z[0,0,:])


        step = asarray([float(self.omfImport.parameters['xstepsize']),
                      float(self.omfImport.parameters['ystepsize']),
                      float(self.omfImport.parameters['zstepsize'])])

        magVec = asarray([self.omfImport.mx,self.omfImport.mx,self.omfImport.mx])


        unitObject = Unit_Cell(Dxyz, unit,value_list,step, magVec,
                               inc_sub = [None,None], rawUnit = self.omfImport)

        return unitObject


class K3DUnit(object):
    '''
    **Overview:**

        This class is for a shape that is entered as a list of polygons and
        points from k3d modeling software.

    **Parameters**

        *filename(str):*
            The name of the file that was exported from the k3d model software.

        *SLD_list([])*
            The list of scattering length densities. there should be the same
            number of SLDs as there are shapes in filename(A^-2).

    '''
    def __init__(self,filename,k3d_scale = 1.0,SLD_list = [None]):

        shape_num_count = -1
        self.k3d_scale = k3d_scale
        self.SLD_list = SLD_list
        x_axis = 0
        y_axis = 1
        z_axis = 2

        openfile = open('C:/Python25/Py_Modeling/K3D data/'+filename +
                        '.txt','r')
        k3d_shapelist = []
        first = True
        counter = 1
        zero_test = 0
        for line in openfile.readlines():
            toks = re.match("([#]+)", line)
            if toks == None:
                if len(string.split(line))== 2:

                    if not first:
                        k3d_shapelist = self.k3d_listform(point_array,
                            poly_array,num_polygons,num_points,k3d_shapelist)

                    num_points, num_polygons = string.split(line)
                    num_points = int(num_points)
                    num_polygons = int(num_polygons)
                    point_array = empty([num_points,3])
                    poly_array = empty ([num_polygons,4])
                    first = False
                    b = 0
                    c = 0

                elif len(string.split(line))== 3:
                    x,y,z =  string.split(line)
                    point_array[b,x_axis] = float(x) * k3d_scale
                    point_array[b,y_axis] = float(y) * k3d_scale
                    point_array[b,z_axis] = float(z) * k3d_scale
                    b += 1

                elif len(string.split(line))== 4:
                    p_one,p_two,p_three,p_four = string.split(line)
                    poly_array[c,0] = int(p_one)
                    poly_array[c,1] = int(p_two)
                    poly_array[c,2] = int(p_three)
                    poly_array[c,3] = int(p_four)
                    c+=1


                elif len(string.split(line))== 0:
                    '''
                    Sometimes the last line is empty. this prevents an empty
                    cell from raising an error.
                    '''
                else:
                    print "ERROR: TO MANY INPUTS FROM ONE ROW"
        k3d_shapelist = self.k3d_listform(point_array,poly_array,num_polygons,
            num_points,k3d_shapelist)

        self.k3d_shapelist = k3d_shapelist

        return

    def k3d_listform(self,point_array,poly_array,num_polygons,num_points,
                     shapelist):
        '''
        **Overview:**

            This method forms the list of features needed to create the
            scattering matrix.

        **Parameters**

            *point_array*
                An array of points in 3d space that make up a feature.

            *poly_array*
                for numbers that represents the points that make up a polyhedran
                face.

            *num_polygons*
                total number of polygons that make up a feature

            *num_points*
                number of points that make up a feature

            *shapelist*
                the list of features that the new feature is being added to

        '''
        feature = self.K3D_Shape(vertices = point_array,edges = poly_array,
                    numpoly = num_polygons,numpoint = num_points)

        shapelist.append(feature)
        return shapelist


    class K3D_Shape_Collection:
        '''
        **Overview:**

            This contains the information for the description of multiple
            features

        **Parameters**

            *feature[x]*
                is an object of type Shape

            *description_list*
                a list of objects of type shape that holds the edges, vertices
                and Sof all of the features

        '''
        def __init__(self, description_list=None, correction_scaling = 1):
            self.description_list = description_list
            self.correction_scaling = correction_scaling



    class K3D_Shape:
        '''
        **Overview:**

            This contains variables to define a shape

        **Parameters**

            *vertices*
                the points in space that make up a shape

            *edges*
                Array containing the thicknesses of all layers in the substrate
                given in the order: Sample Bottom --> Feature/Substrate
                Interface

        '''
        def __init__(self,vertices=None,edges=None, numpoly=None,
                     numpoint=None):
            '''
            assigns the proper values to the class elements

            '''
            self.vertices = vertices
            self.edges = edges
            self.numpoly = numpoly
            self.numpoint = numpoint

            return

    def height(self):
        '''
        **Overview:**

            This module takes in a K3D_Shape object and determines its height in
            real space

        .. Note::
            * features is of type K3D_Shape

        '''
        max_value = 0
        for i in range(len(self.k3d_shapelist)):

            z_values = asarray(self.k3d_shapelist[i].vertices)
            max_tracker = max(z_values[:,2])
            if max_tracker > max_value:
                max_value = max_tracker
        return max_value

    def discritize(self,x_points,y_points,z_points,cell_to_fill, mag_to_fill):
        '''
        **Overview:**

            This method will turn a given K3D file of shapes into a numpy matrix
            of scattering length densities.

        .. Note::
            * This is done for a given Unit_Cell Object
        '''

        cell_to_fill = calculations.K3D_point_test(self,x_points,y_points,
                                                                z_points)
        return cell_to_fill


class GeomUnit(object):

    '''
    **Overview:**

        This is a producer of a Unit_Cell object. Given a Scene of Shape objects
        and other key parameters (defined below), this class will render a three
        dimensional numpy array of the SLD of the unit cell along with the
        magnetic SLD in the case the it is defined.


    **Parameters:**

        *Dxyz:* (float,[3])|(angstroms)
            The x,y and z real space size of the unit cell.

        *n:* (int,[3]|(count)
            The number of elements the x, y and z axis of the unit cell will be
            divided into. This is how course/fine the unit cell is discretized
            into.

        *scene:* (Scene)
            A scene object which holds the collection of shapes to be rendered
            into the unit cell. The renderer renders shapes in the order they
            are provided. Multiple shape are rendered by only filling unit cell
            array values where they have not yet been changed by previous
            shapes. For example, in the case of a core/shell scenario, the
            shapes should be listed so that the core is listed before the shell.

        *inc_sub:* (float,[2]|)
            [SLD of incident media,SLD of Substrate]. This holds the scattering
            length density for the incident and substrate media respectively.
            The attribute does not currently hold the neutron absorption of the
            materials which is negligible.


    **Parameters(Class):**

        *value_list:* (float,(3)|angstroms)

            This is a list of arrays for the x,y and z directions. Each array
            contains the real space distance of the array element from the
            origin. (eg. for 4 points at a step size of .2 angstroms in the x
            direction, value_list[0] = array([0,0.2,0.4,0.6])

        *unit:* (float:3D array| angstroms^2)

            This is the discretized representation of the structural unit cell.
            This is the array for which the scattering is calculated.

        *mag_unit:* (float:3D array|angstroms^2)

            This is the discretized representation of the magnetic unit cell.
            This is for the case of unpolarized neutrons where a value is given
            for the magnetic SLD which differs from the structural SLD.

        *step:* (float:[3]|angstroms)

            This is the reals space step size for the unit cell in the x, y and
            z direction. It is the total real space increment that a single
            array value represents.


    .. Note::
        * The renderer uses the mathematical formulas provided in the
          calculation.py file for each shape class to determine if each point in
          the 3D array falls within the shape.
        * If Dxyz[2] is not defined, the class chooses a value that
          will end just above the top of the tallest feature in the unit cell
          (adds approximately one layer of incident media).
        * Currently, the x,y and z values represent the real space value at
          the beginning of the discretized unit rather than the value of the
          unit at  the center. This must be revised.


    '''
    def __init__(self, Dxyz=[None,None,None], n=[None,None,None],scene = [None],
                 inc_sub = [None,None]):

        self.Dxyz = asarray(Dxyz)
        self.n = asarray(n)
        self.scene = scene


        if size(self.Dxyz) == 2:self.Dxyz = hstack((self.Dxyz,None))
        if self.Dxyz[2] == None:self.Dxyz[2] = scene.query_height()

        self.value_extend()
        self.unit = None
        self.mag_unit = None

        if (inc_sub == [None,None]):
            self.inc_sub = asarray([[0.0,self.step[2],0],[2.7e-6,0,0]])
        else:
            self.inc_sub = asarray([[inc_sub[0],self.step[2],0],
                                    [inc_sub[1],0,0]])
        self.render()

        return


    def value_extend(self):
        '''
        **Overview:**

            This module takes the individual values of step and length and
            creates a list of three arrays [x,y,z] that contains the real space
            value for each discrete piece of the unit cell array.

        '''


        #Measures at the top of the layer
        value_list = list([None]*3)

        self.step = (self.Dxyz/self.n)

        for i in range(3): value_list[i] =(arange(self.n[i]) * self.step[i])


        self.value_list = value_list
        return


    def render(self):
        '''
        **Overview:**

            Uses the discretized method contained in each Shape object in a
            Scene to create a 3D numpy array of SLD values which can be used to
            solve the scattering from.

        .. Note::
            *  This module fills the unit cell array with Shape objects in the
               order that they are listed in the Scene object. Shapes that are
               later in the list will write over those that came earlier. This
               means, for example,in the case of a core-shell sample, the outer
               shell object must be entered before the core.

        '''
        self.unit = zeros(self.n)
        self.mag_unit = zeros(self.n)

        for i in range(size(self.scene.scene_object)):

            sub_unit = zeros(shape(self.unit))

            mag_sub_unit = zeros(shape(self.unit))

            sub_unit,mag_sub_unit = self.scene.scene_object[i].discritize(
              self.value_list[0],self.value_list[1],self.value_list[2],
              sub_unit, mag_sub_unit)

            self.unit[sub_unit != 0.0] = sub_unit[sub_unit != 0.0]
            self.mag_unit[mag_sub_unit != 0.0] = (mag_sub_unit
                                                  [mag_sub_unit != 0.0])

        self.unit[self.unit == 0.0] = self.inc_sub[0,0]
        self.mag_unit[self.mag_unit == 0.0] = self.inc_sub[0,0]

        return


    def buildUnit(self):
        '''
        **Overview:**

            Producer method: This method produces a Unit_Cell object from the
            rendered geomUnit object. Because this is is the original
            development of Unit_Cell, the conversion is pretty simple and most
            of the parameters are in the exact form needed by Unit_Cell.

        .. Note::
            * The GeomUnit object must be rendered first. If it has not been
              rendered the method will do it automatically.

        '''

        if (self.unit == None):
            self.render()

        unitObject = Unit_Cell(self.Dxyz, self.n, self.unit, self.value_list,
                               self.step, self.mag_unit,
                               magVec = [None,None,None],inc_sub = self.inc_sub,
                               rawUnit = self)

        return unitObject

class Scene(object):
    '''

    **Overview:**

        This class is used to aggregate the different Shape objects that form a
        complete unit cell. It is used primarily by GeomUnit as a queue of
        objects that can be rendered into the unit cell array.


    **Parameters:**

        shapelist:(Shape|[]) = a list of Shape objects to be put into the scene.

    .. Note::
        * The Shapes may be entered as a list or as individual items.

    '''
    def __init__(self,shapelist=[]):
        #This if statement makes sure the incoming object is actually a Shape
        #object
        if isinstance(shapelist,Shape):
            shapelist = [shapelist]
        else:
            shapelist = list(shapelist)

        for s in shapelist:
            if not isinstance(s, Shape):
                raise TypeError('Scene needs a list of shapes')
        self.scene_object = shapelist
        return

        if not isinstance(shape, Shape):
            raise TypeError('Scene needs a shape')
        self.scene_object = list([None])
        self.scene_object[0] = shape

        return

    def add_element(self,element):

        '''
        **Overview:**

            Adds a shape object to the Scene object.This may be useful in the
            case where the Scene has been created and the user just wants to add
            one more Shape to it.


        **Parameters:**

            *element*
                A Shape object to add to a scene.
        '''
        self.scene_object.append(element)
        return

    def query_height(self,limit = 'max'):

        '''
        **Overview:**

            This method determines the maximum height of the Scene object (with
            an option to return the minimum if limit is given a value). This is
            not a thickness measurement so the highest Shape in the scene is not
            necessarily the thickest.


        Parameters:

            *limit:* (string)
                Allows the user to filter through the heights of the

        **Returns**

            Shape objects to return only the desired result.

            * 'max' returns the highest shape in the scene.

            * 'min' returns the lowest shape in the scene. This is the location
              of the top of the shape who's top is lowest in real space.

            * 'all' returns an array of all of the heights of all of the shapes.


            .. Note::
                * Each shape in a Scene has its own height method. When this
                  module is called, it searches through the heights of all of
                  the Shape objects in Scene to determine what the height of
                  the unit cell should be. It then records this value in Dxyz
                  as the z value.

        '''
        h_check = zeros(size(self.scene_object))
        for i in range(size(self.scene_object)):
            h_check[i] = self.scene_object[i].height()

        if limit == 'max':
            return max(h_check)
        elif limit == 'min':
            return min(h_check)
        elif limit =='all':
            return h_check
        else: raise "ERROR: Limit value is not an option"


    def query_center(self,limit = 'min'):

        '''
        **Overview:**

            This method returns the vertical component of the center points of
            the Shape objects in the Scene.


        **Parameters:**

            *limit:* (string)
                Allows the user to filter through the vertical center
                values of the Shape objects to return only the desired result.

        **Returns**

            * 'max' returns the highest center value in the scene.
            * 'min' returns the lowest center value in the scene.
              This is the location of the vertical component of the center
              values in real space.
            * 'all' returns an array of all of the centers of all of the shapes.


        '''
        center_collect = zeros([size(self.scene_object),3])

        for i in range(size(self.scene_object)):
            center_collect[i] = self.scene_object[i].center

        if limit == 'min':
            return min(center_collect,axis=0)

        elif limit == 'max':
            return max(center_collect,axis=0)

        elif limit == 'all':
            return center_collect
        else: raise "ERROR: limit value is not an option"



class Shape(object):

    '''
    **Overview:**

        *Abstract Class:* The different possible shape descriptions that can be
        used to build the unit cell.This class allows for the definition of
        shapes for a unit cell with different properties to be treated, on a
        fundamental level, as a geometric structure that can be added to a unit
        cell.

    .. Note::
        * Although all shapes are different, they all have the notion of a
          'bounding box', which is a Parallelapiped shape that can encompass the
          shape. Methods that treat the bounding box of a shape rather than the
          individual attributes that make the shape are held in this class.
    '''

    def on_top_of(self,element):

        '''
        **Overview:**
            This method alters the center value of a shape object so that the
            Shape who's on_top_of method was called is located on top of the
            Shape 'element'.


        **Parameter:**

            *element:* (Shape)
                The Shape object that is being put on top of the selected Shape
                object 'self'(The shapes who's z-component center value is being
                altered).

        '''
        if (element.center[:2] != [None]*2) and (self.center[:2] != [None]*2):
            self.center = copy.copy(element.center)
        self.center[2] = element.center[2] + (self.thickness()/2 +
                                              element.thickness()/2)
        return

    def is_core_of(self,element,offset = [0,0,0]):

        '''
        **Overview:**

            This method is used to place the Shape 'self' into the center of the
            Shape element. This creates a core shell type geometry.


        **Parameter:**

            *element:* (Shape)
                The Shape object that is being embedded into the selected Shape
                object 'self'(The shapes who's center value is being altered).

        .. Note::
            * In the case of use with a Layer Object, the Shape object will only
              be placed in the center in the z-direction because it's location
              in the x-y plane does not make a difference.

        '''
        if (element.center[:2] != [None]*2) and (self.center[:2] != [None]*2):
            self.center = copy.copy(element.center)
            for i in range(3): self.center[i] += offset[i]

        else:
            self.center[2] = copy.copy(element.center[2])
            self.center[2] += offset[2]
        return


class Sphere(Shape):
    '''
    **Overview:**

        Uses the generic formula for a sphere to create a sphere object.

    **Parameters(__init__):**

        *SLD:* (float|angstroms^2)
            The scattering length density of the sphere.

        *r:* (float|angstroms)
            The radius of the sphere.

        *center:* (float,[3]|angstroms)
            The x, y, and z component of the central point of the sphere. In the
            case that the center is set to [None,None,None] the shape will be
            put in the bottom corner of the unit cell (the bounding box will
            start at (0,0,0).

        *Ms:* (float|angstroms)
            The magnetic SLD of the material for this shape.

    '''

    def __init__(self, SLD, r = 1.0, center = [None,None,None], Ms = 0.0):

        self.r = r
        self.SLD = SLD
        self.Ms = Ms

        if (center == [None,None,None]):
            self.center = [r,r,r]
        else:
            self.center = center
        return

    def thickness(self):
        '''
        Overview:

            Returns the total thickness of the sphere.
        '''
        self.thickness = (2*self.r)
        return

    def height(self):
        '''
        Overview:

            Returns the total height of the sphere. This differs from thickness
            which only describes the thickness of the individual sphere whereas
            this method returns the maximum z-value of the shape in the unit
            cell.

        '''
        self.height = (self.center[2] + self.r)
        return

    def length(self):
        '''
        Overview:

            Returns the maximum length of the sphere (x direction)
        '''
        return self.r * 2.0

    def width(self):
        '''
        Overview:

            Returns the maximum width of the sphere (y direction)
        '''
        return self.r * 2.0


    def discritize(self,x_points,y_points,z_points,cell_to_fill, mag_to_fill):
        '''
        **Overview:**

            This module takes in x,y, and z points and fills the matrix array
            with the SLD of the shape for the points that fall within the shape

        **Parameters:**

            *x_points:* (float|angstroms)
                an array of x points to be determined if they fall within the
                sphere.

            *y_points:* (float|angstroms)
                an array of y points to be determined if they  fall within the
                sphere.

            *z_points:* (float|angstroms)
                an array of z points to be determined if they  fall within the
                sphere.

            *cell_to_fill:* (float,array|angstroms)
                This is the SLD matrix of the unit cell. It is filled by the
                render function.

            *mag_to_fill* (float,array|angstroms)
                This is the Ms matrix of the unit cell. It is filled by the
                render function.

        '''

        x_points =  reshape(x_points,[size(x_points),1,1])
        y_points = reshape(y_points,[1,size(y_points),1])
        z_points = reshape(z_points,[1,1,size(z_points)])

        cell_to_fill [calculations.sphere_point_test(self.center,self.r,
                                x_points,y_points,z_points)==True] = self.SLD


        mag_to_fill [calculations.sphere_point_test(self.center,self.r,
                                x_points,y_points,z_points)==True] = self.Ms
        return cell_to_fill, mag_to_fill

class RoundedParPip(Shape):
    '''
    **Overview:**
        It is rarely the case that a sample has totally sharp corners. This
        shape allows the user to determine the extent to which the corners are
        rounded.

    *Parameters(__init__):*
        *SLD:* (float|angstroms^2)
            The scattering length density of the sphere.

        *dim:* (float,[3]|angstroms)
            x, y and z dimensions of the feature.

        *center:* (float,[3]|angstroms)
            The x, y, and z component of the central point of the sphere. In the
            case that the center is set to [None,None,None] the shape will be
            put in the bottom corner of the unit cell (the bounding box will
            start at (0,0,0).

        *Ms:* (float|angstroms)
            The magnetic SLD of the material for this shape.

    '''
    def __init__(self, SLD, dim,  center = [None,None,None],curve = 0.0, Ms = 0.0):
        self.dim = dim
        self.SLD = SLD
        self.curve = 1-curve
        self.Ms = Ms

        if (center == [None,None,None]):
            self.center = [dim[0]/2,dim[1]/2,dim[2]/2]
        else:
            self.center = center

        PFy_one = sqrt(dim[1]**2 + (dim[0]/2.0)**2)
        PFy_two = sqrt((dim[0]/2.0)**2)


        self.ovalMin = asarray([dim[0]/2.0,dim[1]/2.0])

        self.ovalMax = asarray([sqrt(2.0)*self.dim[0]/2.0,sqrt(2.0)*self.dim[1]/2.0])

        self.ovalParam = (abs(self.ovalMax - self.ovalMin)*self.curve) + self.ovalMin

        #self.ovalParam = self.ovalMin
        return
    def thickness(self):
        '''
        **Overview:**

            Returns the total thickness of the parallelapiped.
        '''
        return self.dim[2]

    def height(self):
        '''
        **Overview:**

            Returns the total height of the parallelapiped. This differs from
            thickness which only describes the thickness of the individual
            sphere whereas this method returns the maximum z-value of the shape
            in the unit cell.
        '''
        return self.center[2] + self.dim[2]/2

    def length(self):
        '''
        **Overview:**

            Returns the maximum length of the parallelapiped (x direction)
        '''
        return self.dim[0]

    def width(self):
        '''
        **Overview:**

            Returns the maximum width of the parallelapiped (y direction)
        '''
        return self.dim[1]
    def discritize(self,x_points,y_points,z_points,cell_to_fill, mag_to_fill):
        '''
        **Overview:**

            This module takes in x,y, and z points and fills the matrix array
            with the SLD of the shape for the points that fall within the shape

        **Parameters:**

            *x_points* (float|angstroms)
                an array of x points to be determined if they fall within the
                parallelapiped.

            *y_points* (float|angstroms)
                an array of y points to be determined if they  fall within the
                parallelapiped.

            *z_points(float|angstroms)*
                an array of z points to be determined if they  fall within the
                parallelapiped.

            *cell_to_fill* (float,array|angstroms)
                This is the SLD matrix of the unit cell. It is filled by the
                render function.

            *mag_to_fill(float,array|angstroms)*
                This is the Ms matrix of the unit cell. It is filled by the
                render function.

        '''
        x_points =  reshape(x_points,[size(x_points),1,1])
        y_points = reshape(y_points,[1,size(y_points),1])
        z_points = reshape(z_points,[1,1,size(z_points)])


        cell_to_fill [calculations.parallel_point_test(self.center,
                       self.dim,x_points, y_points,z_points)==True] = self.SLD

        cell_to_fill [calculations.ellipse_point_test(self.center,
                       [self.ovalParam[0],self.ovalParam[1],self.dim[2]],
                       x_points, y_points,z_points)==False] = 0.0

        mag_to_fill [calculations.parallel_point_test(self.center,
                      self.dim,x_points,y_points,z_points)==True] = self.Ms
                      
        mag_to_fill [calculations.ellipse_point_test(self.center,
               [self.ovalParam[0],self.ovalParam[1],self.dim[2]],
               x_points, y_points,z_points)==False] = 0.0

        return cell_to_fill, mag_to_fill

class Parallelapiped(Shape):
    '''
    **Overview:**
        Uses the generic formula for a parallelapiped feature to create a
        parallelapiped object

    *Parameters(__init__):*
        *SLD:* (float|angstroms^2)
            The scattering length density of the sphere.

        *dim:* (float,[3]|angstroms)
            x, y and z dimensions of the feature.

        *center:* (float,[3]|angstroms)
            The x, y, and z component of the central point of the sphere. In the
            case that the center is set to [None,None,None] the shape will be
            put in the bottom corner of the unit cell (the bounding box will
            start at (0,0,0).

        *Ms:* (float|angstroms)
            The magnetic SLD of the material for this shape.

    '''
    def __init__(self, SLD, dim, center = [None,None,None], Ms = 0.0):


        self.dim = dim
        self.SLD = SLD
        self.Ms = Ms

        if (center == [None,None,None]):
            self.center = [dim[0]/2,dim[1]/2,dim[2]/2]
        else:
            self.center = center
        return

    def thickness(self):
        '''
        **Overview:**

            Returns the total thickness of the parallelapiped.
        '''
        return self.dim[2]

    def height(self):
        '''
        **Overview:**

            Returns the total height of the parallelapiped. This differs from
            thickness which only describes the thickness of the individual
            sphere whereas this method returns the maximum z-value of the shape
            in the unit cell.
        '''
        return self.center[2] + self.dim[2]/2

    def length(self):
        '''
        **Overview:**

            Returns the maximum length of the parallelapiped (x direction)
        '''
        return self.dim[0]

    def width(self):
        '''
        **Overview:**

            Returns the maximum width of the parallelapiped (y direction)
        '''
        return self.dim[1]


    def discritize(self,x_points,y_points,z_points,cell_to_fill, mag_to_fill):
        '''
        **Overview:**

            This module takes in x,y, and z points and fills the matrix array
            with the SLD of the shape for the points that fall within the shape

        **Parameters:**

            *x_points* (float|angstroms)
                an array of x points to be determined if they fall within the
                parallelapiped.

            *y_points* (float|angstroms)
                an array of y points to be determined if they  fall within the
                parallelapiped.

            *z_points(float|angstroms)*
                an array of z points to be determined if they  fall within the
                parallelapiped.

            *cell_to_fill* (float,array|angstroms)
                This is the SLD matrix of the unit cell. It is filled by the
                render function.

            *mag_to_fill(float,array|angstroms)*
                This is the Ms matrix of the unit cell. It is filled by the
                render function.

        '''
        x_points =  reshape(x_points,[size(x_points),1,1])
        y_points = reshape(y_points,[1,size(y_points),1])
        z_points = reshape(z_points,[1,1,size(z_points)])

        cell_to_fill [calculations.parallel_point_test(self.center,
                       self.dim,x_points, y_points,z_points)==True] = self.SLD

        mag_to_fill [calculations.parallel_point_test(self.center,
                      self.dim,x_points,y_points,z_points)==True] = self.Ms
        return cell_to_fill, mag_to_fill


class Ellipse(Shape):

    '''
    **Overview:**

        Uses the generic formula for an Ellipse feature to create a Ellipse
        object: :math:`(x^2/a^2) + (y^2/b^2) = 1`. The dim variable will be in
        the form [a,b,z]. This class can also be used to make a cylinder by
        setting dim[0] = dim[1]

    **Parameters(__init__):**
        *SLD:* (float|angstroms^2)
            The scattering length density of the Ellipse.

        *dim:* (float,[3]|angstroms)
            The 'a' component, 'b' component and thickness of the Ellipse
            respectively. 'a' is the radius of the Ellipse in the x direction
            and 'b' is the radius of the ellipsoid in the y direction.

        *center:* (float,[3]|angstroms)
            The x, y, and z component of the central point of the Ellipse. In
            the case that the center is set to [None,None,None] the shape will
            be put in the bottom corner of the unit cell (the bounding box will
            start at (0,0,0).

        *Ms:* (float|angstroms)
            The magnetic SLD of the material for this shape.

    .. Note::
        * This class is different than Ellipsoid which builds a lenticular
          shaped object where as this class produces a pillar shaped object.

    '''

    def __init__(self, SLD, dim, center = [None,None,None], Ms = 0.0):

        self.dim = [dim[0]/2,dim[1]/2,dim[2]]
        self.SLD = SLD
        self.Ms = Ms

        if (center == [None,None,None]):
            self.center = [self.dim[0],self.dim[1],dim[2]/2]
        else:
            self.center = center
        return

    def thickness(self):
        '''
        **Overview:**

            Returns the total thickness of the Ellipse.

        '''

        return self.dim[2]

    def height(self):
        '''
        **Overview:**

            Returns the total height of the ellipsoid. This differs from
            thickness which only describes the thickness of the individual
            Ellipse whereas this method returns the maximum z-value of the shape
            in the unit cell.

        '''
        return self.center[2] + self.dim[2]/2

    def length(self):
        '''
        **Overview:**

            Returns the maximum length of the Ellipse (x direction).

        '''
        return self.dim[0] * 2.0

    def width(self):
        '''
        **Overview:**

            Returns the maximum width of the Ellipse (y direction).

        '''
        return self.dim[1] * 2.0

    def discritize(self,x_points,y_points,z_points,cell_to_fill, mag_to_fill):
        '''
        **Overview:**
            This module takes in x,y, and z points and fills the matrix array
            with the SLD of the shape for the points that fall within the shape.

        **Parameters:**

            *x_points:* (float|angstroms)
                an array of x points to be determined if they fall within the
                Ellipse.

            *y_points:* (float|angstroms)
                an array of y points to be determined if they  fall within the
                Ellipse.

            *z_points:* (float|angstroms)
                an array of z points to be determined if they  fall within the
                Ellipse.

            *cell_to_fill:* (float,array|angstroms)
                This is the SLD matrix of the unit cell. It is filled by the
                render function.

            *mag_to_fill:* (float,array|angstroms)
                This is the Ms matrix of the unit cell. It is filled by the
                render function.

        '''

        x_points =  reshape(x_points,[size(x_points),1,1])
        y_points = reshape(y_points,[1,size(y_points),1])
        z_points = reshape(z_points,[1,1,size(z_points)])

        cell_to_fill [calculations.ellipse_point_test(self.center,
                        self.dim,x_points,y_points,z_points)==True] = self.SLD

        mag_to_fill [calculations.ellipse_point_test(self.center,
                        self.dim,x_points,y_points,z_points)==True] = self.Ms
        return cell_to_fill,mag_to_fill


class Cone(Shape):
    '''
    **Overview:**

        Uses the generic formula for a cone feature to create a cone object.
        Also allows for the creation of a truncated cone by providing a cut-off
        parameter for the thickness.

    **Parameters(__init__):**

        *SLD:* (float| angstroms^2)
            The scattering length density of the cone.

        *dim:* (float,[3]|angstroms)
            The x component, y component and thickness of the cone respectively.
            x is the radius of the cone base in the x direction and b is the
            radius of the cone base in the y direction.

        *stub:* (float|angstroms)
            Provides a hard cut-off for the thickness of the cone. this allows
            for the creation of a truncated cone object who side slope can be
            altered by using different z component values while keeping the stub
            parameter fixed.

        *center:* (float,[3]|angstroms)
            The x, y, and z component of the central point of the cone. In
            the case that the center is set to [None,None,None] the shape will
            be put in the bottom corner of the unit cell (the bounding box will
            start at (0,0,0).

        *Ms:* (float|angstroms)
            The magnetic SLD of the material for this shape.

    '''

    def __init__(self, SLD, dim, stub = None, center = [None,None,None],
                 Ms = 0.0):

        self.dim = [dim[0]/2,dim[1]/2,dim[2]]
        self.SLD = SLD
        self.stub = stub
        self.Ms = Ms

        if (center == [None,None,None]):
            self.center = [self.dim[0],self.dim[1],dim[2]/2]
        else:
            self.center = center

        return

    def thickness(self):
        '''
        **Overview:**

            Returns the total thickness of the cone.

        **Returns:**

            *thickness:* (float|angstroms)
                The total thickness of the cone object (absolute thickness).

        '''

        if self.stub == None: return self.dim[2]
        else: return self.stub

    def height(self):
        '''
        **Overview:**

            Returns the total height of the cone. This differs from thickness
            which only describes the thickness of the individual cone whereas
            this method returns the maximum z-value of the shape in the unit
            cell.

        **Returns:**

            *height:* (float|angstroms)
                The total height of the cone object (measures the top most part
                of the cone in z.)

        '''
        return self.center[2] + self.dim[2]/2

    def length(self):
        '''
        **Overview:**

            Returns the maximum length of the cone (x direction).

        **Returns:**

            *length:* (float|angstroms)
                The total length of the cone object (absolute distance in x)

        '''
        return self.dim[0] * 2.0

    def width(self):
        '''
        **Overview:**

            Returns the maximum width of the cone (y direction).

        **Returns:**

            *width:* (float|angstroms)
                The total width of the cone object (absolute distance in y)

        '''
        return self.dim[1] * 2.0

    def discritize(self,x_points,y_points,z_points,cell_to_fill, mag_to_fill):
        '''
        **Overview:**

            This module takes in x,y, and z points and fills the matrix array
            with the SLD of the shape for the points that fall within the shape

        **Parameters:**

            *x_points:* (float|angstroms)
                An array of x points to be determined if they fall within the
                cone.

            *y_points:* (float|angstroms)
                An array of y points to be determined if they  fall within the
                cone.

            *z_points:* (float|angstroms)
                An array of z points to be determined if they  fall within the
                cone.

            *cell_to_fill:* (float,array|angstroms)
                This is the SLD matrix of the unit cell. It is filled by the
                render function.

            *mag_to_fill*: (float,array|angstroms)
                This is the Ms matrix of the unit cell. It is filled by the
                render function.

        **Returns:**

            *cell_to_fill:* (array|angstroms)
                The discretized unit of the form factor built unit cell.


        '''


        x_points =  reshape(x_points,[size(x_points),1,1])
        y_points = reshape(y_points,[1,size(y_points),1])
        z_points = reshape(z_points,[1,1,size(z_points)])

        cell_to_fill [calculations.cone_point_test(self.center,
                self.dim,self.stub,x_points,y_points,z_points)==True] = self.SLD
        mag_to_fill [calculations.cone_point_test(self.center,
            self.dim,self.stub,x_points,y_points,z_points)==True] = self.Ms
        return cell_to_fill, mag_to_fill

class Pyrimid(Shape):

    '''
    **Overview:**
        Uses the generic formula for a Pyramid feature to create a Pyramid
        object.

    **Parameters(__init__):**
        *SLD:* (float|angstroms^2)
            The scattering length density of the Pyramid.

        *dim:* (float,[3]|angstroms)
            The x component, y component and thickness of the cone respectively.
            x is the length of the Pyramid base and y is the width of the
            Pyramid base.

        *stub:* (float|angstroms)
            provides a hard cut-off for the thickness of the Pyramid. this
            allows for the creation of a trapezoidal object who side slope can
            be altered by using different z component values while keeping the
            stub parameter fixed.

        *center:* (float,[3]|angstroms)
            The x, y, and z component of the central point of the Pyramid. In
            the case that the center is set to [None,None,None] the shape will
            be put in the bottom corner of the unit cell (the bounding box will
            start at (0,0,0).

        *Ms:* (float|angstroms)
            The magnetic SLD of the material for this shape.

    '''

    def __init__(self, SLD, dim, stub = None, center = [None,None,None],
                 Ms = 0.0):

        self.dim = [dim[0]/2,dim[1]/2,dim[2]]
        self.SLD = SLD
        self.stub = stub
        self.Ms = Ms

        if (center == [None,None,None]):
            self.center = [self.dim[0],self.dim[1],dim[2]/2]
        else:
            self.center = center

        return

    def thickness(self):
        '''
        **Overview:**

            Returns the total thickness of the Pyramid.

        '''
        if self.stub == None: return self.dim[2]
        else: return self.stub

    def height(self):
        '''
        **Overview:**

            Returns the total height of the Pyramid. This differs from thickness
            which only describes the thickness of the individual Pyramid whereas
            this method returns the maximum z-value of the shape in the unit
            cell.

        '''
        return self.center[2] + self.dim[2]/2

    def length(self):
        '''
        Overview:

            Returns the maximum length of the Pyramid (x direction)

        '''
        return self.dim[0] * 2.0

    def width(self):
        '''
        Overview:

            Returns the maximum width of the Pyramid (y direction)

        '''
        return self.dim[1] * 2.0

    def discritize(self,x_points,y_points,z_points,cell_to_fill, mag_to_fill):

        '''
        **Overview:**

            This module takes in x,y, and z points and fills the matrix array
            with the SLD of the shape for the points that fall within the shape

        **Parameters:**

            *x_points:* (float|angstroms)
                an array of x points to be determined if they fall within the
                Pyramid.

            *y_points:* (float|angstroms)
                an array of y points to be determined if they  fall within the
                Pyramid.

            *z_points:* (float|angstroms)
                an array of z points to be determined if they  fall within the
                Pyramid.

            *cell_to_fill:* (float,array|angstroms)
                This is the SLD matrix of the unit cell. It is filled by the
                render function.

            *mag_to_fill:* (float,array|angstroms)
                This is the Ms matrix of the unit cell. It is filled by the
                render function.

        '''

        x_points =  reshape(x_points,[size(x_points),1,1])
        y_points = reshape(y_points,[1,size(y_points),1])
        z_points = reshape(z_points,[1,1,size(z_points)])

        cell_to_fill [calculations.pyrimid_point_test(self.center,
                self.dim,self.stub,x_points,y_points,z_points)==True] = self.SLD

        mag_to_fill [calculations.pyrimid_point_test(self.center,
             self.dim,self.stub,x_points,y_points,z_points)==True] = self.Ms
        return cell_to_fill, mag_to_fill

class Layer(Shape):
    '''
    **Overview:**

        Creates an object that extends the length and width of the unit cell but
        is parameterized in the thickness direction.

    **Parameters(__init__):**

        *SLD:* (float|angstroms^2)
            The scattering length density of the Pyramid.

        *thickness_value:* (float|angstroms)
            The thickness of the layer.

        *center:* (float,[3]|angstroms)
            The x, y, and z component of the central point of the layer.
            Although allowed to be provided, the x and y component play no role
            in the layer location. the pertinent parameter here is only the z
            component.

        *Ms:* (float|angstroms)
            The magnetic SLD of the material for this shape.

    '''

    def __init__(self, SLD, thickness_value, center = [None,None,None],
                 Ms = 0.0):
        self.SLD = SLD
        self.thickness_value = thickness_value
        self.Ms = Ms

        if (center == [None,None,None]):
            self.center = [None,None,thickness_value/2]
        else:
            self.center = center

        return

    def thickness(self):
        '''
        **Overview:**

            Returns the total thickness of the layer.

        '''
        return self.thickness_value

    def height(self):
        '''
        **Overview:**

            Returns the total height of the layer. This differs from
            thickness which only describes the thickness of the individual layer
            whereas this method returns the maximum z-value of the shape in the
            unit cell.
        '''
        return self.thickness_value/2 + self.center[2]

    def discritize(self,x_points,y_points,z_points,cell_to_fill, mag_to_fill):
        '''
        **Overview:**
            This module takes in x,y, and z points and fills the matrix array
            with the SLD of the shape for the points that fall within the shape

        **Parameters:**

            *x_points:* (float|angstroms)
                an array of x points to be determined if they fall within the
                layer.

            *y_points:* (float|angstroms)
                an array of y points to be determined if they  fall within the
                layer.

            *z_points:* (float|angstroms)
                an array of z points to be determined if they  fall within the
                layer.

            *cell_to_fill:* (float,array|angstroms)
                This is the SLD matrix of the unit cell. It is filled by the
                render function.

            *mag_to_fill:* (float,array|angstroms)
                This is the Ms matrix of the unit cell. It is filled by the
                render function.

        '''

        cell_to_fill [:,:,calculations.layer_point_test(self.thickness_value,
        self.center[2]-(self.thickness_value/2.0),z_points)==True] = self.SLD


        mag_to_fill [:,:,calculations.layer_point_test(self.thickness_value,
        self.center[2]-(self.thickness_value/2.0),z_points)==True] = self.Ms
        return cell_to_fill,mag_to_fill

class Ellipsoid(Shape):
    '''
    **Overview:**

        Uses the generic formula for a Ellipsoid feature to create a Ellipsoid
        object. This object can be used to create a sphere by setting dim[0] =
        dim[1] = dim[2]

    **Parameters(__init__):**

        *SLD:* (float|angstroms^2)
            The scattering length density of the Ellipsoid.

        *dim:* (float,[3]|angstroms)
            The 'a' component, 'b' component and 'c' component of the Ellipsoid
            respectively. 'a' is the radius of the Ellipsoid in the x direction,
            'b' is the radius of the Ellipsoid in the y direction, and 'c' is
            the radius of the Ellipsoid in the z direction.

        *center:* (float,[3]|angstroms)
            The x, y, and z component of the central point of the ellipsoid. In
            the case that the center is set to [None,None,None] the shape will
            be put in the bottom corner of the unit cell (the bounding box will
            start at (0,0,0).

        *Ms:* (float|angstroms)
            The magnetic SLD of the material for this shape.

    .. Note::
        * This is a lenticular shaped object.

    '''

    def __init__(self, SLD, dim, center = [None,None,None], Ms = 0.0):
        self.SLD = SLD
        self.a = dim[0]
        self.b = dim[1]
        self.c = dim[2]
        self.Ms = Ms

        if (center == [None,None,None]):
            self.center = [self.a, self.b, self.c]
        else:
            self.center = center

        return

    def discritize(self,x_points,y_points,z_points,cell_to_fill, mag_to_fill):
        '''
        **Overview:**

            This module takes in x,y, and z points and fills the matrix array
            with the SLD of the shape for the points that fall within the shape

        **Parameters:**

            *x_points:* (float|angstroms)
                an array of x points to be determined if they fall within the
                Ellipsoid.

            *y_points:* (float|angstroms)
                an array of y points to be determined if they  fall within the
                Ellipsoid.

            *z_points:* (float|angstroms)
                an array of z points to be determined if they  fall within the
                Ellipsoid.

            *cell_to_fill:* (float,array|angstroms)
                This is the SLD matrix of the unit cell. It is filled by the
                render function.

            *mag_to_fill:* (float,array|angstroms)
                This is the Ms matrix of the unit cell. It is filled by the
                render function.

        '''
        x_points =  reshape(x_points,[size(x_points),1,1])
        y_points = reshape(y_points,[1,size(y_points),1])
        z_points = reshape(z_points,[1,1,size(z_points)])

        cell_to_fill [calculations.ellipsoid_point_test(self.center,
        self.a, self.b, self.c, x_points, y_points, z_points)==True] = self.SLD

        mag_to_fill [calculations.ellipsoid_point_test(self.center,
        self.a, self.b, self.c, x_points, y_points, z_points)==True] = self.Ms
        return cell_to_fill,mag_to_fill

    def thickness(self):
        '''
        **Overview:**

            Returns the total thickness of the layer.

        '''
        return self.c*2.0

    def height(self):
        '''
        **Overview:**
            Returns the total height of the layer. This differs from thickness
            which only describes the thickness of the individual layer whereas
            this method returns the maximum z-value of the shape in the unit
            cell.

        '''
        return self.c + self.center[2]

    def width(self):
        '''
        **Overview:**
            Returns the maximum width of the Pyramid (y direction)

        '''
        return self.b * 2.0



class Lattice(object):
    '''
    .. _LatticeClass:

    **Overview:**
        *Abstract Class:* This class is an abstract class which holds objects
        which describe the lattice structure of the repeating feature. These
        objects also hold the calculations required to determine the structural
        contribution to the scattering. Structure factors solved using these
        methods can be solved by integrating over the course q-spacings to
        reduce errors introduced by aliasing. Currently supported classes are:

            *Rectilinear:* A lattice structure that is spaced evenly in the x
            and y direction and is aligned with these directions.

                0|0|0|0|0

                0|0|0|0|0

                0|0|0|0|0

            *Hexagonal:* A lattice structure that is packed in a hexagonal
            ordering.

                _0|0|0|0|0|0

                0|0|0|0|0|0|0

                _0|0|0|0|0|0

    '''

    def rect_ft(self,Q, repeat_mod = [None,None,None]):
        '''
        **Overview:**

            Solves the structure factor for the Q points in the q_space object.
            Using this method solves the structure factor before integrating
            over the q steps. If used explicitly without integrating over the q
            steps aliasing errors can be introduced into the scattering.This is
            especially true where the q-step is course in the qx direction which
            can lead to mismatch in intensities between the negative and
            positive qx diffraction peaks.


        **Parameters:**

            *Q:* (q_space)
                a Q_space object.

            *repeat_mod:* (float,[3],count)
                A repeat modifier for the repeat attribute of a Lattice object.
                This is necessary when the effective repeat of a specific
                lattice type is different than the lattice spacing requested by
                the user.

        .. Note::
            * This calculation should be used in conjunction with integration
              over the x direction.

        '''

        if size(self.repeat) == 2:
            self.repeat = [self.repeat[0],self.repeat[1], 1]

        repeat = [None,None,None]
        for i in range(3):
            if repeat_mod[i] == None:
                repeat[i] = self.repeat[i]
            else:
                repeat[i] = repeat_mod[i]

        struc = ([None,None,None])

        struc[0] = zeros(shape(Q[0]), dtype = 'complex')
        struc[1] = zeros(shape(Q[1]), dtype = 'complex')
        struc[2] = zeros(shape(Q[2]), dtype = 'complex')

        struc[0] = ((sin((Q[0]*self.Dx*repeat[0])/2) / sin((Q[0]*self.Dx)/2)))
        struc[1] = ((sin((Q[1]*self.Dy*repeat[1])/2) / sin((Q[1]*self.Dy)/2)))
        struc[2] = ((sin((Q[2]*self.Dz*repeat[2])/2) / sin((Q[2]*self.Dz)/2)))

        struc[0][Q[0] == 0.0] = repeat[0]
        struc[1][Q[1] == 0.0] = repeat[1]
        struc[2][Q[2] == 0.0] = repeat[2]
        return struc


    def x_calc_sfx(self, qx):
        '''
        **Overview:**
            Used internally to solve the structure factor for a given qx value
            in the qx direction. This is possible because the qx, qy, and qz
            components are separable.

        **Parameters:**

            *qx:* (float)
                a qx value.

        .. Note::
            * This method is used in conjunction with the integration call to
              obtain the structure factor that is returned.

        '''
        if qx == 0:
            return self.repeat[0]
        else:
            return (sin((qx*self.Dx*self.repeat[0])/2.0) /sin((qx*self.Dx)/2.0))


    def x_gauss_sfx(self,qx, args):
        '''
        '''
        sigma = args[0]
        diff_space = args[1]
        orders = int(args[3])

        closePeak = round(qx/diff_space)*diff_space

        gauft = 0.0
        if qx == 0.0:
            return 1.0
        else:

            for m in range(-orders,orders+1):
                gauft += exp(-2 * pi * sigma**2 *(qx -
                                                  (closePeak+diff_space*m))**2)


        return gauft

    def gauss_normalize(self,args):
        sigma = args[0]
        diff_space = args[1]
        qspace = args[2]

        a = -diff_space/2.0
        b = diff_space/2.0
        avg = quad(self.x_gauss_sfx,a,b,args)[0]/diff_space

        return 1.0/avg

    def y_calc_sfx(self, qy):
        '''
        **Overview:**

            Used internally to solve the structure factor for a given qy value
            in the qy direction. This is possible because the qx, qy, and qz
            components are separable.


        **Parameters:**

            *qy:* (float)
                a qy value.

        .. Note::
            * This method is used in conjunction with the integration call to
              obtain the structure factor that is returned.

        '''
        if qy == 0:
            return self.repeat[1]**2
        else:
            return ((sin((qy*self.Dy*self.repeat[1])/2.0) /
                     sin((qy*self.Dy)/2.0)))

    def x_calc_sfx_shift(self, qx):
        '''
        **Overview:**

            Used internally to solve the structure factor for a given qx value
            in the qx direction. This solution applies a 0.5 phase shift to the
            wave solution which can be combined with the unshifted solution to
            give a solution to the scattering from a hexagonal lattice.


        **Parameters:**

            *qx:* (float)
                a qx value.

        .. Note::
            * This method is used in conjunction with the integration call to
              obtain the structure factor that is returned.

        '''
        if qx == 0:
            return self.nx**2 * exp(-1j*qx*(self.Dx/2.0))
        else:
            return ((sin((qx*self.Dx*self.repeat[0])/2.0) /
                     sin((qx*self.Dx)/2.0)))


    def y_calc_sfx_shift(self, qy):
        '''
        **Overview:**

            Used internally to solve the structure factor for a given qy value
            in the qy direction. This solution applies a 0.5 phase shift to the
            wave solution which can be combined with the unshifted solution to
            give a solution to the scattering from a hexagonal lattice.


        **Parameters:**

        *qy:* (float)
            a qy value.

        .. Note::
            * This method is used in conjunction with the integration call to
              obtain the structure factor that is returned

        '''
        if qy == 0:
            return self.ny**2 * exp(-1j*qy*(self.Dy/2.0))
        else:
            return ((sin((qy*self.Dy*self.repeat[0])/2.0) /
                     sin((qy*self.Dy)/2.0)))

    def phase_shift(self,Q):
        '''
        **Overview:**
            Used internally to solve a 0.5 phase shift for both the x and y
            directions.

        **Parameters:**

            *Q:* (q_space)
                a Q_space object.

        .. Note::
            * Used to superimpose two rectilinear lattice structures over each
              other.

        '''
        shift = exp(-1j*Q[0]*(self.Dx/2.0))*exp(-1j*Q[1]*(self.Dy/2.0))

        return shift


class Rectilinear(Lattice):
    '''
    **Overview:**

        A lattice structure that is spaced evenly in the x and y direction and
        is aligned with these directions. This is essentially a girded ordering.
        The ASCII art represents this structure.

    0 = feature

    ~ = spacing

            0~0~0~0~0

            0~0~0~0~0

            0~0~0~0~0

    **Parameters:**

        *Repeat:* (float,[3]|count)
            The number of repeats of the unit cell in the x, y and z direction.

        *Unit:* (Unit_Cell)
            a Unit_Cell object from which the unit cell length, width and height
            parameters can be obtained.

    .. Note::
        * This is a Lattice object and can uses methods from Lattice.

    '''


    def __init__(self,repeat,unit):

        self.Dx = unit.Dxyz[0]
        self.Dy = unit.Dxyz[1]
        self.Dz = unit.Dxyz[2]

        self.repeat = repeat
        return

    def theta_struc_calc(self,theta):
        qx,qz = theta.q_calc(5.0)
        qxlen = shape(qx)
        sigma = self.repeat[0]*self.Dx
        diff_space = 2.0*pi/self.Dx
        structure_factor = zeros(theta.points)

        #m = int(((sqrt(log(gaussCutoff)/(-2.0 * pi * sigma**2)))/Q.q_step[0]-.05)+1)

        orders = 7

        args = [sigma,diff_space,self.Dx,orders]
        norm =self.gauss_normalize(args)

        for i in range(qxlen[0]):
            for ii in range(qxlen[1]):
                structure_factor[i,ii] = self.x_calc_sfx(qx[i,ii])




        print 'STRUCTURE FACTOR WITH COHERENCE LENGTH CALCULATED'
        return structure_factor

    def gauss_struc_calc(self,Q,strucRefract = False):
        '''
        **Overview:**

            This structure calculation applies a gaussian convelution to the
            delta function diffraction peaks produced by the structure factor to
            produce a more accurate theory function. The convelution represents
            The combination of the diffraction from the lattice with the
            coherence length of the probing beam.

        **Parameters:**

            *Q:* (Q_space)
                The scattering produced by the structure factor of the model is
                calculated for the q range supplied by this Q_space object.

        '''
        gaussCutoff = 1.0e-16
        sigma = self.repeat[0]*self.Dx
        diff_space = 2.0*pi/self.Dx
        qspace = Q.q_step[0]
        qxlen = size(Q.q_list[0])
        qylen = size(Q.q_list[1])
        qzlen = size(Q.q_list[2])
        structure_factor = zeros(Q.points,dtype = 'complex')
        #m = int(((sqrt(log(gaussCutoff)/(-2.0 * pi * sigma**2)))/Q.q_step[0]-.05)+1)

        orders = 7.0
        args = [sigma,diff_space,self.Dx,orders]
        structure_factor = zeros(Q.points)
        norm =self.gauss_normalize(args)

        if  (strucRefract==False or Q.qx_refract == None):
            print 'Stucture Factor is not being refracted'
            for i in range(qxlen):
                a = Q.q_list[0][i] - (qspace/2.)
                b = Q.q_list[0][i] + (qspace/2.)
                avg = quad(self.x_gauss_sfx,a,b,args)[0]/qspace

                structure_factor[i,:,:] = avg*norm

        else:
            print 'Stucture Factor is being refracted'
            avg_refract = average(Q.qx_refract,axis = 1)

            for i in range(qxlen):
                print i
                for ii in range(qzlen):
                    a = avg_refract[i,ii] - (qspace/2.)
                    b = avg_refract[i,ii] + (qspace/2.)
                    avg = quad(self.x_gauss_sfx,a,b,args)[0]/qspace
                    structure_factor[i,:,ii] = avg*norm

        print 'STRUCTURE FACTOR WITH COHERENCE LENGTH CALCULATED'
        return asarray(structure_factor,dtype = 'complex')


    def struc_calc(self,Q):
        '''
        **Overview:**

            Returns the structure factor for a rectilinear lattice integrated
            over the qx steps.


        Parameters:

            Q:(q_space) = a Q_space object.

        .. Note::
            * The direct solution is calculated to get the y and z structural
              components. The integrated solution for the qx direction is then
              solved and applied over the direct solution. This is somewhat
              inefficient and can be streamlined.

        '''
        vectorized_q = Q.vectorize()
        struc = self.rect_ft(vectorized_q)
        qxlen = size(Q.q_list[0])
        qzlen = size(Q.q_list[2])

        qxstepsize = Q.q_step[0]

        if Q.qx_refract == None:
            for i in range(qxlen):
                qx = Q.q_list[0][i]
                avg = (quad(self.x_calc_sfx, qx - qxstepsize/2.,
                            qx + qxstepsize/2.)[0]/qxstepsize)
                struc[0][i,0,0] = avg

            structure_factor = struc[0] * struc[1] * struc[2]

        else:
            refractSF = zeros([Q.points[0],1,Q.points[2]])
            avg_refract = average(Q.qx_refract,axis = 1)
            for i in range(qxlen):
                print i
                for iii in range(qzlen):
                    qx = avg_refract[i,iii]

                    refractSF[i,0,iii] = (quad(self.x_calc_sfx,
                                               qx - qxstepsize/2.,
                                               qx + qxstepsize/2.)[0]
                                                / qxstepsize)

            yzstruc = struc[1] * struc[2]
            structure_factor = refractSF * yzstruc

        print 'STRUCTURE FACTOR CALCULATED'

        return structure_factor



class Hexagonal(Lattice):
    '''
    **Overview:**

        A lattice structure that is packed in a hexagonal ordering. This is
        produced by solving the rectilinear structure factor for two phases of
        repeating units and adding these phases together.


    0 = feature phase 1
    O = feature phase 2
    ~ = spacing

    Hexagonal:

    .0~0~0~0~0

    O~O~O~O~O~O

    .0~0~0~0~0

    O~O~O~O~O~O


    **Parameters:**

        *Repeat:* (float,[3]|count)
            The number of repeats of the unit cell in the x, y and z direction.

        *Unit:* (Unit_Cell)
            a Unit_Cell object from which the unit cell length, width and height
            parameters can be obtained.

    .. Note::
        * This is a Lattice object and can uses methods from Lattice.

    '''

    def __init__(self,repeat,unit):

        self.Dx = unit.Dxyz[0]
        self.Dy = unit.Dxyz[1]
        self.Dz = unit.Dxyz[2]

        self.repeat = repeat

        return


    def struc_calc(self,Q):
        '''
        **Overview:**

            This is the calculation of the structure factor for a hexagonal
            lattice. It returns the structure factor integrated over the qx
            direction by solving the scattering for two phases of rectilinear
            scattering and adds the results.

        **Parameters:**

            *Q:* (q_space)
                a Q_space object.

        '''

        vectorized_q = Q.vectorize()

        struc = (self.rect_ft(vectorized_q,[(2*self.repeat[0])+0.5,
                                            ( self.repeat[1])+0.5,None]))

        qxlen = Q.q_list[0]
        qxstepsize = Q.q_step[0]

        qylen = Q.q_list[1]
        qystepsize = Q.q_step[1]

        for i in range(size(qxlen)):

            qx = Q.q_list[0][i]
            avg = (quad(self.x_calc_sfx, qx - qxstepsize/2.0,
                        qx + qxstepsize/2.)[0] / qxstepsize)
            struc[0][i,0,0] = avg



        for i in range(size(qylen)):
            qy = Q.q_list[1][i]
            avg = (quad(self.y_calc_sfx, qy - qystepsize/2.0,
                        qy + qystepsize/2.)[0] / qystepsize)
            struc[1][0,i,0] = avg

        shift = self.phase_shift(vectorized_q)


        struc[0] = struc[0] + (struc[0] * shift)

        struc[1] = struc[1] + (struc[1] * shift)

        structure_factor = (struc[0] * struc[1]) * struc[2]

        return structure_factor

class Space(object):
    '''
    **Overview:**

        *Abstract Class* - This is a an object that holds the information about
        the space that the theory function is being calculated for.

    '''


class Theta_space(Space):
    r'''
    **Overview:**

        In some cases, it may be desirable to calculate the scattering from a
        model in theta space. This object acts like a Q_space object but the
        calculations are carried out in real space not reciprical space.

    **Parameters(__init__):**

        *minimums:* (float,[2]|angstroms)
            The minimum theta values that the user would like solved. The data
            is in the form: [:math:`\theta^{min}_{in}`,
            :math:`\theta^{min}_{out}`]

        *maximums:* (float,[2]|angstroms)
            The maximum theta values that the user would like solved. The data
            is in the form: [:math:`\theta^{max}_{in}`,
            :math:`\theta^{max}_{out}`]

        *points:* (float,[2]|angstroms)
            The number of points that the user would like to calculate for.
            (defined by the minimums and maximums) split into. The data is in
            the form: [:math:`\theta_{in}`, :math:`\theta_{out}`]

    **Parameters(Class):**

        *theta_step:* (float,[2]|degrees)
            Step size in :math:`\theta_{in}` and :math:`\theta_{out}`

        *theta_list:* (float,(2)[array]|degrees)
            The total list of values being solved for in :math:`\theta_{in}` and
            :math:`\theta_{out}`.

    '''
    def __init__(self,minimums, maximums, points):

        minimums = asarray([minimums[0],minimums[1]])

        maximums = asarray([maximums[0],maximums[1]])

        points = asarray([points[0],points[1]])

        self.type = 'Theta_space'
        self.minimums = minimums
        self.maximums = maximums
        self.points = points

        self.theta_list = [None]*2
        self.theta_step = ((asarray(maximums) - asarray(minimums))
                                                            /asarray(points))

        for i in range(2):
            self.theta_list[i] = linspace(self.minimums[i],self.maximums[i],
                                                                self.points[i])

        return


    def vectorize(self,type = 'float',unit = 'deg'):
        '''
        **Overview:**

            Turns the theta information given by a theta_space object into
            vectors to allow for vector math. Uses the numpy reshape
            functionality.


        **Parameters:**

            *type(str):*
                Allows the user to define the type of the numbers that theta is.
                (eg. float, complex)
        '''


        theta = asarray(self.theta_list[0].reshape(
                                              size(self.theta_list[0]),1),type)

        twotheta = asarray(self.theta_list[1].reshape(
                                              1,size(self.theta_list[1])),type)
        if unit == 'rad':
            theta,twotheta = radians(theta),radians(twotheta)

        return theta,twotheta


    def q_calc(self,wl):
        r'''
        **Overview:**

            This calculates the total Q vector based on the given theta values

        **Parameters**

            *wl:* (float|angstroms)
                The wavelength of the probing beam.

        **Return:**

            *q_vector:* (array|angstroms^-1)
                An array of Q vectors calculated for the combination of
                :math:`\theta_{in}` and :math:`theta_{out}` values for this
                object.

        '''

        theta,twotheta = self.vectorize()
        offset = theta - (twotheta/2.0)

        offset = radians(offset)
        twotheta = radians(twotheta)
        qx = zeros(self.points)
        qz = zeros(self.points)

        k0 = 2.0*pi/wl

        qx = 2.0 * (k0 * sin(twotheta / 2.0)) * sin(offset)
        qz = 2.0 * (k0 * sin(twotheta / 2.0)) * cos(offset)
        '''
        for i in range (self.points[0]):
            for ii in range (self.points[1]):
                qx[i,ii] = k0*(cos(twotheta[0,ii]/2.0) - cos(theta[i,0]))
                qz[i,ii] = k0*(sin(twotheta[0,ii]/2.0) + sin(theta[i,0]))

        '''
        return [qx,qz]


class Q_space(Space):
    '''
    **Overview:**

        Holds all of the information for the q-space output for which the
        scattering will be solved. Many of the attributes provided in this class
        make access to information about the scattering easier.


    **Parameters(__init__):**

        *minimums:* (float,[3]|angstroms)
            The minimum q values that the user would like solved. The data is in
            the form: [minimum x, minimum y, minimum z]

        *maximums:* (float,[3]|angstroms)
            The maximums q values that the user would like solved. The data is
            in the form: [maximums x, maximums y, maximums z]

        *points:* (float,[3]|angstroms)
            The number of points that the user would like the provided q space
            (defined by the minimums and maximums) split into. The data is in
            the form: [x points,y points,z points]


    **Parameters(Class):**

        *q_step:* (float,[3]|angstroms^-1)
            The reciprocal space step size for the x,y and z dimensions.

        *q_list:* (float,(3)[array]|angstroms^-1)
            The total list of values being solved for in the x, y and z
            directions.

        *q_refract:* (float,[array]|angstroms^-1)
            When the neutron beam is transmitted through a substrate, the beam
            refracts, altering the effective qx value. This is recorded in this
            variable. Its value is dependent on the ki and ko values for a
            specific qx,qy, qz combination.

        *k_space:* (float,[array]|angstroms)
            This is the equivelent k-space values for the given set of q values.
    '''
    def __init__(self,minimums, maximums, points):
        if size(minimums) == 2:
            minimums = [minimums[0],minimums[0],minimums[1]]

        if size(maximums) == 2:
            maximums = [maximums[0],maximums[0],maximums[1]]

        if size(points) == 2:
            points = [points[0],points[0],points[1]]

        self.type = 'Q_space'
        self.minimums = minimums
        self.maximums = maximums
        self.points = points
        self.q_list = [None]*3
        self.q_step = ((asarray(maximums) - asarray(minimums))/
                                                                asarray(points))
        self.qx_refract = None
        self.kin = None
        self.kout = None

        for i in range(3):
            self.q_list[i] = linspace(self.minimums[i],self.maximums[i],
                                      self.points[i])
        return

    def vectorize(self,type = 'float'):
        '''
        **Overview:**

            Turns the q information given by a q_space object into vectors to
            allow for vector math. Uses the numpy reshape functionality.


        **Parameters:**

            *type(str):*
                Allows the user to define the type of the numbers that q is.
                (eg. float, complex)
        '''

        q = [None,None,None]
        q[0] = asarray(self.q_list[0].reshape(size(self.q_list[0]),1,1),type)
        q[1] = asarray(self.q_list[1].reshape(1,size(self.q_list[1]),1),type)
        q[2] = asarray(self.q_list[2].reshape(1,1,size(self.q_list[2])),type)
        return q

    def normalize(self):
        '''
        **Overview:**

            Creates 3 arrays which contain the qx, qy, and qz value which are
            normalized by the total q magnitude.

        **Returns:**

            (list,3D array|angstroms^-1) The normalized Q values.
        '''
        qx,qy,qz = self.vectorize()

        magQ = sqrt(qx**2 + qy**2 + qz**2)

        qxn = qx/magQ
        qyn = qy/magQ
        qzn = qz/magQ

        return [qxn,qyn,qzn]

    def getExtent(self):
        '''
        **Overview:**

            This method is used to get the minimum and maximum plot area of
            the Q_space object which can be directly fed to a pylab plotting
            object.

        **Returns:**
            (array|angstroms^-1)
                Returns an array in the form
                :math:`[Q^{min}_{x},Q^{max}_{x},Q^{min}_{z},Q^{max}_{z}]`

        '''
        return asarray([self.minimums[0],self.maximums[0],
                        self.minimums[2],self.maximums[2]])

    def getKSpace(self,wavelength):
        '''
        **Overview:**

            This method creates an attribute which holds the equivalent k-space
            values for a given set of Qs.

        **Returns:**
          (array|angstroms)

        '''
        from ..theory.approximations import QxQyQz_to_k
        vecQ = self.vectorize()
        self.kin,self.kout = QxQyQz_to_k(vecQ[0],vecQ[1],vecQ[2],wavelength)
        return


class Beam(object):
    '''
    **Overview:**

            Hold the beam information. These are all of the instrument
            characteristics the have an effect on the scattering.


    **Parameters(__init__):**

        *wavelength:* (float|angstroms)
            For reactor source, the wavelength is used to calculate the
            resolution of the instrument.

        *angluar_div:* (float|degrees)
            The angular divergence of the beam

        *background:* (float|intensity)
            This is the dark counts on the detector

        *resolution:* (float)
            Generally, spallation sources have a resolution that
            they use as a beam parameter.

    .. Note::
              * This class is primarily developed for a reactor source but is
                open to parameters needed for a spallation source.


    '''

    def __init__(self,wavelength = None,angular_div = None,background = None,
                 wavelength_div = None,resolution =None):

        if (wavelength != None) and (wavelength_div != None):
            self.resolution = wavelength_div/wavelength
        elif resolution != None:
            self.resolution = resolution
        else:
            raise "Error: No resolution information given"
        self.angular_div = angular_div
        self.wavelength_div = wavelength_div
        self.background = background
        self.wavelength = wavelength
        return


def _test():

    '''
    a = OOMMFUnit()
    oommfUnit = a.unitBuild(9.8e-6)
    oommfUnit.view()
    return
    '''

    '''
    a = GrayImgUnit()
    '''

    '''
    # test the theta object
    a = Theta_space([0.0,0.0],[4.0,8.0],[100,100])
    #a.vectorize()
    #k = a.vectorize_k(5.0)
    q = a.q_calc(5.0)
    import pylab

    extent = [a.minimums[1],a.maximums[1],a.minimums[0],a.maximums[0]]

    pylab.title('qx')
    pylab.xlabel('2theta')
    pylab.ylabel('theta')
    pylab.imshow(fliplr(q[0]),extent = extent)
    pylab.colorbar()
    pylab.figure()
    pylab.title('qz')
    pylab.xlabel('2theta')
    pylab.ylabel('theta')
    pylab.imshow(fliplr(q[1]),extent = extent)
    pylab.colorbar()
    pylab.show()
    '''

    '''
    space = Q_space([-.0001,-0.001,0.00002],[.0001,0.001,0.04],[50,5,50])
    nqx,nqy,nqz = space.normalize()
    print shape(nqx)
    '''
    #Au = Ellipse(SLD = 4.506842e-6,dim=[3.75e4,3.75e4,600.0])

    Au = RoundedParPip(SLD = 4.506842e-6,dim=[3.75e4,3.75e4,600.0],curve = 0.75)

    #Au = Ellipse(SLD = 4.506842e-6,dim=[5.0e4,5.0e4,550.0])
    Cr = Layer(SLD = 3.01e-7,thickness_value = 27.0)
    Au.on_top_of(Cr)
    scene = Scene([Au,Cr])

    GeoUnit = GeomUnit(Dxyz = [10.0e4,10.0e4,700.0], n = [100,100,100],scene = scene, inc_sub = [0.0,2.7e-6])

    unit = GeoUnit.buildUnit()
    #unit.viewSlice()


if __name__=="__main__":_test()
