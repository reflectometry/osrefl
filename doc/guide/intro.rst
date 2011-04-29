*************
Introduction
*************
This is an instruction manual on how to use the current infrastructure developed to model off-specular neutron scattering data. The manual covers the basic attributes of the software as well as how to build and model a specific sample.


To Begin...
############
This documentation provides instructions on how to write a simple model script. The first requirement is the import statements. Only two import statements are needed for the script to run. They should look like:
::

	import scatter, sample_prep

sample_prep holds all of the code for creating a model. Scatter hold all of the information about the different approximations that can be used.



Creating a Unit Cell
#####################
This section reviews how to create a model to be scattered off of. There are for main model creation tools; GeomUnit, K3D, OOMFUnit, and GrayImgUnit. Each of these can be used to produce a discretized unit cell.

GeomUnit
**********
GeomUnit uses a mathematical description of the shapes in a sample to produce the unit scale to be scattered from. The first step is to create the shapes that are in the sample. A full list is documented in the code and an example is show below:
::

	Au = Ellipse(SLD = 4.506842e-6,dim=[3.75e4,3.75e4,700.0])
	Cr = Layer(SLD = 3.01e-6,thickness_value = 20.0)

Here, a gold ellipse and a chrome adhesion layer are produced.

The modelling software allows the user to specify a centre value for manipulation of shape location in 3D space. To make the shape localisation easier, some tools have been created to orient shapes relative to other shapes. For example:
::

	Au.on_top_of(Cr)

Will place the centre of the gold feature such that the ellipse base is flush with the adhesion layer. Other tools like this help orient the model and are explicitly defined in the documentation. 

Once the shapes are created and oriented appropriately, they can be added to a Scene. A Scene is a container class which holds all of the shapes that make up a model. By allowing for the addition of an arbitrary amount of shapes, arbitrarily complex systems can be produced. The Scene is simply created by:
::

	scene = Scene([Au,Cr])

Now we can produce the GeomUnit object. This objects contains the rest of the pertinent unit cell information such as dimension and discritization count:
::

	GeomUnit = GeomUnit(Dxyz = [10.0e4,10.0e4, 1000.0], n = [50,51,52],scene = scene)

Finally, we need to run a producer command that will tie the GeomUnit object to the infrastructure that handles all discretized unit cells:
::

	unit = GeoUnit.buildUnit()


GrayImagUnit
*************
A unit cell can also be created using the GrayImagUnit. This loader inputs a grey scale .png image file whose grey scale values are related to the SLD of that layer. When an object is created:
::

	a = GrayImgUnit()

A file loader will open asking the user to choose the images file. The file name may also be scripted into the call:
::

	img = sample_prep.GrayImgUnit(newres = numpy.array([150,400]),filename ='/Downloads/sample1_sld.png')

Once the object has been created, the universal 'Unit' must be created. For this, the software needs to know the rest of the unit cell information such as unit cell dimensions, discretizeation count and image scaling factor:
::
 
	unit = img.unitBuild(Dxyz = [8480.0,8480.0,3500.0], scale = 1.0e-5,inc_sub=[0.0,2.0784314e-6])

.. Note::
	* This unit building method assumes the image is extended infinity in the y direction which is the direction into the image, ie. the image is of the x-z plane of the sample and the direction into the image is y.

K3DUnit
********
This unit is created from the `K-3D software <http://www.k-3d.org/>`_. This software allows an output file that contains a list of points and plains that make up the shapes in the 3D model. This loader pares through these shapes using a point tracer method to determine whether or not a point falls inside the polyhedron. Although slow and limited in its modelling capability relatively complicated structures can be created easily using this method.

OOMMFUnit
**********
This unit loader creates a magnetic sample using the magnetic minimization software call `Object Oriented MicroMagnetic Framework <http://math.nist.gov/oommf/>`_. This allows for both the flexibility of a dicritized system with an simple way to produce magnetic structures.



Creating a Model
#################
A unit is only one piece of the information needed to produce a scattering model. The model must also have a Lattice which contains the information about the repeat structure:
::

	lattice = Rectilinear([20,20,1],unit)

A Q_space object which tells the model where to calculate the scattering in reciprocal space:
::

	q_space = Q_space([-.0001,-0.001,0.00002],[.0001,0.001,0.04],[200,50,200])

and a Beam object which provides the model with information about the probing beam:
::

	beam = Beam(5.0,None,None,0.05,None)

Once these objects are created they can be combined to form a Calculator object. This class is made to:

 * Ensure that the user has provided all of the necessary pieces to calculate the scattering.

 * Makes calculating scattering using different theories convenient.

This is created by:
::

	sample = Calculator(lattice,beam,q_space,unit)


Theory Function
#################
Now that the software has everything it needs to calculate off-specular scattering, a modelling formalism must be chosen. The option here can be found elsewhere in the documentation but the modelling itself is easily run by the convention:
::

	sample.BA()

Each theory calculation is a method on the calculator object. The user can now specify if they would like to run a resolution correction on the sample. This is done by:
::

	sample.resolution_correction()



Viewing
########
To view the scattering, the user simply needs to script:
::

	sample.view_uncorrected()

to view the uncorrected scattering or:
::
 
	sample.view_corrected()

for the corrected data. Because there is no set convention for what the user will want to view, the script must have:
::

	show()

to view the output plots.

