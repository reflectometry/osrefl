*************
Introduction
*************
This is an instruction manual on how to use the current infrastructure developed to model off-specular neutron scattering data. The manual covers the basic attributes of the software as well as how to build and model a specific sample.

Installing the Software
########################

There are many scientific libraries which are needed to run this code. All of the libraries are free and can easily be installed simultaniously by going to `pythonxy <http://www.pythonxy.com/>`_ and installing their product. In addition, if a Cuda compatable Nvidia GPU device is available, pyCuda must also be installed which may be downloaded at `pycuda <http://mathema.tician.de/software/pycuda>`_. Once these packages are installed, the osrefl package needs to be installed.

The software package may be downloaded at:
http://www.reflectometry.org/danse/software.html

A link to the source code may also be found on this site.

Windows:
	For a windows install, download the osrefl executable, double click to run it, and follow the on screen instructions.

Linux:
	Download the source code from the link on the reflectometery.org website or go to:

http://danse.us/trac/reflectometry/browser/trunk/osrefl

Once the software is in the desired location, open up a command line prompt (either a terminal in linux or command prompt in windows). Change the directory to the top level (first) osrefl folder and type:
::

	python setup.py build_ext --inplace

This will build the C code in the appropriate place and allow the software to run.

Running Examples
########################

The examples included in this software are located in the top level osrefl folder in a folder called 'examples'. To run the default 
example, change the directory to the top level osrefl folder and use the command:
::

    python osrefl.py

This will run AuFit.py which is the model for an Au pillar system. This can easily be changed for any of the example file by:
::
    python osrefl.py examples/<filename>
    
For example:
::
    python osrefl.py examples/BA_demo.py

This is also the procedure for running your own modeling scripts.

To Begin Modeling
########################

This documentation provides instructions on how to write a simple model script. The first requirement is the import statements. Only two import statements are needed for the script to run. They should look like:
::

	import scatter, sample_prep

sample_prep holds all of the code for creating a model. Scatter hold all of the information about the different approximations that can be used.


Creating a Unit Cell
######################

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

	sample.viewUncor()

to view the uncorrected scattering or:
::
 
	sample.viewUncor()

To view both the corrected and uncorrected plots side-by-side use:

::

	sample.viewCorUncor()

to view the output plots.

Modeling Data
#################
In the examples folder is a python script called AuFit.py. This is an example of how to compare a fit to real data using this software. This will go through the steps taken in this file.

Data Loading
****************
First, a model must be created as was shown in the previous section. The data included for this example was taken from Au pillars on a Si/Cr substrate. The data loading is all completed through GUI interfaces and only requires one line of code in the script. First, the data is loaded using the:

::

	Au_measurments = Data()

call which is found in the osrefl.loaders.andr_load module.

This call will bring up a file selector where multiple .cg1 data files may be loaded and combined. Use the "Choose input files" button to select the files. There is no error checking here to make sure the files are combined nicely so be sure that the selected data files are actually related to a single measurement. "The Main Beam center pixel" button is not used here. Hit the "Save and exit" button. Next, a screen will open to convert the data into qx and qz space plots. Enter the Qx range and the number of points to convert the x axis and the Qz range and number of points to convert the y axis. The X pixel value for Q=0 is the pixel on the detector for which Q=0 and is important for proper conversion. A good check for this is to view the resulting Q plot. The specular scattering should be straight along the Qx=0 line. If it starts to bend at high Qz values, then rerun the script and adjust the value accordingly.
	The next window will be the data selection window. This allows the user to select a specific subset of their data to model. This is important as modeling can be long and areas that don't have data should not have models calculated for it.

Model Building
*******************

	The models are build in the sample way as described in the model building section of this manual. One key additional command that is useful is:

::

	q_space = Au_measurments.space

This command takes the q space values and point count from the selected data q space and uses it as the points to solve the model for. This is convenient for calculating models in the most effcient manner.

Model/Data Interactor
*************************

	There is now a view and GUI interactor for the data and model. This can be used by:

::

	test_data.fitCompare(Au_measurments,titles = ['data','Model Label'])

where the method is run on the model and given the data as a parameter. Other options can be found in the method description in this documentation.

