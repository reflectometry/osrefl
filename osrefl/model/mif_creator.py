# Copyright (C) 2008 University of Maryland
# All rights reserved.
# See LICENSE.txt for details.
# Author: Christopher Metting

#Starting Date:8/11/2010
from numpy import shape,size

class mifData(object):
    def __init__ (self,H_strt_xyz = [0.0,500.0,0.0],H_end_xyz = [0.0,0.0,0.0],
                                                        steps = 5):
        self.H_strt_xyz = H_strt_xyz
        self.H_end_xyz = H_end_xyz
        self.steps = steps

        self.fieldCount = size(self.steps)
        return

def create_mif(unit, mifDataObj = None, filename = None):
    '''
    **Overview:**

        This module creates a mif file of the unit cell so that the magnetic
        character can be solved for in the oommf software. Currently, it only
        inputs the Unit cell information and other magnetic properties like
        magnetic field and demag are not available. This will be added as the
        become needed.

    **Parameters:**

        *filename* (str,filename)
            The location that the user would like to save the .mif file to.

    '''

    if filename == None:
        filename= '/tmp/temp_mif_one.mif'
    from numpy import size
    nofile = unit.mag_unit.flatten()

    if (mifDataObj == None):
        mifDataObj = mifData()

    f = open(filename,'w')
    print >>f, '#MIF 2.1'
    print >>f, '#Mif file created by off-specular modeling software'
    print >>f, '#Description: this file was made by taking a sample that'
    print >>f, '#was created by the off-spec modeling software and turning'
    print >>f, '#the formula used to discritize into a mif file'
    print >>f, ''
    print >>f, 'set pi [expr 4*atan(1.0)]'
    print >>f, 'set mu0 [expr 4*$pi*1e-7]'
    print >>f, ''

    idx = 0
    print >>f, 'set ::data {'
    for i in range(unit.n[2]):
        for i in range(unit.n[1]):
            print>>f, ''
            for i in range(unit.n[0]):
                print >>f,'%.5e' %nofile[idx],'     ',
                idx += 1.0

    print >>f, '   '
    print >>f, '}'
    print >>f, ''
    print >>f, ''
    print >>f, 'RandomSeed 1'
    print >>f, ''
    print >>f, 'Specify Oxs_BoxAtlas:atlas {'
    print >>f, '    xrange {0 ' + str(unit.Dxyz[0]*(1.0e-10)) + '}'
    print >>f, '    yrange {0 ' + str(unit.Dxyz[1]*(1.0e-10)) + '}'
    print >>f, '    zrange {0 ' + str(unit.Dxyz[2]*(1.0e-10)) + '}'
    print >>f, '}'
    print >>f, ''

    print >>f, 'Specify Oxs_RectangularMesh:mesh [subst {'
    print >>f, '    cellsize {' + str(unit.step[0]*(1.0e-10)),
    print >>f, str(unit.step[1]*(1.0e-10)), str(unit.step[2]*(1.0e-10))+ '}'
    print >>f, '    atlas :atlas'
    print >>f, '}]'
    print >>f, ''

    print >>f, 'Specify Oxs_UZeeman [subst {'
    print >>f, '    multiplier [expr 0.001/$mu0]'
    print >>f, '    Hrange {'

    if mifDataObj.fieldCount == 1:
        print >>f, '        {',mifDataObj.H_strt_xyz[0],
        print >>f, mifDataObj.H_strt_xyz[1],mifDataObj.H_strt_xyz[2],
        print >>f, mifDataObj.H_end_xyz[0],mifDataObj.H_end_xyz[1],
        print >>f, mifDataObj.H_end_xyz[2],mifDataObj.steps,'}'

    else:
        for i in range(mifDataObj.fieldCount):
            print >>f, '        {',mifDataObj.H_strt_xyz[i][0], mifDataObj.H_strt_xyz[i][1],mifDataObj.H_strt_xyz[i][2], mifDataObj.H_end_xyz[i][1],mifDataObj.H_end_xyz[i][1],mifDataObj.H_end_xyz[i][2],mifDataObj.steps[i],'}'

    print >>f, '    }'
    print >>f, '}]'
    print >>f, ''

    print >>f, 'Specify Oxs_Demag {}'
    print >>f, ''

    print >>f, 'Specify Oxs_CGEvolve:evolve {}'
    print >>f, ''

    print >>f, 'Specify Oxs_MinDriver {'
    print >>f, '    basename temp_mif'
    print >>f, '    evolver :evolve'
    print >>f, '    stopping_mxHxm 0.1'
    print >>f, '    mesh :mesh'
    print >>f, '    Ms { Oxs_ScriptScalarField {'
    print >>f, '        atlas :atlas'
    print >>f, '        script OffSpecFormula'
    print >>f, '        script_args relpt'
    print >>f, '    }}'
    print >>f, 'm0 { Oxs_RandomVectorField {'
    print >>f, '         min_norm 1.0'
    print >>f, '         max_norm 1.0'
    print >>f, '      }}'
    print >>f, '}'
    print >>f, ''

    print >>f, 'proc OffSpecFormula {x y z} {'
    print >>f, ''

    print >>f, '    set nx',str(unit.n[0])
    print >>f, '    set ny',str(unit.n[1])
    print >>f, '    set nz',str(unit.n[2])

    print >>f, ''
    print >>f, '    set xindex [expr {int(floor($x * $nx))}]'
    print >>f, '    set yindex [expr {int(floor($y * $ny))}]'
    print >>f, '    set zindex [expr {int(floor($z * $nz))}]'
    print >>f, ''
    print >>f, '    set idx [expr {($xindex*$ny*$nz + $yindex*$nz + $zindex)}]'
    print >>f, '    set Ms [lindex $::data $idx]'
    print >>f, '    return $Ms'
    print >>f, '}'
    print >>f, 'Destination archive mmArchive'
    print >>f, 'Schedule DataTable archive Stage 1'
    print >>f, 'Schedule Oxs_MinDriver::Magnetization archive Stage 1'
    f.close()

    print 'Mif File Created at: ' + str(filename)
    return

def _test():

    import sample_prep
    Au = (sample_prep.Parallelapiped(SLD = 4.506842e-6,
                                     Ms = 8.6e5,dim=[5.0e4,5.0e4,2.0e4]))

    Cr = (sample_prep.Layer(SLD = 3.01e-6,Ms = 7.7e8,
                            thickness_value = 1000.0))

    #Au.on_top_of(Cr)
    scene = sample_prep.Scene([Au])

    GeoUnit = (sample_prep.GeomUnit(Dxyz = [10.0e4,10.0e4,2.2e4],
                                    n = [10,10,10],scene = scene))

    unit = GeoUnit.buildUnit()
    mif = mifData()
    #mif = mifData([[0,500,0],[0,500,0],[0,500,0]],[[0,500,0],[0,500,0],[0,500,0]],[1,2,3])
    unit.generateMIF(mif)
if __name__=="__main__":_test()
