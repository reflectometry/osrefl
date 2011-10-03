# This program is in the public domain
# Authors: Paul Kienzle, Christopher Metting
#03/23/2010

import time
import Queue
import threading

from pycuda import gpuarray
import pycuda.driver as cuda
from pycuda.compiler import SourceModule
import numpy

from . import approximations,smba_wave_driver

def magCudaBA_form(cell,Q,lattice,beam,omf,precision='float32', refract = True):
    nqx = numpy.shape(Q.q_list[0])[0]
    nqy = numpy.shape(Q.q_list[1])[0]
    nqz = numpy.shape(Q.q_list[2])[0]

    if refract == True:
        stack = cell.inc_sub
        psi_in_one,psi_in_two,psi_out_one,psi_out_two,qx_refract = (
                    smba_wave_driver.wave(stack, Q.q_list[0], Q.q_list[1],
                             Q.q_list[2],beam.wavelength,cell.step[2],
                             precision=precision))

        Q.qx_refract = qx_refract
    else:
        qx_refract = Q.q_list[0]
        qx_refract = numpy.reshape(qx_refract,[numpy.size(qx_refract),1,1])
        qx_refract = numpy.repeat(qx_refract,nqy,axis = 1)
        qx_refract = numpy.repeat(qx_refract,nqz,axis = 2)

    psi_in_one = numpy.ones([nqx,nqy,nqz],dtype = 'complex')
    psi_in_two = numpy.zeros([nqx,nqy,nqz],dtype = 'complex')
    psi_out_one = numpy.ones([nqx,nqy,nqz],dtype = 'complex')
    psi_out_two = numpy.zeros([nqx,nqy,nqz],dtype = 'complex')

    psi =[psi_in_one,psi_in_two,psi_out_one,psi_out_two]

    magDensity = omf.ConvertRho()
    print '******'
    print magDensity[magDensity!=0.0]
    print '******'
    form_solution = magForm(cell.unit,magDensity, cell.value_list[0],
                            cell.value_list[1], cell.value_list[2],Q.q_list[0],
                            Q.q_list[1],Q.q_list[2],psi,qx_refract, omf.mx,
                            omf.my, omf.mz, precision=precision)

    return form_solution

def cudaBA_form(cell,Q,lattice,beam,precision='float32', refract = True):
    '''
    Overview:
        Uses the cuda code to solve the Born Approximation for the scattering
    off of the given cell. This only returns the form factor. It uses the same
    infrastructure as the cudaSMBA_form but sets the solution to the
    wavefunction to 1 so that the wavefunction perturbation does not take place.


    Parameters:

    cell:(Unit_Cell) = A Unit_Cell object that holds the information needed to
    calculate the form factor.

    Q:(q_space) = A Q_space object that holds all of the information about the
    desired q space output.

    lattice:(Lattice) = A lattice object that holds all of the information
    needed to solve the structure factor of the scattering.

    beam:(Beam) = Holds all of the information about the experimental beam
    needed to apply beam dependent corrections to the data.

    precision:(str|precision) = Allows the user to specify how precision the
    calculation will be. The choices here are generally float32 and float64.


    Note:
    - The precision at float32 precision seems to be accurate enough for
    calculations and can save a significant amount of calculation time.

    '''

    nqx = numpy.shape(Q.q_list[0])[0]
    nqy = numpy.shape(Q.q_list[1])[0]
    nqz = numpy.shape(Q.q_list[2])[0]

    if refract == True:
        stack = cell.inc_sub
        psi_in_one,psi_in_two,psi_out_one,psi_out_two,qx_refract = (
                    smba_wave_driver.wave(stack, Q.q_list[0], Q.q_list[1],
                             Q.q_list[2],beam.wavelength,cell.step[2],
                             precision=precision))
        Q.qx_refract = qx_refract
    else:
        qx_refract = Q.q_list[0]
        qx_refract = numpy.reshape(qx_refract,[numpy.size(qx_refract),1,1])
        qx_refract = numpy.repeat(qx_refract,nqy,axis = 1)
        qx_refract = numpy.repeat(qx_refract,nqz,axis = 2)

    psi_in_one = numpy.ones([nqx,nqy,nqz],dtype = 'complex')
    psi_in_two = numpy.zeros([nqx,nqy,nqz],dtype = 'complex')
    psi_out_one = numpy.ones([nqx,nqy,nqz],dtype = 'complex')
    psi_out_two = numpy.zeros([nqx,nqy,nqz],dtype = 'complex')

    psi =[psi_in_one,psi_in_two,psi_out_one,psi_out_two]

    form_solution = form(cell.unit,cell.value_list[0],cell.value_list[1],
                    cell.value_list[2],Q.q_list[0],Q.q_list[1],Q.q_list[2],
                    psi,qx_refract,precision=precision)

    return form_solution


def cudaSMBA_form(cell,Q,lattice,beam,precision ='float32', refract = True):
    '''
    Overview:
        Uses the cuda code to solve the Substrate Modified Born Approximation
    for the scattering off of the given cell. This only returns the form factor.
    it perturbed the scattering by the solution to the wavefunction for the
    scattering between the incident and substrate media. It produces notable
    features such as horizon effects.


    Parameters:

    cell:(Unit_Cell) = A Unit_Cell object that holds the information needed to
    calculate the form factor.

    Q:(q_space) = A Q_space object that holds all of the information about the
    desired q space output.

    lattice:(Lattice) = A lattice object that holds all of the information
    needed to solve the structure factor of the scattering.

    beam:(Beam) = Holds all of the information about the experimental beam
    needed to apply beam dependent corrections to the data.

    precision:(str|precision) = Allows the user to specify how precision the
    calculation will be. The choices here are generally float32 and float64.

    '''

    nqx = numpy.shape(Q.q_list[0])[0]
    nqy = numpy.shape(Q.q_list[1])[0]
    nqz = numpy.shape(Q.q_list[2])[0]

    if refract == True:
        stack = cell.inc_sub
        psi_in_one,psi_in_two,psi_out_one,psi_out_two,qx_refract = (
                    smba_wave_driver.wave(stack, Q.q_list[0], Q.q_list[1],
                             Q.q_list[2],beam.wavelength,cell.step[2],
                             precision=precision))
        Q.qx_refract = qx_refract
    else:
        qx_refract = Q.q_list[0]
        qx_refract = numpy.reshape(qx_refract,[numpy.size(qx_refract),1,1])
        qx_refract = numpy.repeat(qx_refract,nqy,axis = 1)
        qx_refract = numpy.repeat(qx_refract,nqz,axis = 2)

    #stack = approximations.wavefunction_format(cell.unit,Q.q_step[2])
    stack = cell.inc_sub

    psi_in_one,psi_in_two,psi_out_one,psi_out_two,qx_refract = (
                    smba_wave_driver.wave(stack, Q.q_list[0], Q.q_list[1],
                         Q.q_list[2],beam.wavelength,cell.step[2],
                         precision=precision))

    print 'WAVEFUNCTION CALCULATED'

    Q.qx_refract = qx_refract

    psi =[psi_in_one,psi_in_two,psi_out_one,psi_out_two]

    form_solution = form(cell.unit,cell.value_list[0],cell.value_list[1],
                        cell.value_list[2],Q.q_list[0],Q.q_list[1],Q.q_list[2],
                        psi,qx_refract,precision=precision)


    print 'FORM FACTOR CALCULATED'

    return form_solution



#---------------------------Form calculation----------------------------
def readfile(name):

    file = open(name)
    txt = file.read()
    file.close()
    return txt

def loadformkernelsrc(name, precision='float32', defines={}):
    '''
    Overview:
        Loads the kernel that will be used to solve the form factor for the
    scattering.
    '''
    import os
    src = readfile(os.path.join(os.path.dirname(__file__),name))

    defines = "\n".join(('#define %s %s'%(k,str(defines[k])))
                        for k in sorted(defines.keys()))

    if precision == 'float32':
        typedefs = '''
#define HAVE_CUDA
#include <lib/cufloat.h>

        '''
    else:
        typedefs = '''
#define HAVE_CUDA
#include <lib/cudouble.h>
        '''
    src = defines+typedefs+src

    return SourceModule(src, no_extern_c=True, include_dirs=[os.path.abspath
                                                 (os.path.dirname(__file__))])


def magForm(density, magDensity, x, y, z, Qx, Qy, Qz, psi, qx_refract,mx,my,mz,
            gpu=None, precision='float32'):

    '''
    Overview:
        Sets up the parallelized calculation of the form factor for calculation
    on an nvidia graphics card over all devices. This adds some input and output
    variables on top of the form variables for magnetic calculations.


    Parameters:

    density:(3D array|angstroms^2) = The scattering length density array being
    scattered off of.

    x:(array|angstroms) = An array of real space values that represent the
    location of the discretized units of density in the x direction.

    y:(array|angstroms) = An array of real space values that represent the
    location of the discretized units of density in the y direction.

    z:(array|angstroms) = An array of real space values that represent the
    location of the discretized units of density in the z direction.

    Qx:(array|angstroms^-1) = The array of Qx values that are being solved for.
    This array determines the output.

    Qy:(array|angstroms^-1) = The array of Qy values that are being solved for.
    Although solved for, this is the axis that is integrated over. It simulates
    the contribution of Qy scattering to the data.

    Qz:(array|angstroms^-1) = The array of Qz values that are being solved for.
    This array determines the output.

    gpu:(int) = The number of devices the user wants to utilize when solving
    the wavefunction. This defaults the the maximum number of available
    resources.

    precision:(str|precision) = Allows the user to specify how precision the
    calculation will be. The choices here are generally float32 and float64.

    '''
    density,magDensity,x,y,z,Qx,Qy,Qz,mx,my,mz, qx_refract = [
                                  numpy.asarray(v,precision) for v in
                              density,magDensity,x,y,z,Qx,Qy,Qz,mx,my,mz,
                              qx_refract]

    cplx = 'complex64' if precision=='float32' else 'complex128'

    psi= ([numpy.asarray(v,cplx) for v in psi])


    if gpu is not None:
        numgpus = 1
        gpus = [gpu]

    else:
        numgpus = cuda.Device.count()
        gpus = range(numgpus)

    size = [len(v) for v in Qx,Qy,Qz]

    uuresult = numpy.ones(size,dtype=cplx)
    ddresult = numpy.ones(size,dtype=cplx)
    udresult = numpy.ones(size,dtype=cplx)
    duresult = numpy.ones(size,dtype=cplx)

    work_queue = Queue.Queue()

    for qxi,qx in enumerate(Qx): work_queue.put(qxi)

    threads = [MagFormThread(gpus[k], work_queue, uuresult,ddresult,udresult,duresult,
                      density, magDensity, x, y, z, mx, my, mz, Qx, Qy, Qz,
                      psi,qx_refract)
           for k in range(numgpus)]

    for T in threads: T.start()
    for T in threads: T.join()

    return [uuresult, ddresult, udresult, duresult]

class MagFormThread(threading.Thread):
    '''
    Overview:
        Threads the kernel and the parameters to the GPU card device for a
    parallelized calculation.

    Note:
    -The threading in this object occurs over a single plain. Within the design
    of this use, it threads over each Qx plane (Qy x Qz) and then moves on to
    the next plane. Threading over the full 3D array proved challenging because
    of limited resources on the card.
    '''
    def __init__(self, gpu, work_queue, uuresult,ddresult,udresult,duresult,
                 density, magDensity, x, y, z, mx, my, mz, Qx, Qy, Qz,
                 psi, qx_refract):

        threading.Thread.__init__(self)

        self.born_args = density, magDensity, x, y, z,mx,my,mz, Qx, Qy, Qz, uuresult,ddresult,udresult,duresult, psi, qx_refract
        self.work_queue = work_queue
        self.gpu = gpu
        self.precision = x.dtype

    def run(self):
        '''
        Overview:
            The 'executed' part of the thread object.
        '''
        self.dev = cuda.Device(self.gpu)
        self.ctx = self.dev.make_context()
        self.cudamod = loadformkernelsrc("lib/mag_smba_kernel.c",
                                         precision=self.precision)

        self.cudaBorn = self.cudamod.get_function("cudaBorn")
        self.magFormKernel()
        self.ctx.pop()
        del self.ctx
        del self.dev


    def magFormKernel(self):
        '''
        Overview:
            The kernel that is loaded by run to the device which calculates the
        form factor for a given set of Q space parameters.

        '''
        density, magDensity, x, y, z, mx, my, mz = self.born_args[:8]
        Qx, Qy, Qz, uuresult,ddresult,udresult,duresult = self.born_args[8:15]
        psi, qx_refract = self.born_args[15:]

        nx, ny, nz, nqx, nqy, nqz = [numpy.int32(len(v))
                                     for v in x, y, z, Qx, Qy, Qz]

        fdensity = density.flatten()
        fMagDensity = magDensity.flatten()

        fmx = mx.flatten()
        fmy = my.flatten()
        fmz = mz.flatten()

        cx,cy,cz,cmx,cmy,cmz,cQx,cQy,cQz,cdensity,cmagDensity = [gpuarray.to_gpu(v)
                                          for v in x, y, z, fmx, fmy, fmz, Qx, Qy,
                                          Qz,fdensity,fMagDensity]

        cuuframe = cuda.mem_alloc((uuresult[0].nbytes))
        cddframe = cuda.mem_alloc((ddresult[0].nbytes))
        cudframe = cuda.mem_alloc((udresult[0].nbytes))
        cduframe = cuda.mem_alloc((duresult[0].nbytes))

        n = int(1*nqy*nqz)

        while True:
            try:
                qxi = numpy.int32(self.work_queue.get(block=False))
            except Queue.Empty:
                break

            print "%d of %d on %d\n"%(qxi,nqx,self.gpu),

            cpsi_in_one = gpuarray.to_gpu((psi[0][qxi,:,:]).flatten())
            cpsi_in_two = gpuarray.to_gpu((psi[1][qxi,:,:]).flatten())
            cpsi_out_one = gpuarray.to_gpu((psi[2][qxi,:,:]).flatten())
            cpsi_out_two = gpuarray.to_gpu((psi[3][qxi,:,:]).flatten())

            cqx_refract = gpuarray.to_gpu(qx_refract[qxi,:,:].flatten())

            self.cudaBorn(nx,ny,nz,nqx,nqy,nqz,
                     cdensity,cmagDensity, cmx, cmy, cmz, cx, cy, cz, cQx, cQy, cQz,
                     cpsi_in_one,cpsi_in_two,
                     cpsi_out_one,cpsi_out_two,
                     cqx_refract,
                     qxi,cuuframe,cddframe,cudframe,cduframe,
                     **cuda_partition(n))

            ## Delay fetching result until the kernel is complete
            cuda_sync()

            ## Fetch result back to the CPU
            cuda.memcpy_dtoh(uuresult[qxi], cuuframe)
            cuda.memcpy_dtoh(ddresult[qxi], cddframe)
            cuda.memcpy_dtoh(udresult[qxi], cudframe)
            cuda.memcpy_dtoh(duresult[qxi], cduframe)

        del cx, cy, cz, cQx, cQy, cQz, cdensity, cuuframe,cddframe,cudframe
        del cduframe,cqx_refract, cmx, cmy, cmz
        del cpsi_in_one,cpsi_out_one,cpsi_in_two,cpsi_out_two

def form(density, x, y, z, Qx, Qy, Qz, psi, qx_refract, gpu=None,
         precision='float32'):
    '''
    Overview:
        Sets up the parallelized calculation of the form factor for calculation
    on an nvidia graphics card over all devices.


    Parameters:

    density:(3D array|angstroms^2) = The scattering length density array being
    scattered off of.

    x:(array|angstroms) = An array of real space values that represent the
    location of the discretized units of density in the x direction.

    y:(array|angstroms) = An array of real space values that represent the
    location of the discretized units of density in the y direction.

    z:(array|angstroms) = An array of real space values that represent the
    location of the discretized units of density in the z direction.

    Qx:(array|angstroms^-1) = The array of Qx values that are being solved for.
    This array determines the output.

    Qy:(array|angstroms^-1) = The array of Qy values that are being solved for.
    Although solved for, this is the axis that is integrated over. It simulates
    the contribution of Qy scattering to the data.

    Qz:(array|angstroms^-1) = The array of Qz values that are being solved for.
    This array determines the output.

    gpu:(int) = The number of devices the user wants to utilize when solving
    the wavefunction. This defaults the the maximum number of available
    resources.

    precision:(str|precision) = Allows the user to specify how precision the
    calculation will be. The choices here are generally float32 and float64.

    '''

    density,x,y,z,Qx,Qy,Qz, qx_refract = [numpy.asarray(v,precision) for v in
                              density,x,y,z,Qx,Qy,Qz,qx_refract]

    cplx = 'complex64' if precision=='float32' else 'complex128'

    psi= ([numpy.asarray(v,cplx) for v in psi])


    if gpu is not None:
        numgpus = 1
        gpus = [gpu]

    else:
        numgpus = cuda.Device.count()
        gpus = range(numgpus)

    size = [len(v) for v in Qx,Qy,Qz]

    result = numpy.empty(size,dtype=cplx)

    work_queue = Queue.Queue()

    for qxi,qx in enumerate(Qx): work_queue.put(qxi)

    threads = [FormThread(gpus[k], work_queue, result,
                      density, x, y, z,Qx, Qy, Qz,
                      psi,qx_refract)
           for k in range(numgpus)]

    for T in threads: T.start()
    for T in threads: T.join()
    return result

class FormThread(threading.Thread):
    '''
    Overview:
        Threads the kernel and the parameters to the GPU card device for a
    parallelized calculation.

    Note:
    -The threading in this object occurs over a single plain. Within the design
    of this use, it threads over each Qx plane (Qy x Qz) and then moves on to
    the next plane. Threading over the full 3D array proved challenging because
    of limited resources on the card.
    '''
    def __init__(self, gpu, work_queue, result, density, x, y, z, Qx, Qy, Qz,
                                                            psi, qx_refract):

        threading.Thread.__init__(self)
        self.born_args = density, x, y, z, Qx, Qy, Qz, result, psi, qx_refract
        self.work_queue = work_queue
        self.gpu = gpu
        self.precision = x.dtype

    def run(self):
        '''
        Overview:
            The 'executed' part of the thread object.
        '''
        self.dev = cuda.Device(self.gpu)
        self.ctx = self.dev.make_context()
        self.cudamod = loadformkernelsrc("lib/smba_kernel.c",
                                         precision=self.precision)

        self.cudaBorn = self.cudamod.get_function("cudaBorn")
        self.formkernel()
        self.ctx.pop()
        del self.ctx
        del self.dev


    def formkernel(self):
        '''
        Overview:
            The kernel that is loaded by run to the device which calculates the
        form factor for a given set of Q space parameters.

        '''
        density, x, y, z, Qx, Qy, Qz, result, psi, qx_refract = self.born_args

        nx, ny, nz, nqx, nqy, nqz = [numpy.int32(len(v))
                                     for v in x, y, z, Qx, Qy, Qz]

        fdensity = density.flatten()

        cx,cy,cz,cQx,cQy,cQz,cdensity = [gpuarray.to_gpu(v)
                                          for v in x, y, z, Qx, Qy, Qz,fdensity]

        cframe = cuda.mem_alloc(result[0].nbytes)

        n = int(1*nqy*nqz)

        while True:
            try:
                qxi = numpy.int32(self.work_queue.get(block=False))
            except Queue.Empty:
                break

            print "%d of %d on %d\n"%(qxi,nqx,self.gpu),

            cpsi_in_one = gpuarray.to_gpu((psi[0][qxi,:,:]).flatten())
            cpsi_in_two = gpuarray.to_gpu((psi[1][qxi,:,:]).flatten())
            cpsi_out_one = gpuarray.to_gpu((psi[2][qxi,:,:]).flatten())
            cpsi_out_two = gpuarray.to_gpu((psi[3][qxi,:,:]).flatten())

            cqx_refract = gpuarray.to_gpu(qx_refract[qxi,:,:].flatten())

            self.cudaBorn(nx,ny,nz,nqx,nqy,nqz,
                     cdensity,cx,cy,cz,cQx,cQy,cQz,
                     cpsi_in_one,cpsi_in_two,
                     cpsi_out_one,cpsi_out_two,
                     cqx_refract,
                     qxi,cframe,
                     **cuda_partition(n))

            ## Delay fetching result until the kernel is complete
            cuda_sync()

            ## Fetch result back to the CPU
            cuda.memcpy_dtoh(result[qxi], cframe)

        del cx, cy, cz, cQx, cQy, cQz, cdensity, cframe,cqx_refract
        del cpsi_in_one,cpsi_out_one,cpsi_in_two,cpsi_out_two


def cuda_sync():
    """
    Overview:
        Waits for operation in the current context to complete.
    """
    #return # The following works in C++; don't know what pycuda is doing
    # Create an event with which to synchronize
    done = cuda.Event()

    # Schedule an event trigger on the GPU.
    done.record()

    #line added to not hog resources
    while not done.query(): time.sleep(0.01)

    # Block until the GPU executes the kernel.
    done.synchronize()
    # Clean up the event; I don't think they can be reused.
    del done

def cuda_partition(n):
    '''
    Overview:
        Auto grids the thread blocks to achieve some level of calculation
    efficiency.
    '''
    max_gx,max_gy = 65535,65535
    blocksize = 32

    block = (blocksize,1,1)
    num_blocks = int((n+blocksize-1)/blocksize)
    if num_blocks < max_gx:
        grid = (num_blocks,1)
    else:
        gx = max_gx
        gy = (num_blocks + max_gx - 1) / max_gx
        if gy >= max_gy: raise ValueError("vector is too large")
        grid = (gx,gy)

    return dict(block=block,grid=grid)

def main():
    print 'test here'
    return


if __name__ == "__main__":
    main()
