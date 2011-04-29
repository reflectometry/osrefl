# This program is in the public domain
# Authors: Paul Kienzle, Christopher Metting
#03/23/2010

import Queue
import threading
from pycuda import gpuarray
import pycuda.driver as cuda
from pycuda.compiler import SourceModule
import numpy
import numpy.linalg as linalg
import approximations,sample_prep
import time

cuda.init()

def readfile(name):
    file = open(name)
    txt = file.read()
    file.close()
    return txt

def loadkernelsrc(name, precision='float32', defines={}):
    import os
    src = readfile(os.path.join(os.path.dirname(__file__),name))
    # The following are currently defined by cufloat.h/cudouble.h, so aren't
    # needed here.
    #defines['CUDA_KERNEL'] = 'extern "C" __global__ void'
    #defines['INLINE'] = '__inline__ __host__ __device__'
    defines = "\n".join(('#define %s %s'%(k,str(defines[k]))) 
                        for k in sorted(defines.keys()))
    #define sin __sinf
    #define cos __cosf   
    if precision == 'float32':
        typedefs = '''
#define HAVE_CUDA
#include <cufloat.h>

#define sincos __sincosf
        '''
    else:
        typedefs = '''
#define HAVE_CUDA
#include <cudouble.h>
        '''
    src = defines+typedefs+src
 
    return SourceModule(src, no_extern_c=True, 
                    include_dirs=[os.path.abspath(os.path.dirname(__file__))])



def wave(stack, Qx, Qy, Qz,wavelength, deltaz, gpu=None, precision='float32', proc = 'gpu'):
    
    '''
    Overview:
        Calculates the wavefunction for a given stack of SLD values.
    It solves the specular reflectivity and can be used to determine the
    reflection and transmission from a given 1D stack.
    
    
    Parameters:
    
    stack([3,n]|[angstroms^2, angstroms, absorption?]) = The 1D stack to be
    scattered off of. This is in the form [SLD,thickness,absorption] and there
    is one element for every layer being scattered off of.
    
    Qx:(array|angstroms^-1) = The array of Qx values that are being solved for.
    This array determines the output for the wavefunctions calculated.
    
    Qy:(array|angstroms^-1) = The array of Qy values that are being solved for.
    This array determines the output for the wavefunctions calculated.
    
    Qz:(array|angstroms^-1) = The array of Qz values that are being solved for.
    This array determines the output for the wavefunctions calculated.
    
    wavelength:(array|angstroms) = The wavelength of the probing beam.
    
    gpu:(int) = The number of devices the user wants to utilize when solving
    the wavefunction. This defaults the the maximum number of available 
    resources.
    
    precision:(str|precision) = Allows the user to specify how precision the
    calculation will be. The choices here are generally float32 and float64.
    '''
    
    #Make sure there is a Cuda Device and, if not, use the python calculation
    


    stack, Qx, Qy, Qz = [numpy.asarray(v,precision) for v in 
                              stack, Qx, Qy, Qz]

    
    if precision == 'float32': 
        wavelength = numpy.float32(wavelength)
        deltaz = numpy.float32(deltaz)
    else:
        wavelength = numpy.float64(wavelength)
        deltaz = numpy.float64(deltaz)
        
    cplx = 'complex64' if precision=='float32' else 'complex128'
    
    if gpu is not None:
        numgpus = 1
        gpus = [gpu]
    else:
        numgpus = cuda.Device.count()
        gpus = range(numgpus)
        
    size = [len(v) for v in Qx,Qy,Qz]
    
    psi_in_one = numpy.zeros(size,dtype=cplx)
    psi_in_two = numpy.zeros(size,dtype=cplx)
    psi_out_one = numpy.zeros(size,dtype=cplx)
    psi_out_two = numpy.zeros(size,dtype=cplx)
    
    qx_refract =  numpy.zeros(size,dtype=precision)
    
    if cuda.Device > 0 and proc == 'gpu':
        
        print 'Cuda Device Detected... Calculating wavefunction on GPU'
    
        work_queue = Queue.Queue()
        for qxi,qx in enumerate(Qx): work_queue.put(qxi)
    
        threads = [WaveThread(gpus[k], work_queue, stack, wavelength,deltaz,
                    psi_in_one,psi_in_two,psi_out_one,psi_out_two,qx_refract,
                          Qx, Qy, Qz)
        
               for k in range(numgpus)]
        
        for T in threads: T.start()
        for T in threads: T.join()
        
    else:
        if cuda.Device == 0: print 'No Cuda Device Detected...'
        print 'Calculating wavefunction with Python'
        
        if proc == 'cpu': print 'CPU Calculation Selected'
        
        from approximations import SMBA_wavecalcMultiK,QxQyQz_to_k
        from numpy import array,asarray,reshape
        
        for q in Qx,Qy,Qz: q = asarray(q)
        
        Qx.reshape([len(Qx),1,1])
        
        Qx = array(Qx.reshape(len(Qx),1,1), dtype = precision)
        Qy = array(Qy.reshape(1,len(Qy),1), dtype = precision)
        Qz = array(Qz.reshape(1,1,len(Qz)), dtype = precision)
        
        kin,kout = QxQyQz_to_k(Qx,Qy,Qz,wavelength)
        
        psi_in_one,psi_in_two,psi_out_one,psi_out_two,qx_refract = (
                                   SMBA_wavecalcMultiK(qx_refract,deltaz, 
                                               stack, wavelength,kin, kout))
        
        
        
    return  psi_in_one,psi_in_two,psi_out_one,psi_out_two,qx_refract

class WaveThread(threading.Thread):
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
    
    def __init__(self, gpu, work_queue, stack, wavelength,deltaz,
                      psi_in_one,psi_in_two,psi_out_one,psi_out_two,qx_refract,
                      Qx, Qy, Qz):
        
        threading.Thread.__init__(self)
        self.wave_args = (Qx, Qy, Qz, stack, wavelength,deltaz,psi_in_one,
                          psi_in_two,psi_out_one,psi_out_two,qx_refract)
        self.work_queue = work_queue
        self.gpu = gpu
        self.precision = Qx.dtype
        
    def run(self):

        self.dev = cuda.Device(self.gpu)
        self.ctx = self.dev.make_context()
        self.cudamod = loadkernelsrc("wavefunction_kernel.cc",
                                  precision=self.precision)
        self.cudaWave = self.cudamod.get_function("cudaWave")
        self.wavekernel()
        self.ctx.pop()
        del self.ctx
        del self.dev


    def wavekernel(self):
        '''
        Overview:
            The kernel that is loaded by run to the device which calculates the
        wavefunction for a given set of Q space parameters.

        '''
        
        (Qx, Qy, Qz, stack, wavelength,deltaz,psi_in_one,psi_in_two,psi_out_one,
         psi_out_two,qx_refract) = self.wave_args
        
        nqx, nqy, nqz, nlayer = [numpy.int32(len(v)) for v in (Qx, Qy, 
                                                               Qz, stack[:,0])]

        SLD = numpy.copy(stack[:,0])

        thickness = numpy.copy(stack[:,1])
        mu = numpy.copy(stack[:,2])
        

        cSLD = gpuarray.to_gpu(SLD)
        cthickness = gpuarray.to_gpu(thickness)
        cmu = gpuarray.to_gpu(mu)

        cQx, cQy,cQz = [gpuarray.to_gpu(v) 
                                          for v in Qx, Qy, Qz]
        
        
        cpio = cuda.mem_alloc(psi_in_one[0].nbytes)
        cpit = cuda.mem_alloc(psi_in_two[0].nbytes)
        cpoo = cuda.mem_alloc(psi_out_one[0].nbytes)
        cpot = cuda.mem_alloc(psi_out_two[0].nbytes)

        cqxr = cuda.mem_alloc(qx_refract[0].nbytes)
        
        n = int(1*nqy*nqz)
        
        while True:
            try:
                qxi = numpy.int32(self.work_queue.get(block=False))
            except Queue.Empty:
                break

            print "%d of %d on %d\n"%(qxi,nqx,self.gpu),
            
            
            self.cudaWave(nqx, nqy, nqz, cQx, cQy, cQz, cSLD, cthickness, cmu, 
                          nlayer, wavelength,deltaz, qxi, cpio, cpit, cpoo, 
                          cpot, cqxr, **cuda_partition(n))

            
            ## Delay fetching result until the kernel is complete
            cuda_sync()

            ## Fetch result back to the CPU
            cuda.memcpy_dtoh(psi_in_one[qxi], cpio)

            cuda.memcpy_dtoh(psi_in_two[qxi], cpit)

            cuda.memcpy_dtoh(psi_out_one[qxi], cpoo)

            cuda.memcpy_dtoh(psi_out_two[qxi], cpot)

            cuda.memcpy_dtoh(qx_refract[qxi], cqxr)
            


        del cQx, cQy, cQz, cSLD,cthickness, cmu, cpio,cpit,cpoo ,cpot,cqxr




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
    #max_gx,max_gy = 5,65536
    #blocksize = 3
    block = (blocksize,1,1)
    num_blocks = int((n+blocksize-1)/blocksize)
    if num_blocks < max_gx:
        grid = (num_blocks,1)
    else:
        gx = max_gx
        gy = (num_blocks + max_gx - 1) / max_gx
        if gy >= max_gy: raise ValueError("vector is too large")
        grid = (gx,gy)
    #print "block",block,"grid",grid
    #print "waste",block[0]*block[1]*block[2]*grid[0]*grid[1] - n
    return dict(block=block,grid=grid)


def main():
    
    '''
    This is a test to compare the results of the wave function calculation
    between the Cuda and Python calculations.
    '''
    import wavefunction_BBM, approximations,pylab,numpy
    from pylab import subplot
    '''
    unit_size = (50,50,50)
    unit_metric = (1.0e5,1.0e5,2.5e3)
    feature_size = (25,25,41)
    #Q_size = (100,100,100)
    '''
    Q_size = (10,10,10)
    wavelength = 5.0
    '''
    unit = numpy.zeros(unit_size, 'd')
    unit[0:feature_size[0],0:feature_size[1],0:feature_size[2]] = 4.5e-6
    '''
    img = sample_prep.GrayImgUnit(filename =
          '/home/mettingc/Downloads/sample1_sld.png', 
          newres = numpy.array([200,400]))
    
    unit = img.unitBuild(Dxyz = [8480.0,8480.0,3500.0], 
                         scale = 1.0e-5,inc_sub=[0.0,2.0784314e-6])
    from pylab import imshow,show,colorbar
    from numpy import rot90
    unit.add_media()
 
    imshow(rot90(unit.unit[:,5,:]))

    show()
    unit.viewSlice()

    gpu=None
    #gpu=1  # Use this line if you only want one GPU to be in use

    qx = numpy.linspace(-0.05,0.075,Q_size[0])
    qy = numpy.linspace(-0.0025,0.001,Q_size[1])
    qz = numpy.linspace(0.00002,0.06,Q_size[2])
    
    extent = (-0.05,0.075,0.00002,0.06)
    x = unit.value_list[0]
    y = unit.value_list[1]
    z = unit.value_list[2]
    
    #stack = approximations.wavefunction_format(unit.unit, unit.step[2])
    stack = unit.inc_sub
    '''
    x = (unit_metric[0]/unit_size[0])*numpy.arange(unit_size[0])
    y = (unit_metric[1]/unit_size[1])*numpy.arange(unit_size[1])
    z = (unit_metric[2]/unit_size[2])*numpy.arange(unit_size[2])

    avg = numpy.average(unit, axis = 0)
    avg = numpy.average(avg, axis = 0)
    thickness = numpy.array(numpy.ones(unit_size[2])*(z[1]-z[0]))
    mu = numpy.array(numpy.zeros(unit_size[2]))
    
    stack = numpy.zeros([unit_size[2],3])
    for i in range(unit_size[2]):
        stack[i] = avg[i],thickness[i],mu[i]
    stack = numpy.asarray(stack)
    '''
    t0 = time.time()
    psi_in_one,psi_in_two,psi_out_one,psi_out_two,qx_refract = wave(stack, 
                            qx, qy, qz,wavelength, gpu=gpu, precision='float64')
    print "time",time.time()-t0
    
    
    
    qxv = numpy.reshape(qx,[Q_size[0],1,1])
    qyv = numpy.reshape(qy,[1,Q_size[1],1])
    qzv = numpy.reshape(qz,[1,1,Q_size[2]])
    r = numpy.zeros(Q_size,'complex')
    t0 = time.time()

    ki,kf = approximations.QxQyQz_to_k(qxv, qyv, qzv, wavelength)
    
    
    pio = numpy.zeros(ki.shape,'complex')
    pit =numpy.zeros(ki.shape,'complex')
    poo = numpy.zeros(ki.shape,'complex')
    pot =numpy.zeros(ki.shape,'complex')
    pyrefract =numpy.zeros(ki.shape,'complex')
    t0 = time.time()
    for i in range (Q_size[0]):
        print i
        for ii in range (Q_size[1]):
            for iii in range(Q_size[2]):

                
                pio[i,ii,iii],pit[i,ii,iii],poo[i,ii,iii],pot[i,ii,iii],pyrefract[i,ii,iii] = approximations.SMBA_wavecalc(qx[i],(qz[1]-qz[0]),stack,wavelength, ki[i,ii,iii],kf[i,ii,iii])
    
    print "time",time.time()-t0
    
    print 'done'

    BBM =[pio,pit,poo,pot]
    cuda = [psi_in_one,psi_in_two,psi_out_one,psi_out_two]
    
    '''
    pylab.imshow(numpy.rot90(numpy.sum(pio.real,axis = 1)))
    pylab.colorbar()
    pylab.figure()
    pylab.imshow(numpy.rot90(numpy.sum(psi_in_one.real,axis = 1)))
    pylab.colorbar()
    pylab.show()
    '''
    
    dif = [None,None,None,None]
    difIm = [None,None,None,None]
    intdif = [None,None,None,None]
    intdifIm = [None,None,None,None]
    
    for i in range(4):
        dif[i] = (numpy.abs(cuda[i].real - BBM[i].real)/(.5*(numpy.abs(cuda[i].real) + 
                                                   numpy.abs(BBM[i].real))))
        dif[i][(BBM[i].real==0.0) & (cuda[i].real==0.0)] = 0.0
        intdif[i] = numpy.log10(numpy.sum(dif[i],axis=1))
        #intdif[i] = (numpy.sum(dif[i],axis=1))

    for i in range(4):
        difIm[i] = (numpy.abs(cuda[i].imag - BBM[i].imag)/(.5*(numpy.abs(cuda[i].imag) + 
                                                   numpy.abs(BBM[i].imag))))
        difIm[i][(BBM[i].imag==0.0) & (cuda[i].imag==0.0)] = 0.0
        intdifIm[i] = numpy.log10(numpy.sum(difIm[i],axis=1))
        #intdif[i] = (numpy.sum(difIm[i],axis=1))
    
    print '--------psi values----------'
    print pit
    print '-----------------------------'
    print psi_in_two
    print '=================================='
    
    print '--------real diff values----------'
    print pit[dif[1] > 1e-17]
    print '------------------'
    print psi_in_two[dif[1] > 1e-17]
    print '=================================='
    print '--------imag diff values----------'
    print pit[difIm[1] > 1e-17]
    print '------------------'
    print psi_in_two[difIm[1] > 1e-17]
    print '___________________________________'
    print 'the maximum of the difference is: '
    print 'pio: ',numpy.max(dif[0])
    print 'pit: ',numpy.max(dif[1])
    print 'poo: ',numpy.max(dif[2])
    print 'pot: ',numpy.max(dif[3])
    print '___________________________________'
 
    pylab.colorbar()
       
    subplot(2,2,1)
    pylab.title('pio')
    pylab.imshow(numpy.flipud(intdif[0].T),extent = extent, aspect = 'auto')

    subplot(2,2,2)
    pylab.title('pit')
    pylab.imshow(numpy.flipud(intdif[1].T),extent = extent, aspect = 'auto')

    subplot(2,2,3)
    pylab.title('poo')
    pylab.imshow(numpy.flipud(intdif[2].T),extent = extent, aspect = 'auto')

    subplot(2,2,4)
    pylab.title('pot')
    pylab.imshow(numpy.flipud(intdif[3].T),extent = extent, aspect = 'auto')

    
    pylab.figure()
    pylab.colorbar()  
    subplot(2,2,1)
    pylab.title('pio')
    pylab.imshow(numpy.flipud(intdifIm[0].T),extent = extent, aspect = 'auto')

    subplot(2,2,2)
    pylab.title('pit')
    pylab.imshow(numpy.flipud(intdifIm[1].T),extent = extent, aspect = 'auto')

    subplot(2,2,3)
    pylab.title('poo')
    pylab.imshow(numpy.flipud(intdifIm[2].T),extent = extent, aspect = 'auto')

    subplot(2,2,4)
    pylab.title('pot')
    pylab.imshow(numpy.flipud(intdifIm[3].T),extent = extent, aspect = 'auto')

    pylab.show()
    
# Note: may want to try putting const arrays in texture memory rather
# than global memory; claims that card can do better memory caching
# hence maybe get better performance.
if __name__ == "__main__":
    main()
