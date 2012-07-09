CUDA_KERNEL

cudaDWBA_part1(const Real qx, const Real qy,
	 const Real xstep, const Real ystep,
	 const int nx, const int ny,
	 const Real xx[], const Real yy[],
	 const Real rtor[][][],
	 const Real Vfac,
	 const bool lattice,
	 Cplx output[]
)
{

	const int x = threadIdx.x + blockIdx.x * blockDim.x;
	const int y = threadIdx.y + blockIdx.y * gridDim.y;
	int offset = x + y * blockDim.x * gridDim.x;
	if (offset >= nx * ny) return;
	
	Cplx laux;
	Cplx lauy;

	Cplx ftwRef;

	const Cplx I(0.0,1.0);

	if (qx[offset] != 0 )
		laux = ((-1I / qx[offset]) * exp(1I * qx[offset] * xstep[offset] - 1.0));
	if (qx[offset] != 0 )
		lauy = ((-1I / qy[offset]) * exp(1I * qy[offset] * ystep[offset] - 1.0));
	
	for(int i = 0; i < nx; ++i)
		for(int j = 0; j < ny; ++j)
			ftwRef += rtor[i][j][offset] * exp(1I * qx[offset] * xx[offset]) * exp(1I * qy[offset] * yy[offset]);


	output[offset] = ftwRef * Vfac * laux * lauy; 


}

