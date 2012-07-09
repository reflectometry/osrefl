CUDA_KERNEL

cudaDWBA_part1(const Real qx, const Real qy,
	 const Real xstep, const Real ystep,
	 const int nx, const int ny,
	 const Real xx[], const Real yy[],
	 const Real rtor[1000][1000][1000],
	 const Real Vfac,
	 const bool lattice,
	 Cplx output[])
{

	const int x = threadIdx.x + blockIdx.x * blockDim.x;
	const int y = threadIdx.y + blockIdx.y * gridDim.y;
	int offset = x + y * blockDim.x * gridDim.x;
	if (offset >= ny) return;
	
	Cplx laux;
	Cplx lauy;

	Cplx ftwRef;

	const Cplx I(0.0,1.0);

	if (qx != 0 )
		laux = ((-1 * I) / qx) * exp(I * qx * (xstep - 1));
	if (qx != 0 )
		lauy = ((-1 * I) / qy) * exp(I * qy * (ystep - 1));
	
	for(int i = 0; i < nx; ++i)
		for(int j = 0; j < ny; ++j)
			ftwRef += rtor[i][j][offset] * exp(I * qx * xx[offset]) * exp(I * qy * yy[offset]);


	output[offset] = ftwRef * Vfac * laux * lauy; 


}

