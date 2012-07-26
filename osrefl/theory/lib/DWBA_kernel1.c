// Set Constants

const Real pi(3.14159265358979323846264338327950288);
const Cplx I(0.0,1.0);
const Real m (1.674e-27);
const Real h_bar (6.62607e-14);
const Real Vfac (-m/(2 * pi * h_bar * h_bar));



//const Real x[MAX_DIM][MAX_DIM][MAX_DIM], 
//const Real y[MAX_DIM][MAX_DIM][MAX_DIM], 
//const Real z[MAX_DIM][MAX_DIM][MAX_DIM],

// Calculates the C and D values for the given SLD and k values
// returns a one dimensional array where [0] = C and [1] = D
__device__ Cplx[] 
DWBA_wave_function(Real SLDArray[][], Real k[][][], 
		   int i, int ii, int iii)
{
	Cplx wave_function[];

	Real SLDinc;
	Real SLDsub;
	Real kz;

	Real SLD;
	Real thickness;

	Cplx k0z;
	Cplx nzsub;
	Cplx nz0;
	Cplx nz;
	Cplx kzl;

	Cplx C1;
	Cplx C2;
	Cplx C3;
	Cplx C4;

	Cplx r;

	Cplx c;
	Cplx d;
	
	Cplx p;
	Cplx pp;

	int z_interface;

	SLDinc = SLDArray[0][0];
	SLDsub = SLDArray[zsize][0];
	kz = k[i][ii][iii];

	SLD = SLDArray[iii][0];
	thickness = SLDArray[iii][1];

	k0z = sqrt(kz * kz + 4 * pi * SLDinc);
	nzsub = sqrt(1 - 4 * pi * SLDsub);
	nz0 = sqrt(1 - 4 * pi * SLDinc / (k0z * k0z));
	nz = sqrt(1 - 4 * pi * SLD / (k0z * k0z));
	kzl = k0z * nz;

	C1 = cos(kzl * thickness);
	C2 = (1 / nz) * sin(kzl * thickness);
	C3 = -nz * sin(kzl * thickness);
	C4 = cos(kzl * thickness);

	r = (C1 + I * nz0 * C2 + 1 / (I * nzsub) * -C3 - I * nz0 * C4) /
	    (-C1 + I * nz0 * C2 + 1 / (I * nzsub) * C3 - I * nz0 * C4);

	if(iii == 0)
	{
		c = 1.0;
		d = r;
	} 
	else if(iii == zsize) 
	{
		if(nzsub.real() != 0.0)
			c = 1.0 + r;
		else
			c = 0.0;
		d = 0.0;		
	}
	else
	{
		p = 1.0 + r;
		pp = I * k[i][ii][0];
		z_interface = (iii - 1) * thickness;
		
		for(int j = 0; 1 < j && j < iii; ++i)
		{
			p = (p * C1) + (pp * C2 / k0z);
			pp = (p * C3 * k0z) + (pp * C4);
		}
		
		c = 0.5 * exp(-I * kzl * z_interface) * 
			 (p + pp / (I * kzl));
		d = 0.5 * exp(I * kzl * z_interface) * 
			 (p - pp / (I * kzl));
	}

	wave_function[0] = c;
	wave_function[1] = d;

	return wave_function;

}


CUDA_KERNEL

// Set maximum Q space x, y, and z value (or dimension) 

const int MAX_DIM = 1000;

// Cuda DWBA implementation that takes in various parameters 
// returns through variable scatOut
cudaDWBA(const Real RTOR[MAX_DIM][MAX_DIM][MAX_DIM],
	 const int csx, const int csy, const int csz,
	 const Real SLDArray[MAX_DIM][],
	 const Real kin[MAX_DIM][MAX_DIM][MAX_DIM], 
	 const Real kout[MAX_DIM][MAX_DIM][MAX_DIM],
	 const int xsize, const int ysize, const int zsize,
	 Cplx scatOut[MAX_DIM][MAX_DIM][MAX_DIM])
{
	// 3) Variable Declarations

	Cplx wave_function[];
	Cplx pio;
	Cplx pit;
	Cplx poo;
	Cplx pot;

	// 4) Calculate indices	

	const int i = blockDim.x * blockIdx.x + threadIdx.x;
	const int ii = blockDim.y * blockIdx.y + threadIdx.y;
	const int iii = blockDim.z * blockIdx.z + threadIdx.z;

	if(i >= xsize) return;
	if(ii >= ysize) return;
	if(iii >= zsize) return;

	// 5) calculate c & d with DWBA_wave_function

	wave_function = DWBA_wave_function(SLDArray, kin, i, ii, iii);
	pio = wave_function[0];
	pit = wave_function[1];

	wave_function = DWBA_wave_function(SLDArray, kout, i, ii, iii);
	poo = wave_function[0];
	pot = wave_function[1];

	// 6) rest of DWBA calculation

	

	// 7) set scatOut and memcopy from host

}



