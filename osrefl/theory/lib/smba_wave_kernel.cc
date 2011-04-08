//This Kernel solves the wavefunction
const Real PI(3.14159265358979323846264338327950288);

__device__ void rt_calc(Cplx k, const Real SLD[], const int nlayer,
		const Real thickness[], Cplx* t, Cplx* r, Cplx* kz_tran, int invert)
{

	const Cplx I(0.0,1.0);


	Cplx temp_v00(0.0);
	Cplx temp_v01(0.0);
	Cplx temp_v10(0.0);
	Cplx temp_v11(0.0);

	Cplx ml00(0.0);
	Cplx ml01(0.0);
	Cplx ml10(0.0);
	Cplx ml11(0.0);

	Cplx M00(1.0);
	Cplx M01(0.0);

	Cplx M10(0.0);
	Cplx M11(1.0);

	Cplx kz(0.0);

	Cplx knotz(0.0);

	Cplx n(0.0);
	Cplx no(0.0);
	Cplx nf(0.0);

	int iter = 0;
	int layer;

	if (invert == 1)
	{
		iter = -1;
		layer = nlayer-2;
		knotz = (sqrt(pow(k,2) + 4.0 * PI * SLD[nlayer-1]));

		no = (sqrt(1 - 4.0 * PI * SLD[nlayer-1]/pow(knotz,2)));
		nf = (sqrt(1 - 4.0 * PI * SLD[0]/pow(knotz,2)));
	}

	else
	{
		iter = 1;
		layer = 1;
		knotz = (sqrt(pow(k,2) + 4.0 * PI * SLD[0]));

		no = (sqrt(1 - 4.0 * PI * SLD[0]/pow(knotz,2)));
		nf = (sqrt(1 - 4.0 * PI * SLD[nlayer-1]/pow(knotz,2)));

	}


	for (int z=1; z < nlayer-1; z++) {


		n = sqrt(1.0 - (4.0 *PI* SLD[layer]/pow(knotz,2)));
		kz = n * k;

		ml00 = ml11 = 1*cos(kz*thickness[layer]);
		ml01 = (1/n) * sin(kz*thickness[layer]);
		ml10 = -n * sin(kz*thickness[layer]);

		temp_v00 = (ml00*M00) + (ml10*M01);
		temp_v10 = (ml00*M10) + (ml10*M11);

		temp_v01 = (ml01*M00) + (ml11*M01);
		temp_v11 = (ml01*M10) + (ml11*M11);

		M00 = temp_v00;
		M01 = temp_v01;
		M10 = temp_v10;
		M11 = temp_v11;

		layer += iter;
	}


	Cplx frac(0.0,0.0);
	Cplx num(0.0,0.0);
	Cplx denom(0.0,0.0);

	frac = (-1*I) / nf;

	num = M11 + (I * no * M01) + frac * (-M10 - I * no * M00);
	denom = -M11 + I*no*M01 + frac*(M10 - I * no * M00);

	*r = num/denom;

	if (nf.real() == 0.0) *t = 0.0;
	else *t = 1.0 + *r;

	*kz_tran = nf * knotz;

	return;
}

CUDA_KERNEL


cudaWave(const int nqx, const int nqy, const int nqz,
		const Real Qx[], const Real Qy[], const Real Qz[], const Real SLD[],
		const Real thickness[], const Real mu[],const int nlayer,
		const Real wavelength, const int qxi,
		Cplx pio[],Cplx pit[],Cplx poo[],Cplx pot[], Real qxr[])
		
{


		const int idx = (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x + threadIdx.x;
		if (idx >= nqy*nqz) return;

		const int qzi = idx%nqz;
		const int qyi = (idx/nqz)%nqy;
		const Cplx I(0.0,1.0);

		Cplx t(0,0);
		Cplx r(0,0);
		Cplx kz_tran(0,0);


		Real qzstep = Qz[1]-Qz[0];

		Real Qmag = sqrt(pow(Qx[qxi],2) + pow(Qy[qyi],2) + pow(Qz[qzi],2));

		if (Qmag == 0.0) Qmag = 1.0;

		Real knot = (2.0 * PI)/wavelength;

		Cplx twoth = 2.0 * asin(Qmag/(2.0*knot));

		Cplx tilt = atan2(Qx[qxi],Qz[qzi]);


		Cplx th_in = (twoth/2.0) + tilt;
		Cplx th_out = (twoth/2.0) - tilt;


		Cplx kz_in = knot * sin(th_in);
		Cplx kz_out = -knot * sin(th_out);



		if  (kz_in.real() < 0.0)
		{
			rt_calc(-kz_in, &SLD[0], nlayer, &thickness[0], &t, &r, &kz_tran, 1);
			pio[idx] = t*exp(-I*kz_tran*qzstep);
			pit[idx] = 0.0;
			qxr[idx] = Qx[qxi] + wavelength * SLD[nlayer-1];

		}

		else
		{
			rt_calc(kz_in,&SLD[0], nlayer, &thickness[0], &t, &r, &kz_tran, 0);
			pio[idx] = 1.0*exp(I*kz_in*qzstep);
			pit[idx] = r*exp(-I*kz_in*qzstep);
			qxr[idx] = Qx[qxi];

		}


		if (kz_out.real() < 0.0)
		{
			rt_calc(-kz_out, &SLD[0], nlayer, &thickness[0], &t, &r, &kz_tran, 0);
			poo[idx] = 1.0*exp(-I*kz_out*qzstep);
			pot[idx] = r * exp(I*kz_out*qzstep);

		}

		else
		{
			rt_calc(kz_out, &SLD[0], nlayer, &thickness[0], &t, &r, &kz_tran, 0);
			poo[idx] = t*exp(-I*kz_tran*qzstep);
			pot[idx] = 0.0;
			qxr[idx] = Qx[qxi] - wavelength * SLD[nlayer-1];

		}
}


