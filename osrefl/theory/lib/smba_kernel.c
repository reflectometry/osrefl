CUDA_KERNEL

cudaBorn(int nx, int ny, int nz, int nqx, int nqy, int nqz,
		const Real density[],
		const Real x[], const Real y[], const Real z[],
		const Real Qx[], const Real Qy[], const Real Qz[],
		const Cplx psi_in_one[],
		const Cplx psi_in_two[],
		const Cplx psi_out_one[],
		const Cplx psi_out_two[],
		const Real qx_refract[],
		int qxi,
		Cplx result[])
{

    const int idx = (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x + threadIdx.x;
    if (idx >= nqy*nqz) return;

    const int qzi = idx%nqz;
    const int qyi = (idx/nqz)%nqy;



    const Cplx I(0.0,1.0);

    Cplx scat(0.0,0.0);
    Cplx ft(0.0,0.0);

    Cplx scat_PioPoo;
    Cplx scat_PioPot;
    Cplx scat_PitPoo;
    Cplx scat_PitPot;


    int densityidx = 0;

    for (int xi=0; xi < nx; xi++) {
        for (int yi=0; yi < ny; yi++) {
            for (int zi=0; zi < nz; zi++) {

				//ft = density[densityidx]*exp(I*Qx[qxi]*x[xi])*exp(I*Qy[qyi]*y[yi])*exp(I*Qz[qzi]*z[zi]);
				ft = density[densityidx]*exp(I*qx_refract[idx]*x[xi])*exp(I*Qy[qyi]*y[yi])*exp(I*Qz[qzi]*z[zi]);

                scat_PioPoo = psi_in_one[idx] * ft * psi_out_one[idx];
                scat_PioPot = psi_in_one[idx] * ft * psi_out_two[idx];
                scat_PitPoo = psi_in_two[idx] * ft * psi_out_one[idx];
                scat_PitPot = psi_in_two[idx] * ft * psi_out_two[idx];

                scat += (scat_PioPoo + scat_PioPot + scat_PitPoo + scat_PitPot);

                densityidx++;
            }
       }
    }

	if (Qx[qxi] == 0.0)
	{
		scat *= (x[1]-x[0]);

	}
	else
	{
		scat *= (-I/Qx[qxi])*(1-exp(I*Qx[qxi]*(x[1]-x[0])));

	}
	if (Qy[qyi] == 0.0)
	{
		scat *= (y[1]-y[0]);

	}
	else
	{
		scat *= (-I/Qy[qyi])*(1-exp(I*Qy[qyi]*(y[1]-y[0])));


	}
	if (Qz[qzi] == 0.0)
	{
		scat *= (z[1]-z[0]);

	}
	else
	{
		scat *= (-I/Qz[qzi])*(1-exp(I*Qz[qzi]*(z[1]-z[0])));

	}

    result[idx] = scat;
}


