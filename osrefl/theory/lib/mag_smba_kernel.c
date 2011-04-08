CUDA_KERNEL

cudaBorn(int nx, int ny, int nz, int nqx, int nqy, int nqz,
		const Real density[], const Real magDensity[],
		const Real mx[],const Real my[],const Real mz[],
		const Real x[], const Real y[], const Real z[],
		const Real Qx[], const Real Qy[], const Real Qz[],
		const Cplx psi_in_one[],
		const Cplx psi_in_two[],
		const Cplx psi_out_one[],
		const Cplx psi_out_two[],
		const Real qx_refract[],
		int qxi,
		Cplx resultuu[],Cplx resultdd[],Cplx resultud[],Cplx resultdu[])
{

    const int idx = (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x + threadIdx.x;
    if (idx >= nqy*nqz) return;

    const int qzi = idx%nqz;
    const int qyi = (idx/nqz)%nqy;

    const Cplx I(0.0,1.0);

    Cplx ruu;
    Cplx rdd;
    Cplx rud;
    Cplx rdu;



    Cplx ft[] = {0+0*I,0+0*I,0+0*I,0+0*I};
    Cplx scat[]= {0+0*I,0+0*I,0+0*I,0+0*I};

    Cplx scat_PioPoo;
    Cplx scat_PioPot;
    Cplx scat_PitPoo;
    Cplx scat_PitPot;

    Real qMag;
    Real QdotM;
    Real qnx;
    Real qny;
    Real qnz;

    Real qcx;
    Real qcy;
    Real qcz;

    qMag = sqrt(pow(Qx[qxi],2)+pow(Qy[qyi],2)+pow(Qz[qzi],2));
    qnx = Qx[qxi]/qMag;
    qny = Qy[qyi]/qMag;
    qnz = Qz[qzi]/qMag;


    int densityidx = 0;

    for (int xi=0; xi < nx; xi++) {
        for (int yi=0; yi < ny; yi++) {
            for (int zi=0; zi < nz; zi++) {

            	QdotM = (qnx * mx[densityidx])+(qny * my[densityidx])+(qnz * mz[densityidx]);

            	qcx = mx[densityidx] - (qnx*QdotM);
            	qcy = my[densityidx] - (qny*QdotM);
            	qcz = mz[densityidx] - (qnz*QdotM);

            	ruu = density[densityidx] + qcx*magDensity[densityidx];
            	rdd = density[densityidx] + qcx*magDensity[densityidx];

            	rud = (qcy + I * qcz) * magDensity[densityidx];
            	rdu = (qcy - I * qcz) * magDensity[densityidx];

				ft[0] = ruu*exp(I*qx_refract[idx]*x[xi])*exp(I*Qy[qyi]*y[yi])*exp(I*Qz[qzi]*z[zi]);
				ft[1] = rdd*exp(I*qx_refract[idx]*x[xi])*exp(I*Qy[qyi]*y[yi])*exp(I*Qz[qzi]*z[zi]);
				ft[2] = rud*exp(I*qx_refract[idx]*x[xi])*exp(I*Qy[qyi]*y[yi])*exp(I*Qz[qzi]*z[zi]);
				ft[3] = rdu*exp(I*qx_refract[idx]*x[xi])*exp(I*Qy[qyi]*y[yi])*exp(I*Qz[qzi]*z[zi]);



				for (int cs=0; cs < 4;cs++){

					scat_PioPoo = psi_in_one[idx] * ft[cs] * psi_out_one[idx];
					scat_PioPot = psi_in_one[idx] * ft[cs] * psi_out_two[idx];
					scat_PitPoo = psi_in_two[idx] * ft[cs] * psi_out_one[idx];
					scat_PitPot = psi_in_two[idx] * ft[cs] * psi_out_two[idx];

					scat[cs] += (scat_PioPoo + scat_PioPot + scat_PitPoo + scat_PitPot);
				}

                densityidx++;
            }
       }
    }

    for (int cs=0; cs < 4;cs++)
    {
		if (Qx[qxi] == 0.0)
		{
			scat[cs] *= (x[1]-x[0]);

		}
		else
		{
			scat[cs] *= (-I/Qx[qxi])*(1-exp(I*Qx[qxi]*(x[1]-x[0])));

		}
		if (Qy[qyi] == 0.0)
		{
			scat[cs] *= (y[1]-y[0]);

		}
		else
		{
			scat[cs] *= (-I/Qy[qyi])*(1-exp(I*Qy[qyi]*(y[1]-y[0])));


		}
		if (Qz[qzi] == 0.0)
		{
			scat[cs] *= (z[1]-z[0]);

		}
		else
		{
			scat[cs] *= (-I/Qz[qzi])*(1-exp(I*Qz[qzi]*(z[1]-z[0])));

		}
    }

    resultuu[idx] = scat[0];
    resultdd[idx] = scat[1];
    resultud[idx] = scat[2];
    resultdu[idx] = scat[3];



}


