// This program is public domain.

#include <iostream>
#include <fstream>
#include <complex>
#include <cstdlib>
#include "reflcalc.h"

extern "C" void
Cr4xb(const int &N, const double D[], const double SIGMA[],
      const int &IP,
      const double RHO[], const double IRHO[],
      const double RHOM[], const Cplx U1[], const Cplx U3[],
      const double &AGUIDE, const double &KZ,
      Cplx &YA, Cplx &YB, Cplx &YC, Cplx &YD)

{
/*
C A,B,C,D is the result for the ++, +-, -+ and -- cross sections
C
C Notes:
C
C 1. If Q is negative then the beam is assumed to come in from the
C bottom of the sample, and all the layers are reversed.
C
C 2. The fronting and backing materials are assumed to be semi-infinite,
C so depth is ignored for the first and last layer.
C
C 3. Absorption is ignored for the fronting material, or the backing
C material for negative Q.  For negative Q, the beam is coming in
C through the side of the substrate, and you will need to multiply
C by a substrate absorption factor depending on the path length through
C the substrate.  For neutron reflectivity, this is approximately
C constant for the angles we need to consider.
C
C 4. This subroutine does not deal with any component of sample moment
C that may lie out of the plane of the film.  Such a perpendicular
C component will cause a neutron precession, therefore an additional
C spin flip term.  If reflectivity data from a sample with an
C out-of-plane moment is modeled using this subroutine, one will
C obtain erroneous results, since all of the spin flip scattering
C will be attributed to in-plane moments perpendicular to the neutron.

C Created Jan 13, 2015 by Brian B. Maranville
*/

//     paramters
      int I,L,LP,STEP;

//    variables calculating S1, S3, and exponents
      double E0;
      Cplx S1L,S3L,S1LP,S3LP,ES1L,ES3L,ENS1L,ENS3L,ES1LP,ES3LP,ENS1LP,ENS3LP;
      Cplx FS1S1, FS3S1, FS1S3, FS3S3;

//    completely unrolled matrices for B=A*B update
      Cplx DELTA,U1L,U3L,U1LP,U3LP;
      Cplx DU1U1, DU1U3, DU3U1, DU3U3;
      Cplx Z;
      Cplx A11,A12,A13,A14,A21,A22,A23,A24;
      Cplx A31,A32,A33,A34,A41,A42,A43,A44;
      Cplx B11,B12,B13,B14,B21,B22,B23,B24;
      Cplx B31,B32,B33,B34,B41,B42,B43,B44;
      Cplx C1,C2,C3,C4;

//    variables for translating resulting B into a signal
      Cplx DETW;

//    constants
      const Cplx CR(1.0,0.0);
      const Cplx CI(0.0,1.0);
      const double PI4=12.566370614359172e-6;
//    Check for KZ near zero.  If KZ < 0, reverse the indices
      if (KZ<=-1.e-10) {
         L=N-1;
         STEP=-1;
      } else if (KZ>=1.e-10) {
         L=0;
         STEP=1;
      } else {
         YA = -1.;
         YB = 0.;
         YC = 0.;
         YD = -1.;
         return;
      }

//     B = I
      B11=Cplx(1.0,0.0);
      B12=Cplx(0.0,0.0);
      B13=Cplx(0.0,0.0);
      B14=Cplx(0.0,0.0);
      B21=Cplx(0.0,0.0);
      B22=Cplx(1.0,0.0);
      B23=Cplx(0.0,0.0);
      B24=Cplx(0.0,0.0);
      B31=Cplx(0.0,0.0);
      B32=Cplx(0.0,0.0);
      B33=Cplx(1.0,0.0);
      B34=Cplx(0.0,0.0);
      B41=Cplx(0.0,0.0);
      B42=Cplx(0.0,0.0);
      B43=Cplx(0.0,0.0);
      B44=Cplx(1.0,0.0);

//    Changing the target KZ is equivalent to subtracting the fronting
//    medium SLD.
      if (IP > 0) {
        // IP = 1 specifies polarization of the incident beam I+
        E0 = KZ*KZ + PI4*(RHO[L]+RHOM[L]);
      } else {
        // IP = 0 specifies polarization of the incident beam I-
        E0 = KZ*KZ + PI4*(RHO[L]-RHOM[L]);
      }
      
      Z = 0.0;
      
      if (N>1) {
        // chi in layer 1
        LP = L + STEP;
        S1L = -sqrt(CR*PI4*(RHO[L]+RHOM[L])-E0 + CI*PI4*IRHO[L]);
        S3L = -sqrt(CR*PI4*(RHO[L]-RHOM[L])-E0 + CI*PI4*IRHO[L]);
        S1LP = -sqrt(CR*PI4*(RHO[LP]+RHOM[LP])-E0 + CI*PI4*IRHO[LP]);
        S3LP = -sqrt(CR*PI4*(RHO[LP]-RHOM[LP])-E0 + CI*PI4*IRHO[LP]);
        U1L = U1[L];
        U3L = U3[L];
        U1LP = U1[LP];
        U3LP = U3[LP];
        
        DELTA = 0.5*CR / (U3LP - U1LP);
        
        FS1S1 = S1L/S1LP;
        FS1S3 = S1L/S3LP;
        FS3S1 = S3L/S1LP;
        FS3S3 = S3L/S3LP;
         
        B11 = DELTA *  U3LP * (1.0 + FS1S1);
        B12 = DELTA *  U3LP * (1.0 - FS1S1);
        B13 = DELTA *  -1.0 * (1.0 + FS3S1);
        B14 = DELTA *  -1.0 * (1.0 - FS3S1);
        
        B21 = DELTA *  U3LP * (1.0 - FS1S1);
        B22 = DELTA *  U3LP * (1.0 + FS1S1);
        B23 = DELTA *  -1.0 * (1.0 - FS3S1);
        B24 = DELTA *  -1.0 * (1.0 + FS3S1);
        
        B31 = DELTA * -U1LP * (1.0 + FS1S3);
        B32 = DELTA * -U1LP * (1.0 - FS1S3);
        B33 = DELTA *   1.0 * (1.0 + FS3S3);
        B34 = DELTA *   1.0 * (1.0 - FS3S3);
        
        B41 = DELTA * -U1LP * (1.0 - FS1S3);
        B42 = DELTA * -U1LP * (1.0 + FS1S3);
        B43 = DELTA *   1.0 * (1.0 - FS3S3);
        B44 = DELTA *   1.0 * (1.0 + FS3S3);
        
        Z += D[LP];
        L = LP;
      }
      
//    Process the loop once for each interior layer, either from
//    front to back or back to front.
      for (I=1; I < N-1; I++) {
        LP = L + STEP;
        S1L = S1LP; // copy from the layer before
        S3L = S3LP; //
        S1LP = -CR*sqrt(CR*PI4*(RHO[LP]+RHOM[LP])-E0 + CI*PI4*IRHO[LP]);
        S3LP = -CR*sqrt(CR*PI4*(RHO[LP]-RHOM[LP])-E0 + CI*PI4*IRHO[LP]);
        U1L = U1[L];
        U3L = U3[L];
        U1LP = U1[LP];
        U3LP = U3[LP];
        
        DELTA = 0.5*CR / (U3LP - U1LP);
        DU1U1 = -U1LP + U1L;
        DU1U3 = -U1LP + U3L;
        DU3U1 = U3LP - U1L;
        DU3U3 = U3LP - U3L;        

        ES1L = exp(S1L*Z);
        ENS1L = CR / ES1L;
        ES1LP = exp(S1LP*Z);
        ENS1LP = CR / ES1LP;
        ES3L = exp(S3L*Z);
        ENS3L = CR / ES3L;
        ES3LP = exp(S3LP*Z);
        ENS3LP = CR / ES3LP;
        
        FS1S1 = S1L/S1LP;
        FS1S3 = S1L/S3LP;
        FS3S1 = S3L/S1LP;
        FS3S3 = S3L/S3LP;
        
        A11 = DELTA * ES1L  * ENS1LP * DU3U1 * (1.0 + FS1S1);
        A12 = DELTA * ENS1L * ENS1LP * DU3U1 * (1.0 - FS1S1);
        A13 = DELTA * ES3L  * ENS1LP * DU3U3 * (1.0 + FS3S1);
        A14 = DELTA * ENS3L * ENS1LP * DU3U3 * (1.0 - FS3S1);
        
        A21 = DELTA * ES1L  * ES1LP  * DU3U1 * (1.0 - FS1S1);
        A22 = DELTA * ENS1L * ES1LP  * DU3U1 * (1.0 + FS1S1);
        A23 = DELTA * ES3L  * ES1LP  * DU3U3 * (1.0 - FS3S1);
        A24 = DELTA * ENS3L * ES1LP  * DU3U3 * (1.0 + FS3S1);
        
        A31 = DELTA * ES1L  * ENS3LP * DU1U1 * (1.0 + FS1S3);
        A32 = DELTA * ENS1L * ENS3LP * DU1U1 * (1.0 - FS1S3);
        A33 = DELTA * ES3L  * ENS3LP * DU1U3 * (1.0 + FS3S3);
        A34 = DELTA * ENS3L * ENS3LP * DU1U3 * (1.0 - FS3S3);
        
        A41 = DELTA * ES1L  * ES3LP * DU1U1 * (1.0 - FS1S3);
        A42 = DELTA * ENS1L * ES3LP * DU1U1 * (1.0 + FS1S3);
        A43 = DELTA * ES3L  * ES3LP * DU1U3 * (1.0 - FS3S3);
        A44 = DELTA * ENS3L * ES3LP * DU1U3 * (1.0 + FS3S3);


//    Matrix update B=A*B
        C1=A11*B11+A12*B21+A13*B31+A14*B41;
        C2=A21*B11+A22*B21+A23*B31+A24*B41;
        C3=A31*B11+A32*B21+A33*B31+A34*B41;
        C4=A41*B11+A42*B21+A43*B31+A44*B41;
        B11=C1;
        B21=C2;
        B31=C3;
        B41=C4;

        C1=A11*B12+A12*B22+A13*B32+A14*B42;
        C2=A21*B12+A22*B22+A23*B32+A24*B42;
        C3=A31*B12+A32*B22+A33*B32+A34*B42;
        C4=A41*B12+A42*B22+A43*B32+A44*B42;
        B12=C1;
        B22=C2;
        B32=C3;
        B42=C4;

        C1=A11*B13+A12*B23+A13*B33+A14*B43;
        C2=A21*B13+A22*B23+A23*B33+A24*B43;
        C3=A31*B13+A32*B23+A33*B33+A34*B43;
        C4=A41*B13+A42*B23+A43*B33+A44*B43;
        B13=C1;
        B23=C2;
        B33=C3;
        B43=C4;

        C1=A11*B14+A12*B24+A13*B34+A14*B44;
        C2=A21*B14+A22*B24+A23*B34+A24*B44;
        C3=A31*B14+A32*B24+A33*B34+A34*B44;
        C4=A41*B14+A42*B24+A43*B34+A44*B44;
        B14=C1;
        B24=C2;
        B34=C3;
        B44=C4;
        
        Z += D[LP];
        L = LP;
      }
//    Done computing B = A(N)*...*A(2)*A(1)*I

      DETW=(B44*B22-B24*B42);

//    Calculate reflectivity coefficients specified by POLSTAT
      YA = (B24*B41-B21*B44)/DETW;
      YB = (B21*B42-B41*B22)/DETW;
      YC = (B24*B43-B23*B44)/DETW;
      YD = (B23*B42-B43*B22)/DETW;

}


extern "C" void
Cc4xb(const int &N, const double D[], const double SIGMA[],
      const int &IP,
      const double RHO[], const double IRHO[],
      const double RHOM[], const Cplx U1[], const Cplx U3[],
      const double &AGUIDE, const double &KZ,
      Cplx &C1, Cplx &C2, Cplx &C3, Cplx &C4)

{
/*
C Calculate the C-coefficients of the wavefunction amplitude in 
C each layer.  The B-matrix is passed in empty and filled.
C Notes:
C
C 1. If Q is negative then the beam is assumed to come in from the
C bottom of the sample, and all the layers are reversed.
C
C 2. The fronting and backing materials are assumed to be semi-infinite,
C so depth is ignored for the first and last layer.
C
C 3. Absorption is ignored for the fronting material, or the backing
C material for negative Q.  For negative Q, the beam is coming in
C through the side of the substrate, and you will need to multiply
C by a substrate absorption factor depending on the path length through
C the substrate.  For neutron reflectivity, this is approximately
C constant for the angles we need to consider.
C
C 4. Magnetic scattering is ignored for the fronting and backing.
C
C 5. This subroutine does not deal with any component of sample moment
C that may lie out of the plane of the film.  Such a perpendicular
C component will cause a neutron presession, therefore an additional
C spin flip term.  If reflectivity data from a sample with an
C out-of-plane moment is modeled using this subroutine, one will
C obtain erroneous results, since all of the spin flip scattering
C will be attributed to in-plane moments perpendicular to the neutron.


C $Log$
C Created Jan 13, 2015 by Brian Maranville
C
*/

//     paramters
      int I,L,LP,STEP;

//    variables calculating S1, S3, and exponents
      double E0;
      Cplx S1L,S3L,S1LP,S3LP,ES1L,ES3L,ENS1L,ENS3L,ES1LP,ES3LP,ENS1LP,ENS3LP;
      Cplx FS1S1, FS3S1, FS1S3, FS3S3;

//    completely unrolled matrices for B=A*B update
      Cplx DELTA,U1L,U3L,U1LP,U3LP;
      Cplx DU1U1, DU1U3, DU3U1, DU3U3;
      Cplx Z;
      Cplx A11,A12,A13,A14,A21,A22,A23,A24;
      Cplx A31,A32,A33,A34,A41,A42,A43,A44;
      Cplx B11,B12,B13,B14,B21,B22,B23,B24;
      Cplx B31,B32,B33,B34,B41,B42,B43,B44;
      Cplx D1,D2,D3,D4;
      Cplx YA, YB, YC, YD;
      
      int B_offset = 0;
      Cplx * B = new Cplx[N*16];
      // the size of B is N * 4 * 4, and the order is 
      // B[layer n][row i][col j] = B[n*4*4 + i*4 + j]

//    variables for translating resulting B into a signal
      Cplx DETW;

//    constants
      const Cplx CR(1.0,0.0);
      const Cplx CI(0.0,1.0);
      const double PI4=12.566370614359172e-6;
//    Check for KZ near zero.  If KZ < 0, reverse the indices
      if (KZ<=-1.e-10) {
         L=N-1;
         STEP=-1;
      } else if (KZ>=1.e-10) {
         L=0;
         STEP=1;
      } else {
         YA = -1.;
         YB = 0.;
         YC = 0.;
         YD = -1.;
         return;
      }

//     B = I
      B11=Cplx(1.0,0.0);
      B12=Cplx(0.0,0.0);
      B13=Cplx(0.0,0.0);
      B14=Cplx(0.0,0.0);
      B21=Cplx(0.0,0.0);
      B22=Cplx(1.0,0.0);
      B23=Cplx(0.0,0.0);
      B24=Cplx(0.0,0.0);
      B31=Cplx(0.0,0.0);
      B32=Cplx(0.0,0.0);
      B33=Cplx(1.0,0.0);
      B34=Cplx(0.0,0.0);
      B41=Cplx(0.0,0.0);
      B42=Cplx(0.0,0.0);
      B43=Cplx(0.0,0.0);
      B44=Cplx(1.0,0.0);

//    Changing the target KZ is equivalent to subtracting the fronting
//    medium SLD.
      if (IP > 0) {
        // IP = 1 specifies polarization of the incident beam I+
        E0 = KZ*KZ + PI4*(RHO[L]+RHOM[L]);
      } else {
        // IP = 0 specifies polarization of the incident beam I-
        E0 = KZ*KZ + PI4*(RHO[L]-RHOM[L]);
      }
      
      Z = 0.0;
      
      if (N>1) {
        // chi in layer 1
        LP = L + STEP;
        S1L = -sqrt(CR*PI4*(RHO[L]+RHOM[L])-E0 + CI*PI4*IRHO[L]);
        S3L = -sqrt(CR*PI4*(RHO[L]-RHOM[L])-E0 + CI*PI4*IRHO[L]);
        S1LP = -sqrt(CR*PI4*(RHO[LP]+RHOM[LP])-E0 + CI*PI4*IRHO[LP]);
        S3LP = -sqrt(CR*PI4*(RHO[LP]-RHOM[LP])-E0 + CI*PI4*IRHO[LP]);
        U1L = U1[L];
        U3L = U3[L];
        U1LP = U1[LP];
        U3LP = U3[LP];
        
        DELTA = 0.5*CR / (U3LP - U1LP);
        
        FS1S1 = S1L/S1LP;
        FS1S3 = S1L/S3LP;
        FS3S1 = S3L/S1LP;
        FS3S3 = S3L/S3LP;
         
        B11 = DELTA *  U3LP * (1.0 + FS1S1); B[B_offset++] = B11;
        B12 = DELTA *  U3LP * (1.0 - FS1S1); B[B_offset++] = B12;
        B13 = DELTA *  -1.0 * (1.0 + FS3S1); B[B_offset++] = B13;
        B14 = DELTA *  -1.0 * (1.0 - FS3S1); B[B_offset++] = B14;
        
        B21 = DELTA *  U3LP * (1.0 - FS1S1); B[B_offset++] = B21;
        B22 = DELTA *  U3LP * (1.0 + FS1S1); B[B_offset++] = B22;
        B23 = DELTA *  -1.0 * (1.0 - FS3S1); B[B_offset++] = B23;
        B24 = DELTA *  -1.0 * (1.0 + FS3S1); B[B_offset++] = B24;
        
        B31 = DELTA * -U1LP * (1.0 + FS1S3); B[B_offset++] = B31;
        B32 = DELTA * -U1LP * (1.0 - FS1S3); B[B_offset++] = B32;
        B33 = DELTA *   1.0 * (1.0 + FS3S3); B[B_offset++] = B33;
        B34 = DELTA *   1.0 * (1.0 - FS3S3); B[B_offset++] = B34;
        
        B41 = DELTA * -U1LP * (1.0 - FS1S3); B[B_offset++] = B41;
        B42 = DELTA * -U1LP * (1.0 + FS1S3); B[B_offset++] = B42;
        B43 = DELTA *   1.0 * (1.0 - FS3S3); B[B_offset++] = B43;
        B44 = DELTA *   1.0 * (1.0 + FS3S3); B[B_offset++] = B44;
        
        Z += D[LP];
        L = LP;
      }
      
//    Process the loop once for each interior layer, either from
//    front to back or back to front.
      for (I=1; I < N-1; I++) {
        LP = L + STEP;
        S1L = S1LP; // copy from the layer before
        S3L = S3LP; //
        S1LP = -CR*sqrt(CR*PI4*(RHO[LP]+RHOM[LP])-E0 + CI*PI4*IRHO[LP]);
        S3LP = -CR*sqrt(CR*PI4*(RHO[LP]-RHOM[LP])-E0 + CI*PI4*IRHO[LP]);
        U1L = U1[L];
        U3L = U3[L];
        U1LP = U1[LP];
        U3LP = U3[LP];
        
        DELTA = 0.5*CR / (U3LP - U1LP);
        DU1U1 = -U1LP + U1L;
        DU1U3 = -U1LP + U3L;
        DU3U1 = U3LP - U1L;
        DU3U3 = U3LP - U3L;        

        ES1L = exp(S1L*Z);
        ENS1L = CR / ES1L;
        ES1LP = exp(S1LP*Z);
        ENS1LP = CR / ES1LP;
        ES3L = exp(S3L*Z);
        ENS3L = CR / ES3L;
        ES3LP = exp(S3LP*Z);
        ENS3LP = CR / ES3LP;
        
        FS1S1 = S1L/S1LP;
        FS1S3 = S1L/S3LP;
        FS3S1 = S3L/S1LP;
        FS3S3 = S3L/S3LP;
        
        A11 = DELTA * ES1L  * ENS1LP * DU3U1 * (1.0 + FS1S1);
        A12 = DELTA * ENS1L * ENS1LP * DU3U1 * (1.0 - FS1S1);
        A13 = DELTA * ES3L  * ENS1LP * DU3U3 * (1.0 + FS3S1);
        A14 = DELTA * ENS3L * ENS1LP * DU3U3 * (1.0 - FS3S1);
        
        A21 = DELTA * ES1L  * ES1LP  * DU3U1 * (1.0 - FS1S1);
        A22 = DELTA * ENS1L * ES1LP  * DU3U1 * (1.0 + FS1S1);
        A23 = DELTA * ES3L  * ES1LP  * DU3U3 * (1.0 - FS3S1);
        A24 = DELTA * ENS3L * ES1LP  * DU3U3 * (1.0 + FS3S1);
        
        A31 = DELTA * ES1L  * ENS3LP * DU1U1 * (1.0 + FS1S3);
        A32 = DELTA * ENS1L * ENS3LP * DU1U1 * (1.0 - FS1S3);
        A33 = DELTA * ES3L  * ENS3LP * DU1U3 * (1.0 + FS3S3);
        A34 = DELTA * ENS3L * ENS3LP * DU1U3 * (1.0 - FS3S3);
        
        A41 = DELTA * ES1L  * ES3LP * DU1U1 * (1.0 - FS1S3);
        A42 = DELTA * ENS1L * ES3LP * DU1U1 * (1.0 + FS1S3);
        A43 = DELTA * ES3L  * ES3LP * DU1U3 * (1.0 - FS3S3);
        A44 = DELTA * ENS3L * ES3LP * DU1U3 * (1.0 + FS3S3);


//    Matrix update B=A*B
        D1=A11*B11+A12*B21+A13*B31+A14*B41;
        D2=A21*B11+A22*B21+A23*B31+A24*B41;
        D3=A31*B11+A32*B21+A33*B31+A34*B41;
        D4=A41*B11+A42*B21+A43*B31+A44*B41;
        B11=D1;
        B21=D2;
        B31=D3;
        B41=D4;

        D1=A11*B12+A12*B22+A13*B32+A14*B42;
        D2=A21*B12+A22*B22+A23*B32+A24*B42;
        D3=A31*B12+A32*B22+A33*B32+A34*B42;
        D4=A41*B12+A42*B22+A43*B32+A44*B42;
        B12=D1;
        B22=D2;
        B32=D3;
        B42=D4;

        D1=A11*B13+A12*B23+A13*B33+A14*B43;
        D2=A21*B13+A22*B23+A23*B33+A24*B43;
        D3=A31*B13+A32*B23+A33*B33+A34*B43;
        D4=A41*B13+A42*B23+A43*B33+A44*B43;
        B13=D1;
        B23=D2;
        B33=D3;
        B43=D4;

        D1=A11*B14+A12*B24+A13*B34+A14*B44;
        D2=A21*B14+A22*B24+A23*B34+A24*B44;
        D3=A31*B14+A32*B24+A33*B34+A34*B44;
        D4=A41*B14+A42*B24+A43*B34+A44*B44;
        B14=D1;
        B24=D2;
        B34=D3;
        B44=D4;
        
        Z += D[LP];
        L = LP;
      }
//    Done computing B = A(N)*...*A(2)*A(1)*I

      DETW=(B44*B22-B24*B42);

//    Calculate reflectivity coefficients specified by POLSTAT
      YA = (B24*B41-B21*B44)/DETW; // r++
      YB = (B21*B42-B41*B22)/DETW; // r+-
      YC = (B24*B43-B23*B44)/DETW; // r-+
      YD = (B23*B42-B43*B22)/DETW; // r--
      
//    Now get the C values per layer:
      delete[] B;

}

extern "C" void
magnetic_amplitude(const int layers,
                      const double d[], const double sigma[],
                      const double rho[], const double irho[],
                      const double rhoM[], const Cplx u1[], const Cplx u3[],
                      const double Aguide,
                      const int points, const double KZ[], const int rho_index[],
                      Cplx Ra[], Cplx Rb[], Cplx Rc[], Cplx Rd[])
{
  Cplx dummy1,dummy2;
  int ip;
  if (rhoM[0] == 0.0 && rhoM[layers-1] == 0.0) {
    ip = 1; // calculations for I+ and I- are the same in the fronting and backing.
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int i=0; i < points; i++) {
      const int offset = layers*(rho_index != NULL?rho_index[i]:0);
      Cr4xb(layers,d,sigma,ip,rho+offset,irho+offset,rhoM,u1,u3,
            Aguide,KZ[i],Ra[i],Rb[i],Rc[i],Rd[i]);
    }
  } else {
    ip = 1; // plus polarization
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int i=0; i < points; i++) {
      const int offset = layers*(rho_index != NULL?rho_index[i]:0);
      Cr4xb(layers,d,sigma,ip,rho+offset,irho+offset,rhoM,u1,u3,
            Aguide,KZ[i],Ra[i],Rb[i],dummy1,dummy2);
    }
    ip = 0; // minus polarization
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int i=0; i < points; i++) {
      const int offset = layers*(rho_index != NULL?rho_index[i]:0);
      Cr4xb(layers,d,sigma,ip,rho+offset,irho+offset,rhoM,u1,u3,
            Aguide,KZ[i],dummy1,dummy2,Rc[i],Rd[i]);
    }
  }
}

extern "C" void
magnetic_layer_amplitude(const int layers, const int ip,
                      const double d[], const double sigma[],
                      const double rho[], const double irho[],
                      const double rhoM[], const Cplx u1[], const Cplx u3[],
                      const double Aguide,
                      const int points, const double KZ[], const int rho_index[],
                      Cplx C1[], Cplx C2[], Cplx C3[], Cplx C4[])
{
  // returns all 4 C values for each layer.  r+ is C2 for layer 0, and r- is C4.
  // Also, t+ is C1 for layer N and t- is C3.  
  // C1 for layer 0 is 1.0 for I+, and C3 for layer 0 is I-
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for (int i=0; i < points; i++) {
    const int offset = layers*(rho_index != NULL?rho_index[i]:0);
    Cc4xb(layers,d,sigma,ip,rho+offset,irho+offset,rhoM,u1,u3,
          Aguide,KZ[i],C1[i],C2[i],C3[i],C4[i]);
  }
}


// $Id: magnetic_C.cc 2015-01-13 bbm $
