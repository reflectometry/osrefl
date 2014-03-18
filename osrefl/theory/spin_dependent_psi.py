from numpy import *

EPSILON = 1e-10

def calculateRB(kz, dz, rhoN, rhoM, mx, my, mz, A, B, C):
    """\
    Calculation of reflectivity in magnetic sample in framework that 
    also yields the wavefunction in each layer for DWBA.  
    Arguments:
    
     kz is the incoming momentum along the z direction
     dz is the thickness of each layer (top and bottom thickness ignored)
     rhoN is nuclear scattering length density array
     rhoM is magnetic scattering length density array
     mx, my, mz are components of unit vector along M for layer
    ###################################################################
    #  all of dz, rhoN, rhoM, mx, my, mz should have length equal to  #
    #  the number of layers in the sample including fronting and      #
    #  substrate.                                                     #
    ###################################################################
     A, B, C are components of unit vector describing quantization direction:
       A is along H_external, B and C are perpendicular to A and each other and
       B x C = A
    """
    # sld is array([[sld, thickness, mu], [...], ...])
    # ordered from top (vacuum usually) to bottom (substrate)
    N = len(rhoN) 
    PI4 = complex(pi * 4.0)
    
    if abs(kz) < EPSILON:
        # we're at zero:
        YA = complex(-1.0);
        YB = complex(0.0);
        YC = complex(0.0);
        YD = complex(-1.0);
        return [YA, YB, YC, YD];
    
    
    B = zeros((N,4,4)) # one matrix per layer
    
    newB = array([[1, 0, 0, 0], 
                  [0, 1, 0, 0], 
                  [0, 0, 1, 0], 
                  [0, 0, 0, 1]], dtype='complex')
                  
    KSQREL = kz**2 + PI4*rhoN[0] # fronting medium removed from effective kz
    
    S1 = sqrt(PI4*(rhoN + rhoM)-KSQREL)
    S3 = sqrt(PI4*(rhoN - rhoM)-KSQREL)
    mu_1 = (complex(1.0)  + mx + 1j * mz - my) / (complex(1.0)  + mx - 1j * mz + my)
    mu_3 = (complex(-1.0) + mx + 1j * mz - my) / (complex(-1.0) + mx - 1j * mz + my)
    
    l=0
    step = 1
    
    for I in range(N-1):
        # chi in layer l
        lp = l+step
        chi_l = matrix(\
            [[1,              1,             1,              1            ],
             [mu_1[l],        mu_1[l],       mu_3[l],        mu_3[l]      ],
             [S1[l],         -S1[l],         S3[l],         -S3[l]        ],
             [mu_1[l]*S1[l], -mu_1[l]*S1[l], mu_3[l]*S3[l], -mu_3[l]*S3[l]]])
        
        # inverse chi at layer l+1
        chi_inv_lp = matrix(\
            [[ mu_3[lp], -1,  mu_3[lp]/S1[lp], -1/S1[lp] ],
             [ mu_3[lp], -1, -mu_3[lp]/S1[lp],  1/S1[lp] ],
             [-mu_1[lp],  1, -mu_1[lp]/S3[lp],  1/S3[lp] ],
             [-mu_1[lp],  1,  mu_1[lp]/S3[lp], -1/S3[lp] ]])
        chi_inv_lp *= 1/(2.0*(mu_3[lp] - mu_1[lp]))
        
        S_l_vec = matrix([[exp( 1j * S1[l] * z),
                           exp(-1j * S1[l] * z),
                           exp( 1j * S3[l] * z),
                           exp(-1j * S3[l] * z)]])
                          
        S_lp_inv_vec = matrix([[exp(-1j * S1[lp] * z)],
                               [exp( 1j * S1[lp] * z)],
                               [exp(-1j * S3[lp] * z)],
                               [exp( 1j * S3[lp] * z)]])                         
        
        l += step
        sld_NL = sld_n[L][0];
        sld_NLi = sld_n[L][2]
        sld_ML = sld_m[L][0];
        sld_NLn = sld_n[Ln];
        sld_MLn = sld_m[Ln];
        mu_L = mu_1[L]; # need to fill this above!
        mu_Ln = mu_1[Ln];
        inv_mu_Ln = 1.0/mu_Ln;       
        mu_ratio = mu_L/inv_mu_Ln;

        S1L = sqrt((PI4*(sld_NL + sld_ML)-KSQREL + 1j*PI4*sld_L.sldi));
        muS1L = Cplx.multiply(mu_L, S1L);
        S3L = Cplx.sqrt(new Cplx(PI4*(sld_L.sld - sld_L.sldm)-KSQREL,  PI4*sld_L.sldi));
        muS3L = Cplx.multiply(mu_L, S3L);
        
        S1Ln = (Cplx.sqrt(new Cplx(PI4*(sld_Ln.sld + sld_Ln.sldm)-KSQREL,  PI4*sld_Ln.sldi)));
        inv_S1Ln = S1Ln.inverse();
        inv_mu_S1Ln = Cplx.multiply(inv_mu_Ln, inv_S1Ln);
        S3Ln = (Cplx.sqrt(new Cplx(PI4*(sld_Ln.sld - sld_Ln.sldm)-KSQREL,  PI4*sld_Ln.sldi)));
        inv_S3Ln = S3Ln.inverse();
        inv_mu_S3Ln = Cplx.multiply(inv_mu_Ln, inv_S3Ln);
        
