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
    
    m_tot = sqrt(mx**2 + my**2 + mz**2)
    m_tot_nz = (m_tot != 0.0)
    
    m_ip = sqrt(mx**2 + my**2) # in-plane projection of m unit vector
    m_ip_nz = (m_ip != 0.0) # mask for nonzero
    
    mx_ip = zeros_like(mx) 
    my_ip = zeros_like(my) 
    mx_ip[m_ip_nz] = mx[m_ip_nz] / m_ip[m_ip_nz]
    my_ip[m_ip_nz] = my[m_ip_nz] / m_ip[m_ip_nz]
    
    rhoM_ip = zeros_like(rhoM)
    rhoM_ip[m_tot_nz] = rhoM[m_tot_nz] * m_ip[m_tot_nz] / m_tot[m_tot_nz]
    
    KSQREL = kz**2 + PI4*rhoN[0] # fronting medium removed from effective kz
    
    S1 = sqrt(PI4*(rhoN + rhoM_ip)-KSQREL)
    S3 = sqrt(PI4*(rhoN - rhoM_ip)-KSQREL)
    mu = (complex(1.0)  + mx_ip + 1j * my_ip) / (complex(1.0)  + mx_ip - 1j * my_ip)
    
    #mu_1 = (complex(1.0)  + mx + 1j * mz - my) / (complex(1.0)  + mx - 1j * mz + my)
    #mu_3 = (complex(-1.0) + mx + 1j * mz - my) / (complex(-1.0) + mx - 1j * mz + my)
    
    l=0
    step = 1
    z = complex(0)
    
    for I in range(N-1):
        # chi in layer l
        lp = l+step
        chi_l = matrix(\
            [[1,              1,             1,              1            ],
             [mu[l],          mu[l],        -mu[l],         -mu[l]        ],
             [S1[l],         -S1[l],         S3[l],         -S3[l]        ],
             [mu[l]*S1[l],   -mu[l]*S1[l],  -mu[l]*S3[l],    mu[l]*S3[l]] ])
        
        # inverse chi at layer l+1
        chi_inv_lp = matrix(\
            [[ 1,  1/mu[lp],  1/S1[lp],  1/(mu[lp]*S1[lp]) ],
             [ 1,  1/mu[lp], -1/S1[lp], -1/(mu[lp]*S1[lp]) ],
             [ 1, -1/mu[lp],  1/S3[lp], -1/(mu[lp]*S3[lp]) ],
             [ 1, -1/mu[lp], -1/S3[lp],  1/(mu[lp]*S3[lp]) ]])
        chi_inv_lp *= 0.25
        
        S_l_vec = matrix([[exp( 1j * S1[l] * z),
                           exp(-1j * S1[l] * z),
                           exp( 1j * S3[l] * z),
                           exp(-1j * S3[l] * z)]])
                          
        S_lp_inv_vec = matrix([[exp(-1j * S1[lp] * z)],
                               [exp( 1j * S1[lp] * z)],
                               [exp(-1j * S3[lp] * z)],
                               [exp( 1j * S3[lp] * z)]])                         
        
        A = S_lp_inv_vec * chi_inv_lp * chi_l * S_l_vec
        
        newB = A * newB
        B[I] = newB.copy()
        z += dz[l]
    
    newB = unitary_LAB_SAM_LAB(newB, AGUIDE)
    
    
def get_U_sam_lab = function(AGUIDE) {
    C = complex(cos(AGUIDE/2.0*pi/180.))
    IS = 1j * sin(AGUIDE/2.0*pi/180.)
    U = matrix([ [C , IS, 0 , 0 ],
                 [IS, C , 0 , 0 ],
                 [0 , 0 , C , IS],
                 [0 , 0 , IS, C ] ])
    return U

def get_Uinv_sam_lab(AGUIDE):
    C = complex(cos(AGUIDE/2.0*pi/180.))
    NS = 1j * -sin(AGUIDE/2.0*pi/180.)
    Uinv = matrix([ [C , NS, 0 , 0 ],
                    [NS, C , 0 , 0 ],
                    [0 , 0 , C , NS],
                    [0 , 0 , NS, C ] ])
    return Uinv

def unitary_LAB_SAM_LAB(A, AGUIDE):
    """ perform rotation of coordinate system on one side of matrix
    and perform inverse on the other side """
    U = get_U_sam_lab(AGUIDE)
    Uinv = get_U_sam_lab(AGUIDE)
    CST =  (U * A) * Uinv
    return CST
